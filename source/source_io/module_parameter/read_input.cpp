#include "read_input.h"


#include <algorithm>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <array>
#include <vector>
#include <cassert>
#include <cctype>
#include <limits>
#include "source_base/formatter.h"
#include "source_base/global_file.h"
#include "source_base/global_function.h"
#include "source_base/tool_quit.h"
#include "source_base/tool_title.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <cerrno>

namespace ModuleIO
{

std::string longstring(const std::vector<std::string>& words)
{
    return FmtCore::join(" ", words);
}

bool assume_as_boolean(const std::string& val)
{
    const std::string val_ = FmtCore::lower(val);

    const std::array<std::string, 7> t_ = {"true", "1", "t", "yes", "y", "on", ".true."};
    const std::array<std::string, 7> f_ = {"false", "0", "f", "no", "n", "off", ".false."};
    // This will work because std::array<T, N>::size() is a constexpr function
    // Ouch it is of C++17 standard...
    // static_assert(t_.size() == f_.size(), "t_ and f_ must have the same lengths");
#ifdef __DEBUG // C++11 can do this
    assert(t_.size() == f_.size());
#endif

    if (std::find(t_.begin(), t_.end(), val_) != t_.end())
    {
        return true;
    }
    else if (std::find(f_.begin(), f_.end(), val_) != f_.end())
    {
        return false;
    }
    else
    {
        std::string warnmsg = "Bad boolean parameter ";
        warnmsg.append(val);
        warnmsg.append(", please check the input parameters in file INPUT");
        ModuleBase::WARNING_QUIT("Input", warnmsg);
    }
}

std::string to_dir(const std::string& str)
{
    std::string str_dir = str;
    if (str_dir.empty())
    {
        return "./";
    }
    else if (str_dir.back() != '/')
    {
        str_dir += "/";
    }
    
    return str_dir;
}

void read_information(std::stringstream& ifs, std::vector<std::string>& output, const std::string& delimiters)
{
    std::string line;
    getline(ifs, line);

    std::istringstream iss(line);
    std::string word;
    while (iss >> word)
    {
        if (delimiters.find(word[0]) != std::string::npos)
        {
            break;
        }
        output.push_back(word);
    }
}

bool ReadInput::check_mode = false;

bool filter_nonascii_and_comment(std::ifstream& ifs,
                       std::stringstream& out_ascii_stream)
{
	if (!ifs.is_open()) 
	{
		if (!ifs) return false;
    }

    std::streampos old_pos = ifs.tellg();
    ifs.clear();
    ifs.seekg(0, std::ios::beg);

	char c = '\0';
	while (ifs.get(c)) 
	{
		// If comment start, skip until end of line (but keep the newline)
		if (c == '#') 
		{
			char d = '\0';
			bool newline_found = false;
			while (ifs.get(d)) 
			{
				if (d == '\n' || d == '\r') 
				{
					// preserve line break in output
					out_ascii_stream.put('\n');
					// If CRLF, consume the LF after CR (already wrote a single '\n')
					if (d == '\r' && ifs.peek() == '\n') 
					{
						ifs.get(d); // consume '\n'
					}
					newline_found = true;
					break;
				}
            }
			if (!newline_found) 
			{
                // reached EOF while skipping comment
                break;
            }
            continue;
        }

        unsigned char uc = static_cast<unsigned char>(c);
		if (uc <= 0x7F) 
		{
			// ASCII character
            out_ascii_stream.put(c);
        }
		else 
		{
			// replace non-ASCII with space character
			out_ascii_stream.put(' ');
        }
    }

    // recover ifstream state and position
    ifs.clear();
    ifs.seekg(old_pos, std::ios::beg);

    return true;
}


ReadInput::ReadInput(const int& rank)
{
    this->rank = rank;

    this->item_h0();
}

void ReadInput::read_parameters(Parameter& param, const std::string& filename_in)
{
    ModuleBase::TITLE("ReadInput", "read_parameters");

    // 1. only rank 0 read the input file
    if (this->rank == 0)
    {
        // We can also easily add other input file formats here
        this->read_txt_input(param, filename_in);
    }

    // 2. check the number of atom types from STRU file
    // set the global directories
    this->set_global_dir(param.inp, param.sys); 
    if (this->check_ntype_flag && this->rank == 0)
    {
        check_ntype(param.globalv.global_in_stru, param.input.ntype);
    }

    // 3. broadcast input parameters
    // It must be after the check_ntype, because some parameters need to be filled due to ntype
    for (auto& bcastfunc: this->bcastfuncs)
    {
        bcastfunc(param);
    }

    // 4. set the globalv parameters, some parameters in different processes are different. e.g. rank, log_file
    this->set_globalv(param.inp, param.sys);

    // 5. check the value of the parameters
    // It must be after the check_ntype, because some parameters need to be checked according to ntype
    // It must be after the set_globalv, because some parameters need to be checked according to param.sys
    if (this->rank == 0)
    {
        for (auto& input_item: this->input_lists)
        {
            Input_Item* checkvalue_item = &(input_item.second);
            if (checkvalue_item->check_value != nullptr)
            {
                checkvalue_item->check_value(*checkvalue_item, param);
            }
        }
    }

    if (this->check_mode)
    {
        std::cout << "----------------------------------------------------------" << std::endl;
        std::cout << "  INPUT parameters have been successfully checked!" << std::endl;
        std::cout << "----------------------------------------------------------" << std::endl;
        exit(0);
        return;
    }
}

void ReadInput::create_directory(const Parameter& param)
{
    ModuleBase::TITLE("ReadInput", "create_directory");

    ModuleBase::Global_File::make_h0_output_dir(this->rank,
                                                param.globalv.global_out_dir,
                                                param.globalv.log_file);
    struct stat st;
    if (stat(PARAM.inp.read_file_dir.c_str(), &st) != 0 || !S_ISDIR(st.st_mode))
    {
        ModuleBase::WARNING_QUIT("ReadInput", "please set right files directory for reading in.");
    }
    return;
}

void ReadInput::write_parameters(const Parameter& param, const std::string& filename_out)
{
    if (this->rank == 0)
    {
        this->write_txt_input(param, filename_out);
    }
}

void ReadInput::read_txt_input(Parameter& param, const std::string& filename)
{
    ModuleBase::TITLE("ReadInput", "read_txt_input");

    std::stringstream ascii_stream;

	std::ifstream ifs(filename.c_str(), std::ios::in);

	if (!ifs)
	{
		std::cout << " Can't find the INPUT file." << std::endl;
		ModuleBase::WARNING_QUIT("Input::Init", "Error during readin parameters.", 1);
	}

	ifs.clear();
	ifs.seekg(0);

	filter_nonascii_and_comment(ifs, ascii_stream);
	ifs.clear();

	// file close after reading

    int ierr = 0;
    ascii_stream.rdstate();
    while (ascii_stream.good())
    {
        std::string word;
        ascii_stream >> word;
        ascii_stream.ignore(150, '\n');
        if (word == "INPUT_PARAMETERS")
        {
            ierr = 1;
            break;
        }
    }

    if (ierr == 0)
    {
        std::cout << " Error parameter list. "
                  << " The parameter list always starts with key word "
                     "'INPUT_PARAMETERS'. "
                  << std::endl;
        ModuleBase::WARNING_QUIT("Input", 
            "Bad parameter, please check the input parameters in file INPUT", 1);
    }

    ascii_stream.rdstate();
    while (ascii_stream.good())
    {
        std::string word; // temporary variable to store the keyword read-in
        ascii_stream >> word;
        if (ascii_stream.eof()) { break; }
        word = FmtCore::lower(word); // the lowercase of the keyword
        auto it = std::find_if(input_lists.begin(), input_lists.end(),
            [&word](const std::pair<std::string, Input_Item>& item) { return item.first == word; });
        if (it != this->input_lists.end()) // find the keyword
        {
            Input_Item* p_item = &(it->second);
            this->readvalue_items.push_back(p_item);
            if(p_item->is_read())
            {
                std::string warningstr = "The parameter " + p_item->label + " has been read twice.";
                ModuleBase::WARNING_QUIT("ReadInput", warningstr);
            }
            // qianrui delete '/' 2024-07-10, because path has '/' head.
            read_information(ascii_stream, p_item->str_values, "#!");
        }
        else // otherwise, it should be a comment or an unrecognized parameter
        {
            if (word[0] != '#' && word[0] != '/' && word[0] != '!') // if not recognized
            {
                std::cout << " THE PARAMETER NAME '" << word << "' IS INCORRECT!" << std::endl;
                ModuleBase::WARNING_QUIT("Input",
                    "Bad parameter, please check the input parameters in file INPUT", 1);
            }
            // otherwise, it is a comment. However, ...
            // but it is not always to be shorter than 150 characters
            // we can use ignore to skip the rest of the line
            ascii_stream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }

        ascii_stream.rdstate();
        if (ascii_stream.eof())
        {
            break;
        }
        else if (ascii_stream.bad())
        {
            ModuleBase::WARNING_QUIT("Input", 
                                    " Bad input parameters. ", 1);
        }
        else if (ascii_stream.fail())
        {
            ascii_stream.clear();
            ModuleBase::WARNING_QUIT("Input", 
                                    " fail to read parameters. ", 1);
        }
    }

    // 1) read the value of the parameters
    for (auto& readvalue_item: this->readvalue_items)
    {
        readvalue_item->read_value(*readvalue_item, param);
    }

    // 2) reset this value when some conditions are met
    //    e.g. if (calulation_type == "nscf") then set "init_chg" to "file".
    for (auto& input_item: this->input_lists)
    {
        Input_Item* resetvalue_item = &(input_item.second);
        if (resetvalue_item->reset_value != nullptr) 
		{
			resetvalue_item->reset_value(*resetvalue_item, param);
        }
    }
}

void ReadInput::write_txt_input(const Parameter& param, const std::string& filename)
{
    ModuleBase::TITLE("ReadInput", "write_txt_input");
    std::ofstream ofs(filename.c_str(), std::ios::out);
    ofs << "INPUT_PARAMETERS" << std::endl;
    ofs << std::setiosflags(std::ios::left);

    ofs << "# H0-only parameters" << std::endl;
    for (auto& item: this->input_lists)
    {
        Input_Item* p_item = &(item.second);
        if (p_item->get_final_value == nullptr) {
            continue;
}
        p_item->get_final_value(*p_item, param);
        ModuleBase::GlobalFunc::OUTP(ofs, p_item->label, p_item->final_value.str(), p_item->annotation);
    }
}

void ReadInput::check_ntype(const std::string& fn, int& param_ntype)
{
    std::ifstream ifa(fn.c_str(), std::ios::in);
    if (!ifa)
    {
        GlobalV::ofs_warning << fn;
        ModuleBase::WARNING_QUIT("ReadInput::check_ntype", "Can not find the file: " + fn);
    }

    int ntype_stru = 0;
    std::string temp;
    if (ModuleBase::GlobalFunc::SCAN_BEGIN(ifa, "ATOMIC_SPECIES"))
    {
        while (true)
        {
            ModuleBase::GlobalFunc::READ_VALUE(ifa, temp);
            if (temp == "LATTICE_CONSTANT" || temp == "NUMERICAL_ORBITAL" || temp == "NUMERICAL_DESCRIPTOR"
                || temp == "PAW_FILES" || ifa.eof())
            {
                break;
            }
            else if (!temp.empty() && std::isalpha(static_cast<unsigned char>(temp[0])))
            {
                ntype_stru += 1;
            }
        }
    }

    if (ntype_stru <= 0)
    {
        ModuleBase::WARNING_QUIT("ReadInput::check_ntype",
                                 "Failed to detect valid ntype from STRU: no valid ATOMIC_SPECIES entries were found.");
    }

    if (param_ntype < 0)
    {
        ModuleBase::WARNING_QUIT("ReadInput::check_ntype", "The ntype in INPUT should not be less than 0.");
    }
    else if (param_ntype != 0 && param_ntype != ntype_stru)
    {
        ModuleBase::WARNING_QUIT("ReadInput::check_ntype",
                                 "The ntype in INPUT is not equal to the ntype "
                                 "counted in STRU, check it.");
    }
    else if (param_ntype == 0)
    {
        param_ntype = ntype_stru;
        GlobalV::ofs_running << " 'ntype' is automatically set to " << param_ntype << std::endl;
    }
}

void ReadInput::add_item(const Input_Item& item)
{
    // Normally only rank 0 reads the input file
    // But rank -1 is used for help system mode where items should also be added
    if (this->rank == 0 || this->rank == -1)
    {
        this->input_lists.push_back(make_pair(item.label, item));
    }
}

} // namespace ModuleIO
