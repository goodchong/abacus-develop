#include "xc_functional.h"

#include "source_base/global_function.h"
#include "source_base/tool_title.h"

#include <algorithm>
#include <cctype>
#include <sstream>

XC_Functional::XC_Functional() {}

XC_Functional::~XC_Functional() {}

std::vector<int> XC_Functional::func_id(1);
int XC_Functional::func_type = 0;

void XC_Functional::set_xc_type(const std::string xc_func_in)
{
    ModuleBase::TITLE("XC_Functional", "set_xc_type");

    std::string xc_func = xc_func_in;
    std::transform(xc_func.begin(), xc_func.end(), xc_func.begin(),
                   [](const unsigned char c) { return static_cast<char>(std::toupper(c)); });

    func_id.clear();
    if (xc_func == "LDA" || xc_func == "PZ" || xc_func == "SLAPZNOGXNOGC")
    {
        func_id.push_back(XC_LDA_X);
        func_id.push_back(XC_LDA_C_PZ);
        func_type = 1;
    }
    else if (xc_func == "PWLDA")
    {
        func_id.push_back(XC_LDA_X);
        func_id.push_back(XC_LDA_C_PW);
        func_type = 1;
    }
    else if (xc_func == "PBE" || xc_func == "SLAPWPBXPBC")
    {
        func_id.push_back(XC_GGA_X_PBE);
        func_id.push_back(XC_GGA_C_PBE);
        func_type = 2;
    }
    else if (xc_func == "PBESOL")
    {
        func_id.push_back(XC_GGA_X_PBE_SOL);
        func_id.push_back(XC_GGA_C_PBE_SOL);
        func_type = 2;
    }
    else if (xc_func == "REVPBE")
    {
        func_id.push_back(XC_GGA_X_PBE_R);
        func_id.push_back(XC_GGA_C_PBE);
        func_type = 2;
    }
    else if (xc_func == "WC")
    {
        func_id.push_back(XC_GGA_X_WC);
        func_id.push_back(XC_GGA_C_PBE);
        func_type = 2;
    }
    else if (xc_func == "BLYP")
    {
        func_id.push_back(XC_GGA_X_B88);
        func_id.push_back(XC_GGA_C_LYP);
        func_type = 2;
    }
    else if (xc_func == "BP")
    {
        func_id.push_back(XC_GGA_X_B88);
        func_id.push_back(XC_GGA_C_P86);
        func_type = 2;
    }
    else if (xc_func == "PW91")
    {
        func_id.push_back(XC_GGA_X_PW91);
        func_id.push_back(XC_GGA_C_PW91);
        func_type = 2;
    }
    else if (xc_func == "HCTH")
    {
        func_id.push_back(XC_GGA_X_HCTH_A);
        func_id.push_back(XC_GGA_C_HCTH_A);
        func_type = 2;
    }
    else if (xc_func == "OLYP")
    {
        func_id.push_back(XC_GGA_X_OPTX);
        func_id.push_back(XC_GGA_C_LYP);
        func_type = 2;
    }
    else
    {
        ModuleBase::WARNING_QUIT(
            "XC_Functional::set_xc_type",
            "Only the built-in LDA and GGA functionals supported by the H0 executable are available; got '"
                + xc_func + "'.");
    }
}

std::string XC_Functional::output_info()
{
    std::ostringstream output;
    output << " XC (built-in):";
    for (std::vector<int>::const_iterator id = func_id.begin(); id != func_id.end(); ++id)
    {
        output << ' ' << *id;
    }
    return output.str();
}
