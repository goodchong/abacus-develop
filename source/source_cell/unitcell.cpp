#include <cstdlib>
#include <cstring> // Peize Lin fix bug about strcmp 2016-08-02

#include "source_base/constants.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "unitcell.h"
#include "bcast_cell.h"
#include "source_base/tool_quit.h"
#include "source_base/output.h"

#include "source_cell/read_stru.h"
#include "source_base/atom_in.h"
#include "source_base/element_elec_config.h"
#include "source_base/global_file.h"
#include "source_base/parallel_common.h"

#ifdef __MPI
#include "mpi.h"
#endif

UnitCell::UnitCell()
{
    itia2iat.create(1, 1);
}

UnitCell::~UnitCell()
{
    if (set_atom_flag)
    {
        delete[] atoms;
    }
}


void UnitCell::print_cell(std::ofstream& ofs) const {

    ModuleBase::GlobalFunc::OUT(ofs, "print_unitcell()");

    ModuleBase::GlobalFunc::OUT(ofs, "latName", latName);
    ModuleBase::GlobalFunc::OUT(ofs, "ntype", ntype);
    ModuleBase::GlobalFunc::OUT(ofs, "nat", nat);
    ModuleBase::GlobalFunc::OUT(ofs, "lat0", lat0);
    ModuleBase::GlobalFunc::OUT(ofs, "lat0_angstrom", lat0_angstrom);
    ModuleBase::GlobalFunc::OUT(ofs, "tpiba", tpiba);
    ModuleBase::GlobalFunc::OUT(ofs, "omega", omega);

    output::printM3(ofs, "Lattices Vector (R) : ", latvec);
    output::printM3(ofs, "Reciprocal lattice Vector (G): ", G);
    output::printM3(ofs, "GGT : ", GGT);

    ofs << std::endl;
    return;
}


void UnitCell::set_iat2itia() {
    assert(nat > 0);
    delete[] iat2it;
    delete[] iat2ia;
    this->iat2it = new int[nat];
    this->iat2ia = new int[nat];
    int iat = 0;
    for (int it = 0; it < ntype; it++) {
        for (int ia = 0; ia < atoms[it].na; ia++) {
            this->iat2it[iat] = it;
            this->iat2ia[iat] = ia;
            ++iat;
        }
    }
    return;
}


//==============================================================
// Calculate various lattice related quantities for given latvec
//==============================================================
void UnitCell::setup_cell(const std::string& fn,
                          std::ofstream& log,
                          const int nspin,
                          const std::string& orbital_dir,
                          const bool noncolin)
{
    ModuleBase::TITLE("UnitCell", "setup_cell");

    assert(ntype > 0);

    // (1) init *Atom class array.
    this->atoms = new Atom[this->ntype]; // atom species.
    this->set_atom_flag = true;

    bool ok = true;
    bool ok2 = true;

    // (3) read in atom information
    this->atom_mass.resize(ntype);
    this->atom_label.resize(ntype);
    this->pseudo_fn.resize(ntype);
    this->pseudo_type.resize(ntype);
    this->orbital_fn.resize(ntype);

    if (GlobalV::MY_RANK == 0)
    {
        // open "atom_unitcell" file.
        std::ifstream ifa(fn.c_str(), std::ios::in);
        if (!ifa)
        {
            GlobalV::ofs_warning << fn;
            ok = false;
        }

        if (ok)
        {
            log << "\n\n";
            log << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
            log << " |                                                                    |" << std::endl;
            log << " |                        #Setup Unitcell#                            |" << std::endl;
            log << " | From the input file and the structure file we know the number of   |" << std::endl;
            log << " | different elments in this unitcell, then we list the detail        |" << std::endl;
            log << " | information for each element, especially the zeta and polar atomic |" << std::endl;
            log << " | orbital number for each element. The total atom number is counted. |" << std::endl;
            log << " | We calculate the nearest atom distance for each atom and show the  |" << std::endl;
            log << " | Cartesian and Direct coordinates for each atom. We list the file   |" << std::endl;
            log << " | address for atomic orbitals. The volume and the lattice vectors    |" << std::endl;
            log << " | in real and reciprocal space is also shown.                        |" << std::endl;
            log << " |                                                                    |" << std::endl;
            log << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
            log << "\n";

            log << " READING UNITCELL INFORMATION" << std::endl;
            //========================
            // call read_atom_species
            //========================
            const bool read_atom_species
                = unitcell::read_atom_species(ifa, log, *this);
            //========================
            // call read_lattice_constant
            //========================
            const bool read_lattice_constant = unitcell::read_lattice_constant(ifa, log ,this->lat);
            ok2 = unitcell::read_atom_positions(*this,
                                                ifa,
                                                log,
                                                GlobalV::ofs_warning,
                                                nspin,
                                                orbital_dir,
                                                noncolin);
        }
    }
#ifdef __MPI
    Parallel_Common::bcast_bool(ok);
    Parallel_Common::bcast_bool(ok2);
#endif
    if (!ok) {
        ModuleBase::WARNING_QUIT(
            "UnitCell::setup_cell",
            "Can not find the file containing atom positions.!");
    }
    if (!ok2) {
        ModuleBase::WARNING_QUIT("UnitCell::setup_cell",
                                 "Something wrong during read_atom_positions.");
    }
#ifdef __MPI
    unitcell::bcast_unitcell(*this, nspin);
#endif

    //========================================================
    // Calculate unit cell volume
    // the reason to calculate volume here is
    // Firstly, latvec must be read in.
    //========================================================
    assert(lat0 > 0.0);
    this->omega = latvec.Det() * this->lat0 * this->lat0 * this->lat0;


    if (this->omega < 0)
    {
        std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
        std::cout << " Warning: The lattice vector is left-handed; a right-handed vector is prefered." << std::endl;
        std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
        GlobalV::ofs_warning <<
        "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
        GlobalV::ofs_warning <<
        " Warning: The lattice vector is left-handed; a right-handed vector is prefered." << std::endl;
        GlobalV::ofs_warning <<
        "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
        this->omega = std::abs(this->omega);
    }
    else if (this->omega == 0)
    {
        ModuleBase::WARNING_QUIT("setup_cell", "The volume is zero.");
    }
    else
    {
        ModuleBase::GlobalFunc::OUT(log, "Cell volume (Bohr^3)", this->omega);
        ModuleBase::GlobalFunc::OUT(log, "Cell volume (A^3)", this->omega * pow(ModuleBase::BOHR_TO_A, 3));
    }

    //==========================================================
    // Calculate recip. lattice vectors and dot products
    // latvec have the unit of lat0, but G has the unit 2Pi/lat0
    //==========================================================
    this->GT = latvec.Inverse();
    this->G = GT.Transpose();
    this->GGT = G * GT;
    this->invGGT = GGT.Inverse();

    log << std::endl;
    output::printM3(log,
                    "Lattice vectors: (Cartesian coordinate: in unit of a_0)",
                    latvec);
    output::printM3(
        log,
        "Reciprocal vectors: (Cartesian coordinate: in unit of 2 pi/a_0)",
        G);

    //===================================
    // set index for iat2it, iat2ia
    //===================================
    this->set_iat2itia();

    return;
}


void UnitCell::set_iat2iwt(const int& npol_in)
{
#ifdef __DEBUG
    assert(npol_in == 1 || npol_in == 2);
    assert(this->nat > 0);
    assert(this->ntype > 0);
#endif
    this->iat2iwt.resize(this->nat);
    this->npol = npol_in;
    int iat = 0;
    int iwt = 0;

    for (int it = 0; it < this->ntype; it++)
    {
        for (int ia = 0; ia < atoms[it].na; ia++)
        {
            this->iat2iwt[iat] = iwt;
            iwt += atoms[it].nw * this->npol;
            ++iat;
        }
    }
    return;
}



void UnitCell::setup(const int ntype_in)
{
    this->ntype = ntype_in;
    this->magnet.start_mag.resize(ntype_in, 0.0);
}


void UnitCell::compare_atom_labels(const std::string& label1, const std::string& label2) const
{
    if (label1!= label2) //'!( "Ag" == "Ag" || "47" == "47" || "Silver" == Silver" )'
    {
        atom_in ai;
        if (!(std::to_string(ai.atom_Z[label1]) == label2
              ||                                  // '!( "Ag" == "47" )'
              ai.atom_symbol[label1] == label2 || // '!( "Ag" == "Silver" )'
              label1 == std::to_string(ai.atom_Z[label2])
              || // '!( "47" == "Ag" )'
              label1 == std::to_string(ai.symbol_Z[label2])
              ||                                  // '!( "47" == "Silver" )'
              label1 == ai.atom_symbol[label2] || // '!( "Silver" == "Ag" )'
              std::to_string(ai.symbol_Z[label1])
                  == label2)) // '!( "Silver" == "47" )'
        {
            std::string stru_label = "";
            std::string psuedo_label = "";
			for (int ip = 0; ip < label1.length(); ip++)
			{
				if (!(isdigit(label1[ip]) || label1[ip] == '_'))
				{
					stru_label += label1[ip];
				}
				else
				{
					break;
				}
			}
			stru_label[0] = toupper(stru_label[0]);

			for (int ip = 0; ip < label2.length(); ip++)
			{
				if (!(isdigit(label2[ip]) || label2[ip] == '_'))
				{
					psuedo_label += label2[ip];
				}
				else
				{
					break;
				}
			}
            psuedo_label[0] = toupper(psuedo_label[0]);

            if (!(stru_label == psuedo_label
                  || //' !("Ag1" == "ag_locpsp" || "47" == "47" || "Silver" ==
                     //Silver" )'
                  std::to_string(ai.atom_Z[stru_label]) == psuedo_label
                  || // ' !("Ag1" == "47" )'
                  ai.atom_symbol[stru_label] == psuedo_label
                  || // ' !("Ag1" == "Silver")'
                  stru_label == std::to_string(ai.atom_Z[psuedo_label])
                  || // ' !("47" == "Ag1" )'
                  stru_label == std::to_string(ai.symbol_Z[psuedo_label])
                  || // ' !("47" == "Silver1" )'
                  stru_label == ai.atom_symbol[psuedo_label]
                  || // ' !("Silver1" == "Ag" )'
                  std::to_string(ai.symbol_Z[stru_label])
                      == psuedo_label)) // ' !("Silver1" == "47" )'

            {
                std::string atom_label_in_orbtial
                    = "atom label in orbital file ";
                std::string mismatch_with_pseudo
                    = " mismatch with pseudo file of ";
                ModuleBase::WARNING_QUIT("UnitCell::read_pseudo",
                                         atom_label_in_orbtial + label1
                                             + mismatch_with_pseudo + label2);
            }
        }
    }
}
