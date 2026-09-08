#ifndef UNITCELL_H
#define UNITCELL_H

#include <memory>
#include "source_base/global_function.h"
#include "source_cell/atom_spec.h"
#include "source_cell/magnetism.h"
#include "source_cell/nonlocal_info_base.h"
#include "source_cell/unitcell_data.h"

// provide the basic information about unitcell.
class UnitCell {
  public:
    double get_lat0() const {
        return lat0;
    }

    double get_omega() const {
        return omega;
    }

    const ModuleBase::Matrix3& get_latvec() const {
        return latvec;
    }

    int get_natom() const {
        return nat;
    }

    int get_na(int i) const {
        return atoms[i].na;
    }

    int get_ntype() const {
        return ntype;
    }

    ModuleBase::Vector3<double> get_tau(int i, int j) const {
        return atoms[i].tau[j];
    }

    Atom* atoms = nullptr;
    bool set_atom_flag = false;                     // added on 2009-3-8 by mohan
    Magnetism magnet;                               // magnetism Yu Liu 2021-07-03

    Lattice lat;
    std::string& Coordinate = lat.Coordinate;
    std::string& latName = lat.latName;
    double& lat0 = lat.lat0;
    double& lat0_angstrom = lat.lat0_angstrom;
    double& tpiba = lat.tpiba;
    double& tpiba2 = lat.tpiba2;
    double& omega = lat.omega;
    ModuleBase::Matrix3& latvec = lat.latvec;
    ModuleBase::Vector3<double>&a1 = lat.a1, &a2 = lat.a2, &a3 = lat.a3;
    ModuleBase::Vector3<double>& latcenter = lat.latcenter;
    ModuleBase::Matrix3& G = lat.G;
    ModuleBase::Matrix3& GT = lat.GT;
    ModuleBase::Matrix3& GGT = lat.GGT;
    ModuleBase::Matrix3& invGGT = lat.invGGT;

    Statistics st;
    int& ntype = st.ntype;
    int& nat = st.nat;
    int*& iat2it = st.iat2it;
    int*& iat2ia = st.iat2ia;
    ModuleBase::IntArray& itia2iat = st.itia2iat;
    int& namax = st.namax;

    // ========================================================
    // iat2iwt is the atom index iat to the first global index for orbital of
    // this atom the size of iat2iwt is nat, the value should be
    // sum_{i=0}^{iat-1} atoms[it].nw * npol where the npol is the number of
    // polarizations, 1 for non-magnetic(NSPIN=1 or 2), 2 for magnetic(only
    // NSPIN=4) this part only used for Atomic Orbital based calculation
    // ========================================================
  public:
    // indexing tool for find orbital global index from it,ia,iw
    template <typename Tiait>
    inline Tiait
        itiaiw2iwt(const Tiait& it, const Tiait& ia, const Tiait& iw) const {
        return Tiait(this->iat2iwt[this->itia2iat(it, ia)] + iw);
    }
    // initialize iat2iwt
    void set_iat2iwt(const int& npol_in);
    // get iat2iwt
    inline const int* get_iat2iwt() const { return iat2iwt.data(); }
    // get npol
    inline const int& get_npol() const { return npol; }

  private:
    std::vector<int> iat2iwt; // iat ==> iwt, the first global index for orbital of this atom
    int npol = 1; // number of spin polarizations, initialized in set_iat2iwt
                  // ----------------- END of iat2iwt part -----------------

  public:
    //========================================================
    // indexing tools for ia and it
    // return true if the last out is reset
    //========================================================
    template <typename Tiat, typename Tiait>
    inline bool iat2iait(const Tiat iat, Tiait* ia, Tiait* it) const {
        if (iat >= nat) {
            *ia = 0;
            *it = ntype;
            return false;
        }
        *ia = (Tiait)iat2ia[iat];
        *it = (Tiait)iat2it[iat];
        return true;
    }

    template <typename Tiat, typename Tiait>
    inline bool ijat2iaitjajt(const Tiat ijat,
                              Tiait* ia,
                              Tiait* it,
                              Tiait* ja,
                              Tiait* jt) const {
        Tiat iat = ijat / nat;
        Tiat jat = ijat % nat;
        iat2iait(iat, ia, it);
        iat2iait(jat, ja, jt);
        return true;
    }

    template <typename Tiait>
    inline bool step_it(Tiait* it) const {
        if (++(*it) >= ntype) {
            *it = 0;
            return true;
        }
        return false;
    }

    template <typename Tiait>
    inline bool step_ia(const Tiait it, Tiait* ia) const {
        if (++(*ia) >= atoms[it].na) {
            *ia = 0;
            return true;
        }
        return false;
    }

    template <typename Tiait>
    inline bool step_iait(Tiait* ia, Tiait* it) const {
        if (step_ia(*it, ia)) {
            return step_it(it);
        }
        return false;
    }

    template <typename Tiait>
    inline bool
        step_jajtiait(Tiait* ja, Tiait* jt, Tiait* ia, Tiait* it) const {
        if (step_iait(ja, jt)) {
            return step_iait(ia, it);
        }
        return false;
    }

    // get tau for atom iat
    inline const ModuleBase::Vector3<double>& get_tau(const int& iat) const {
        return atoms[iat2it[iat]].tau[iat2ia[iat]];
    }

    // calculate vector between two atoms with R cell
    inline const ModuleBase::Vector3<double>
        cal_dtau(const int& iat1,
                 const int& iat2,
                 const ModuleBase::Vector3<int>& R) const {
        return get_tau(iat2) + double(R.x) * a1 + double(R.y) * a2
               + double(R.z) * a3 - get_tau(iat1);
    }

    //============================================================
    // meshx : max number of mesh point in pseudopotential file
    // meshx : max number of mesh points in pseudopotential file
    //============================================================
    int meshx = 0;

  public:
    UnitCell();
    ~UnitCell();
    void print_cell(std::ofstream& ofs) const;

    std::vector<double>      atom_mass;
    std::vector<std::string> atom_label;
    std::vector<std::string> pseudo_fn;
    std::vector<std::string> pseudo_type;

    std::vector<std::string> orbital_fn;  // filenames of orbitals, liuyu add 2022-10-19
    void set_iat2itia();

    void setup_cell(const std::string& fn,
                    std::ofstream& log,
                    int nspin,
                    const std::string& orbital_dir,
                    bool noncolin);

    /**
     * @brief Pointer to non-local pseudopotential information.
     *
     * This pointer is set during LCAO initialization and provides access
     * to non-local projector data. It is null for non-LCAO calculations.
     */
    std::unique_ptr<NonlocalInfoBase> infoNL;

    void setup(int ntype_in);

    /// @brief check consistency between two atom labels from STRU and pseudo or
    /// orb file
    void compare_atom_labels(const std::string& label1, const std::string& label2) const;
};

#endif // unitcell class
