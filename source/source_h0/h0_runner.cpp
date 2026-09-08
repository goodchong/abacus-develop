#include "h0_runner.h"

#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#include "source_base/tool_title.h"
#include "source_basis/module_ao/ORB_read.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_basis/module_nao/two_center_bundle.h"
#include "source_cell/module_neighbor/sltk_atom_arrange.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_estate/cal_ux.h"
#include "source_estate/module_charge/charge.h"
#include "source_estate/module_pot/potential_new.h"
#include "source_estate/param_update.h"
#include "source_estate/read_pseudo.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_io/module_hs/write_h0.h"
#include "source_io/module_parameter/input_parameter.h"
#include "source_io/module_parameter/parameter.h"
#include "source_lcao/module_gint/gint.h"
#include "source_lcao/module_gint/gint_info.h"
#include "source_lcao/module_gint/gint_interface.h"
#include "source_lcao/module_hcontainer/hcontainer.h"
#include "source_lcao/module_operator_lcao/ekinetic.h"
#include "source_lcao/module_operator_lcao/nonlocal.h"
#include "source_lcao/lcao_h0_basis.h"
#include "source_lcao/record_adj.h"
#include "source_pw/module_pwdft/parallel_grid.h"
#include "source_pw/module_pwdft/setup_pwrho.h"
#include "source_pw/module_pwdft/structure_factor.h"
#include "source_pw/module_pwdft/vl_pw.h"

#ifdef __MPI
#include "source_base/module_external/blacs_connector.h"
#endif

#include <complex>
#include <string>
#include <vector>

namespace ModuleH0
{

struct H0Runner::Impl
{
    bool pw_rho_allocated = false;
    ModulePW::PW_Basis* pw_rho = nullptr;
    ModulePW::PW_Basis* pw_rhod = nullptr;
    ModulePW::PW_Basis_Big* pw_big = nullptr;
    Parallel_Grid pgrid;
    Structure_Factor structure_factor;
    Charge charge;
    pseudopot_cell_vl local_potential;
    std::unique_ptr<elecstate::Potential> potential;
    double etxc = 0.0;
    double vtxc = 0.0;

    Parallel_Orbitals para_v;
    TwoCenterBundle two_center_bundle;
    LCAO_Orbitals orbitals;
    Grid_Driver grid_driver;
    Record_adj record_adj;
    std::unique_ptr<ModuleGint::GintInfo> gint_info;
};

namespace
{

void validate_xc()
{
    const int functional_type = XC_Functional::get_func_type();
    if (functional_type != 1 && functional_type != 2)
    {
        ModuleBase::WARNING_QUIT("H0Runner",
                                 "Only built-in LDA and GGA functionals are supported for H0 output.");
    }
}

template <typename TK, typename TR>
hamilt::HContainer<TR> build_fixed_h0(const UnitCell& ucell,
                                      const Grid_Driver& grid_driver,
                                      const Parallel_Orbitals& para_v,
                                      const LCAO_Orbitals& orbitals,
                                      const TwoCenterBundle& two_center_bundle)
{
    hamilt::HContainer<TR> kinetic_h(&para_v);
    GlobalV::ofs_running << " H0 component: allocate T\n" << std::flush;
    hamilt::EKinetic<hamilt::OperatorLCAO<TK, TR>> kinetic(&kinetic_h,
                                                           &ucell,
                                                           orbitals.cutoffs(),
                                                           &grid_driver,
                                                           two_center_bundle.kinetic_orb.get());
    GlobalV::ofs_running << " H0 component: integrate T\n" << std::flush;
    kinetic.contributeHR();

    hamilt::HContainer<TR> nonlocal_h(&para_v);
    GlobalV::ofs_running << " H0 component: allocate V_nl\n" << std::flush;
    hamilt::Nonlocal<hamilt::OperatorLCAO<TK, TR>> nonlocal(&nonlocal_h,
                                                            &ucell,
                                                            orbitals.cutoffs(),
                                                            &grid_driver,
                                                            two_center_bundle.overlap_orb_beta.get());
    GlobalV::ofs_running << " H0 component: integrate V_nl\n" << std::flush;
    nonlocal.contributeHR();

    GlobalV::ofs_running << " H0 component: union T and V_nl\n" << std::flush;
    kinetic_h.add_value_union(nonlocal_h);
    return kinetic_h;
}

void add_local_potential(hamilt::HContainer<double>& h0,
                         ModuleGint::GintInfo& gint_info,
                         const elecstate::Potential& potential,
                         const int spin)
{
    const std::vector<int>& local_ijr = gint_info.get_ijr_info();
    hamilt::HContainer<double> local_h(h0.get_paraV(), nullptr, &local_ijr);
    ModuleGint::cal_gint_vl(potential.get_eff_v(spin), &local_h);
    h0.add_value_union(local_h);
}

void add_local_potential(hamilt::HContainer<std::complex<double>>& h0,
                         ModuleGint::GintInfo& gint_info,
                         const elecstate::Potential& potential)
{
    GlobalV::ofs_running << " H0 component: allocate spinor local matrix\n" << std::flush;
    const std::vector<int>& local_ijr = gint_info.get_ijr_info();
    hamilt::HContainer<std::complex<double>> local_h(h0.get_paraV(), nullptr, &local_ijr);
    std::vector<const double*> effective_potential(4, nullptr);
    for (int is = 0; is < 4; ++is)
    {
        effective_potential[is] = potential.get_eff_v(is);
    }
    GlobalV::ofs_running << " H0 component: integrate spinor local matrix\n" << std::flush;
    ModuleGint::cal_gint_vl(effective_potential, &local_h);
    GlobalV::ofs_running << " H0 component: union spinor local matrix\n" << std::flush;
    h0.add_value_union(local_h);
}

} // namespace

H0Runner::H0Runner() : impl_(new Impl) {}

H0Runner::~H0Runner()
{
    impl_->potential.reset();
    pw::teardown_pwrho(impl_->pw_rho_allocated,
                       impl_->pw_rho,
                       impl_->pw_rhod);
    impl_->pw_big = nullptr;
}

void H0Runner::initialize(UnitCell& ucell, const Input_para& inp)
{
    ModuleBase::TITLE("H0Runner", "initialize");
    ModuleBase::timer::start("H0Runner", "initialize");

    const auto atoms_info
        = elecstate::read_pseudo(GlobalV::ofs_running,
                                 ucell,
                                 inp.pseudo_dir,
                                 inp.dft_functional,
                                 inp.lspinorb,
                                 inp.pseudo_rcut,
                                 inp.soc_lambda,
                                 inp.nspin,
                                 PARAM.globalv.npol,
                                 PARAM.globalv.two_fermi,
                                 inp.nelec,
                                 inp.nupdown);
    elecstate::ParamUpdater::update_from_atoms_info(atoms_info);

    pw::setup_pwrho(ucell,
                    impl_->pw_rho_allocated,
                    impl_->pw_rho,
                    impl_->pw_rhod,
                    impl_->pw_big,
                    inp);

    impl_->structure_factor.set(impl_->pw_rhod, inp.nbspline);
    impl_->pgrid.init(impl_->pw_rhod->nx,
                      impl_->pw_rhod->ny,
                      impl_->pw_rhod->nz,
                      impl_->pw_rhod->nplane,
                      impl_->pw_rhod->nrxx,
                      impl_->pw_big->nbz,
                      impl_->pw_big->bz);
    impl_->structure_factor.setup(&ucell, impl_->pgrid, impl_->pw_rhod);

    XC_Functional::set_xc_type(ucell.atoms[0].ncpp.xc_func);
    validate_xc();
    GlobalV::ofs_running << XC_Functional::output_info() << std::endl;

    if (inp.h0_type == "full")
    {
        impl_->charge.set_rhopw(impl_->pw_rhod);
        impl_->charge.allocate(inp.nspin);
        impl_->charge.init_rho(ucell,
                               impl_->pgrid,
                               impl_->structure_factor.strucFac);
        impl_->charge.check_rho();
    }

    LCAO_domain::init_basis_lcao(impl_->para_v,
                                 inp.lcao_ecut,
                                 inp.lcao_dk,
                                 ucell,
                                 impl_->two_center_bundle,
                                 impl_->orbitals);

    impl_->local_potential.init_vloc(ucell, impl_->pw_rho);
    impl_->potential.reset(new elecstate::Potential(impl_->pw_rhod,
                                                    &ucell,
                                                    &impl_->local_potential.vloc,
                                                    &impl_->structure_factor,
                                                    &impl_->etxc,
                                                    &impl_->vtxc));

    ModuleBase::timer::end("H0Runner", "initialize");
}

void H0Runner::run(UnitCell& ucell)
{
    ModuleBase::TITLE("H0Runner", "run");
    ModuleBase::timer::start("H0Runner", "run");

    const double search_radius
        = atom_arrange::set_sr_NL(GlobalV::ofs_running,
                                  "ie",
                                  impl_->orbitals.get_rcutmax_Phi(),
                                  ucell.infoNL->get_rcutmax_Beta(),
                                  false);
    atom_arrange::search(PARAM.globalv.search_pbc,
                         GlobalV::ofs_running,
                         impl_->grid_driver,
                         ucell,
                         search_radius,
                         false);

    impl_->gint_info.reset(new ModuleGint::GintInfo(impl_->pw_big->nbx,
                                                    impl_->pw_big->nby,
                                                    impl_->pw_big->nbz,
                                                    impl_->pw_rho->nx,
                                                    impl_->pw_rho->ny,
                                                    impl_->pw_rho->nz,
                                                    0,
                                                    0,
                                                    impl_->pw_big->nbzp_start,
                                                    impl_->pw_big->nbx,
                                                    impl_->pw_big->nby,
                                                    impl_->pw_big->nbzp,
                                                    impl_->orbitals.Phi,
                                                    ucell,
                                                    impl_->grid_driver));
    ModuleGint::Gint::set_gint_info(impl_->gint_info.get());

    impl_->record_adj.for_2d(ucell,
                             impl_->grid_driver,
                             impl_->para_v,
                             impl_->orbitals.cutoffs());
    GlobalV::ofs_running << " H0 phase: initialize spin frame\n" << std::flush;
    elecstate::cal_ux(ucell, PARAM.inp.nspin);

    std::vector<std::string> potential_components(1, "local");
    if (PARAM.inp.h0_type == "full")
    {
        potential_components.push_back("hartree");
        potential_components.push_back("xc");
        impl_->charge.set_rho_core(ucell,
                                   impl_->structure_factor.strucFac,
                                   impl_->local_potential.numeric);
        impl_->charge.renormalize_rho();
    }
    GlobalV::ofs_running << " H0 phase: construct local effective potential\n" << std::flush;
    impl_->potential->pot_register(potential_components);
    impl_->potential->init_pot(PARAM.inp.h0_type == "full" ? &impl_->charge : nullptr);

    const int output_nspin = (PARAM.inp.nspin == 2 && PARAM.inp.h0_type == "full") ? 2 : 1;
    GlobalV::ofs_running << " H0 phase: construct T and V_nl\n" << std::flush;
    if (PARAM.inp.nspin < 4)
    {
        hamilt::HContainer<double> fixed_h
            = build_fixed_h0<std::complex<double>, double>(ucell,
                                                            impl_->grid_driver,
                                                            impl_->para_v,
                                                            impl_->orbitals,
                                                            impl_->two_center_bundle);
        for (int is = 0; is < output_nspin; ++is)
        {
            hamilt::HContainer<double> h0(fixed_h);
            h0.add_value_intersection(fixed_h);
            add_local_potential(h0, *impl_->gint_info, *impl_->potential, is);
            ModuleIO::write_h0(h0,
                               ucell,
                               impl_->para_v,
                               PARAM.globalv.global_out_dir,
                               PARAM.inp.h0_type,
                               PARAM.inp.nspin,
                               output_nspin,
                               is,
                               PARAM.inp.h0_sparse_threshold,
                               PARAM.inp.h0_precision,
                               GlobalV::MY_RANK);
        }
    }
    else
    {
        hamilt::HContainer<std::complex<double>> h0
            = build_fixed_h0<std::complex<double>, std::complex<double>>(ucell,
                                                                          impl_->grid_driver,
                                                                          impl_->para_v,
                                                                          impl_->orbitals,
                                                                          impl_->two_center_bundle);
        GlobalV::ofs_running << " H0 component: project V_loc/V_H/V_xc (spinor)\n" << std::flush;
        add_local_potential(h0, *impl_->gint_info, *impl_->potential);
        GlobalV::ofs_running << " H0 phase: write spinor CSR\n" << std::flush;
        ModuleIO::write_h0(h0,
                           ucell,
                           impl_->para_v,
                           PARAM.globalv.global_out_dir,
                           PARAM.inp.h0_type,
                           PARAM.inp.nspin,
                           1,
                           0,
                           PARAM.inp.h0_sparse_threshold,
                           PARAM.inp.h0_precision,
                           GlobalV::MY_RANK);
    }

    GlobalV::ofs_running << " H0 matrix construction and export completed; exiting one-shot workflow.\n";
    ModuleBase::timer::end("H0Runner", "run");
}

void H0Runner::finalize()
{
#ifdef __MPI
    Cblacs_exit(1);
#endif
}

} // namespace ModuleH0
