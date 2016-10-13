/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 */

/**
 * @file   wimdiscr.hpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Mon Aug  3 11:53:19 2015
 */


#ifndef __WIMDISCR_H
#define __WIMDISCR_H 1

#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <boost/multi_array.hpp>
#include <boost/array.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/any.hpp>
#include <boost/format.hpp>
#include <boost/mpi/timer.hpp>
#include <InterpFromGridToMeshx.h>
#include <date.hpp>
#include <iomanip>
#include <omp.h>

#ifdef __cplusplus
extern "C"
{
#endif
#include<RTparam_outer.h>
#include <mapx.h>
#ifdef __cplusplus
}
#endif

#define PI M_PI

namespace Wim
{

namespace po = boost::program_options;
namespace fs = boost::filesystem;

template<typename T=float> class WimDiscr
{
	typedef T value_type;
    typedef size_t size_type;
	typedef boost::multi_array<value_type, 2> array2_type;
    typedef boost::multi_array<value_type, 3> array3_type;
    typedef boost::multi_array<value_type, 4> array4_type;
    typedef typename array2_type::index index;

public:

    typedef struct BreakInfo
    {
        // information needed for breaking
        value_type conc;   // concentration
        value_type thick;  // thickness
        value_type mom0;
        value_type mom2;
        value_type var_strain;
        mutable value_type dfloe;
        mutable bool broken;
    } BreakInfo;

    typedef struct WimGrid
    {
        // information describing wim grid
        int nx;
        int ny;
        value_type dx;
        value_type dy;
        std::vector<value_type> X;
        std::vector<value_type> Y;
    } WimGrid;


public:

    WimDiscr()
        :
        vm(),
        nx(),
        ny()
    {}

    WimDiscr(po::variables_map const& vmIn)
        :
        vm(vmIn),
        nx(vm["wim.nx"].template as<int>()),
        ny(vm["wim.ny"].template as<int>())
    {}

    void gridProcessing();
    void saveGrid();
    void readGridFromFile();
    void readDataFromFile(std::string const& filein);
    void exportResults(std::string const& output_type, value_type const& t_out) const;
    void saveLog(value_type const& t_out) const;
    void init();

    void assign(std::vector<value_type> const& ice_c = std::vector<value_type>(),
                std::vector<value_type> const& ice_h = std::vector<value_type>(),
                std::vector<value_type> const& n_floes = std::vector<value_type>(),
                std::vector<value_type> const& swh_in = std::vector<value_type>(),
                std::vector<value_type> const& mwp_in = std::vector<value_type>(),
                std::vector<value_type> const& mwd_in = std::vector<value_type>(),
                bool step = false);

    void timeStep(bool step = false);

    //void doBreaking(value_type const& mom0, value_type const& mom2, value_type const& var_strain);
    void doBreaking(BreakInfo const& breakinfo);

    void setMesh(std::vector<value_type> const& m_rx, std::vector<value_type> const& m_ry,
                 std::vector<value_type> const& m_conc, std::vector<value_type> const& m_thick, std::vector<value_type> const& m_dfloe);
    void clearMesh();

    WimGrid wimGrid() const { return wim_grid; }

    void test(value_type* toto);

    void run(std::vector<value_type> const& ice_c = std::vector<value_type>(),
             std::vector<value_type> const& ice_h = std::vector<value_type>(),
             std::vector<value_type> const& n_floes = std::vector<value_type>(),
             std::vector<value_type> const& swh_in = std::vector<value_type>(),
             std::vector<value_type> const& mwp_in = std::vector<value_type>(),
             std::vector<value_type> const& mwd_in = std::vector<value_type>(),
             bool step = false);

    //===========================================================================
    //FSD: Dmax -> <D^moment> conversion
    void floeScaling(
          value_type const& dmax, int const& moment, value_type& dave);
    void floeScalingSmooth(
          value_type const& dmax, int const& moment, value_type& dave);
    //===========================================================================

    //===========================================================================
    //advection/attenuation
    void advAttenSimple(
          array3_type& Sdir, array2_type& Sfreq,
          array2_type& taux_omega,array2_type& tauy_omega,
          array2_type const& ag2d_eff);
    void advAttenIsotropic(array3_type& Sdir, array2_type& Sfreq,
          array2_type& taux_omega,array2_type& tauy_omega,
          array2_type const& ag2d_eff);
    void waveAdvWeno(
          array2_type& h, array2_type const& u, array2_type const& v);
    void weno3pdV2(
          array2_type const& gin, array2_type const& u, array2_type const& v,
          array2_type const& scuy, array2_type const& scvx,
          array2_type const& scp2i, array2_type const& scp2,
          array2_type& saoout);
    void padVar(array2_type const& u, array2_type& upad);
    //===========================================================================


    void calcMWD();
    void idealWaveFields(array2_type& wave_mask,value_type const xfac);
    void idealIceFields(array2_type& ice_mask,value_type const xfac);
    //void getWimCenters(value_type& x,value_type& y,value_type const& rotangle);

    value_type thetaDirFrac(value_type const& th1_, value_type const& dtheta_, value_type const& mwd_);
    value_type thetaInRange(value_type const& th_, value_type const& th1, bool const& close_on_right=false);

    array2_type getX() const { return X_array; }
    array2_type getY() const { return Y_array; }
    array2_type getSCUY() const { return SCUY_array; }
    array2_type getSCVX() const { return SCVX_array; }
    array2_type getSCP2() const { return SCP2_array; }
    array2_type getSCP2I() const { return SCP2I_array; }
    array2_type getLANDMASK() const { return LANDMASK_array; }

    std::string getWimGridFilename() const { return wim_gridfile; }
    //std::vector<int> getWimShape();

    std::vector<value_type> getTaux() const { return tau_x; }
    std::vector<value_type> getTauy() const { return tau_y; }
    std::vector<value_type> getNFloes() const { return nfloes; }


private:

    po::variables_map vm;
    int nx, ny, nxext, nyext, nbdy, nbdx, nghost;
    int wim_itest, wim_jtest;
    array2_type X_array, Y_array, SCUY_array, SCVX_array,
                SCP2_array, SCP2I_array, LANDMASK_array;

    value_type cfl, dom, guess, Hs_inc, Tp_inc, mwd_inc, Tmin, Tmax, gravity, om;
    value_type xmax, ym, x0, y0, dx, dy, x_edge, unifc, unifh,
               dfloe_pack_init, dfloe_pack_thresh, amin, amax;
    value_type rhowtr, rhoice, poisson, dmin, xi, fragility,
               young, visc_rp, kice, kwtr, int_adm, modT, argR, argT, rhoi, rho, rhow;
    value_type fmin, fmax, df, epsc, sigma_c, vbf, vb, flex_rig_coeff;
    value_type dt,duration;

    int nwavedirn, nwavefreq, advdim, ncs ,nt;
    bool ref_Hs_ice, atten, useicevel, steady, breaking, dumpDiag;
    bool docoupling;
    std::string scatmod, advopt, fsdopt;
    std::string wim_gridfile;
    std::vector<value_type> wavedir, wt_simp, wt_om, freq_vec, vec_period, wlng, ag, ap;

    array2_type steady_mask, wave_mask, ice_mask, wtr_mask,
                icec, iceh, swh_in_array,mwp_in_array,mwd_in_array,
                dave, atten_dim, damp_dim, ag2d_eff_temp, mwd, Hs, Tp;
    array3_type ag_eff, ap_eff, wlng_ice, atten_nond, damping, disp_ratio, sdf3d_dir_temp;
    array4_type sdf_dir, sdf_inc;

    array2_type S_freq, taux_om, tauy_om;
    array2_type hp;
    array2_type Fdmax, Ftaux, Ftauy, Fhs, Ftp;

    //array2_type dfloe, nfloes, tau_x, tau_y;
    std::vector<value_type> dfloe, nfloes, tau_x, tau_y;
    std::vector<value_type> mesh_x, mesh_y, mesh_conc, mesh_thick, mesh_dfloe, broken;

    WimGrid wim_grid;

    boost::mpi::timer chrono;
    std::string init_time_str;
    int cpt;

};

} // namespace Wim

#endif
