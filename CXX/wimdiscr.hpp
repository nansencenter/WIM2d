/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

#ifndef __WIMDISCR_H
#define __WIMDISCR_H 1

#include <iostream>
#include <fstream>
#include <math.h>
#include <algorithm>

#include <boost/multi_array.hpp>
#include <boost/array.hpp>
#include <boost/program_options.hpp>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>

#include <boost/any.hpp>
#include <boost/strong_typedef.hpp>
#include <boost/format.hpp>

#ifdef __cplusplus
extern "C"
{
#endif
#include<RTparam_outer.h>
#ifdef __cplusplus
}
#endif


#define PI M_PI

namespace po = boost::program_options;
namespace fs = boost::filesystem;

namespace WIM2D
{
template<typename T=float> class WimDiscr
{

	typedef T value_type;
    typedef size_t size_type;
	typedef boost::multi_array<value_type, 2> array2_type;
    typedef boost::multi_array<value_type, 3> array3_type;
    typedef boost::multi_array<value_type, 4> array4_type;
    typedef typename array2_type::index index;

public:

    WimDiscr (po::variables_map const& vmin)
        :
        vm(vmin),
        nx(vm["nx"].template as<int>()),
        ny(vm["ny"].template as<int>())
    {}

    void wimGrid();
	void readFile(std::string const& filein);
    void writeFile(size_type const& timestp) const;
    void wimInit();
    void wimStep();
    void wimRun();
    void floeScaling(value_type const& dmax, value_type& dave);
    void advAttenSimple(array3_type& Sdir, array2_type& Sfreq,array2_type& taux_om,array2_type& tauy_om, array2_type& ag2d_eff);
    void advAttenIsotropic(array3_type& Sdir, array2_type& Sfreq,array2_type& taux_om,array2_type& tauy_om, array2_type& ag2d_eff);
    void waveAdvWeno(array2_type& h, array2_type const& u, array2_type const& v);
    void weno3pdV2(array2_type const& gin, array2_type const& u, array2_type const& v, array2_type const& scuy,
                   array2_type const& scvx, array2_type const& scp2i, array2_type const& scp2, array2_type& saoout);
    void padVar(array2_type const& u, array2_type& upad);

    array2_type getX() const { return X_array; }
    array2_type getY() const { return Y_array; }
    array2_type getSCUY() const { return SCUY_array; }
    array2_type getSCVX() const { return SCVX_array; }
    array2_type getSCP2() const { return SCP2_array; }
    array2_type getSCP2I() const { return SCP2I_array; }
    array2_type getLANDMASK() const { return LANDMASK_array; }



private:

    po::variables_map const& vm;
    size_type nx, ny, nxext, nyext, nbdy;
    array2_type X_array, Y_array, SCUY_array, SCVX_array, SCP2_array, SCP2I_array, LANDMASK_array;

    value_type cfl, dom, guess, Hs_inc, Tp_inc, mwd_inc, Tmin, Tmax, gravity, om;
    value_type xmax, ym, x0, y0, dx, dy, x_edge, unifc, unifh, dfloe_pack_init, amin, amax, dt;
    value_type rhowtr, rhoice, poisson, dmin, xi, fragility, young, visc_rp, kice, kwtr, int_adm, modT, argR, argT, rhoi, rho, rhow;
    value_type fmin, fmax, df;

    int nwavedirn, nwavefreq, advdim, ncs;
    bool ref_Hs_ice, atten, icevel, steady, breaking;
    std::string scatmod, advopt;
    std::vector<value_type> wavedir, wt_simp, wt_om, freq_vec, vec_period, wlng, ag, ap;

    array2_type wave_mask2, wave_mask, ice_mask, wtr_mask, icec, iceh, dfloe, atten_dim, damp_dim, ag2d_eff_temp, tau_x, tau_y, mwd, Hs, Tp;
    array3_type ag_eff, ap_eff, wlng_ice, atten_nond, damping, disp_ratio, sdf3d_dir_temp;
    array4_type sdf_dir, sdf_inc;


    array2_type S_freq, taux_om, tauy_om, tmp1, mom0, mom2, mom0w, mom2w, var_strain;
    array2_type uwave, vwave, temp;

    array2_type hp;

    array2_type ful, fuh, fvl, fvh, gt, sao;

    array2_type u_pad, v_pad, scp2_pad, scp2i_pad, scuy_pad, scvx_pad, h_pad;

};


template<typename T>
void
WimDiscr<T>::wimGrid()
{

    std::string str = vm["outparentdir"].template as<std::string>();

    char * senv = ::getenv( "WIM2D_PATH" );
    if ( (str == ".") && (senv != NULL) && (senv[0] != '\0') )
    {
        str = std::string( senv ) + "/CXX";
    }

    fs::path path(str);
    path /= "outputs/binaries";

    if ( !fs::exists(path) )
        fs::create_directories(path);


    std::string fileout = (boost::format( "%1%/wim_grid.a" ) % path.string()).str();
    std::fstream out(fileout, std::ios::binary | std::ios::out | std::ios::trunc);

    X_array.resize(boost::extents[nx][ny]);
    Y_array.resize(boost::extents[nx][ny]);
    SCUY_array.resize(boost::extents[nx][ny]);
    SCVX_array.resize(boost::extents[nx][ny]);
    SCP2_array.resize(boost::extents[nx][ny]);
    SCP2I_array.resize(boost::extents[nx][ny]);
    LANDMASK_array.resize(boost::extents[nx][ny]);

    dx = vm["dx"].template as<double>();
    dy = vm["dy"].template as<double>();

    x0 = vm["xmin"].template as<double>();
    y0 = vm["ymin"].template as<double>();


    std::cout<<"x0= "<< x0 <<"\n";
    std::cout<<"y0= "<< y0 <<"\n";

    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            X_array[i][j] = x0 + i*dx;
            Y_array[i][j] = y0 + j*dy;
            SCUY_array[i][j] = dy;
            SCVX_array[i][j] = dx;
            SCP2_array[i][j] = dx*dy;
            SCP2I_array[i][j] = 1./(dx*dy);
        }
    }

    if (out.is_open())
    {
        for (int i = 0; i < X_array.shape()[0]; i++)
            for (int j = 0; j < X_array.shape()[1]; j++)
                out.write((char *)&X_array[i][j], sizeof(int));

        for (int i = 0; i < Y_array.shape()[0]; i++)
            for (int j = 0; j < Y_array.shape()[1]; j++)
                out.write((char *)&Y_array[i][j], sizeof(int));

        for (int i = 0; i < SCUY_array.shape()[0]; i++)
            for (int j = 0; j < SCUY_array.shape()[1]; j++)
                out.write((char *)&SCUY_array[i][j], sizeof(int));

        for (int i = 0; i < SCVX_array.shape()[0]; i++)
            for (int j = 0; j < SCVX_array.shape()[1]; j++)
                out.write((char *)&SCVX_array[i][j], sizeof(int));

        for (int i = 0; i < SCP2_array.shape()[0]; i++)
            for (int j = 0; j < SCP2_array.shape()[1]; j++)
                out.write((char *)&SCP2_array[i][j], sizeof(int));

        for (int i = 0; i < SCP2I_array.shape()[0]; i++)
            for (int j = 0; j < SCP2I_array.shape()[1]; j++)
                out.write((char *)&SCP2I_array[i][j], sizeof(int));

        for (int i = 0; i < LANDMASK_array.shape()[0]; i++)
            for (int j = 0; j < LANDMASK_array.shape()[1]; j++)
                out.write((char *)&LANDMASK_array[i][j], sizeof(int));

        out.close();
    }
    else
    {
        std::cout << "Cannot open " << fileout  << "\n";
        std::cerr << "error: open file " << fileout << " for output failed!" <<"\n";
        std::abort();
    }

}

template<typename T>
void
WimDiscr<T>::readFile (std::string const& filein)
{

    //array2_type X_array(boost::extents[nx][ny],boost::fortran_storage_order());
    X_array.resize(boost::extents[nx][ny]);
    Y_array.resize(boost::extents[nx][ny]);
    SCUY_array.resize(boost::extents[nx][ny]);
    SCVX_array.resize(boost::extents[nx][ny]);
    SCP2_array.resize(boost::extents[nx][ny]);
    SCP2I_array.resize(boost::extents[nx][ny]);
    LANDMASK_array.resize(boost::extents[nx][ny]);

    std::fstream in(filein, std::ios::binary | std::ios::in);

    if (in.is_open())
    {
        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                in.read((char *)&X_array[i][j], sizeof(int));

        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                in.read((char *)&Y_array[i][j], sizeof(int));


        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                in.read((char *)&SCUY_array[i][j], sizeof(int));


        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                in.read((char *)&SCVX_array[i][j], sizeof(int));


        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                in.read((char *)&SCP2_array[i][j], sizeof(int));


        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                in.read((char *)&SCP2I_array[i][j], sizeof(int));


        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                in.read((char *)&LANDMASK_array[i][j], sizeof(int));

        in.close();
    }
    else
    {
        std::cout << "Cannot open " << filein  << "\n";
        std::cerr << "error: open file " << filein << " for input failed!" <<"\n";
        std::abort();
    }
}


template<typename T>
void
WimDiscr<T>::writeFile (size_type const& timestp) const
{

    std::string str = vm["outparentdir"].template as<std::string>();

    char * senv = ::getenv( "WIM2D_PATH" );
    if ( (str == ".") && (senv != NULL) && (senv[0] != '\0') )
    {
        str = std::string( senv ) + "/CXX";
    }

    fs::path path(str);
    path /= "outputs/binaries";

    if ( !fs::exists(path) )
        fs::create_directories(path);


    std::string fileout = (boost::format( "%1%/wim_prog%2%.a" ) % path.string() % timestp).str();
    std::fstream out(fileout, std::ios::binary | std::ios::out | std::ios::trunc);

    if (out.is_open())
    {
        for (int i = 0; i < icec.shape()[0]; i++)
            for (int j = 0; j < icec.shape()[1]; j++)
                out.write((char *)&icec[i][j], sizeof(int));

        for (int i = 0; i < iceh.shape()[0]; i++)
            for (int j = 0; j < iceh.shape()[1]; j++)
                out.write((char *)&iceh[i][j], sizeof(int));

        for (int i = 0; i < dfloe.shape()[0]; i++)
            for (int j = 0; j < dfloe.shape()[1]; j++)
                out.write((char *)&dfloe[i][j], sizeof(int));

        for (int i = 0; i < tau_x.shape()[0]; i++)
            for (int j = 0; j < tau_x.shape()[1]; j++)
                out.write((char *)&tau_x[i][j], sizeof(int));

        for (int i = 0; i < tau_y.shape()[0]; i++)
            for (int j = 0; j < tau_y.shape()[1]; j++)
                out.write((char *)&tau_y[i][j], sizeof(int));

        for (int i = 0; i < Hs.shape()[0]; i++)
            for (int j = 0; j < Hs.shape()[1]; j++)
                out.write((char *)&Hs[i][j], sizeof(int));

        for (int i = 0; i < Tp.shape()[0]; i++)
            for (int j = 0; j < Tp.shape()[1]; j++)
                out.write((char *)&Tp[i][j], sizeof(int));

        out.close();
    }
    else
    {
        std::cout << "Cannot open " << fileout  << "\n";
        std::cerr << "error: open file " << fileout << " for output failed!" <<"\n";
        std::abort();
    }

}

} // namespace WIM2D

#endif
