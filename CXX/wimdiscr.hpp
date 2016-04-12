/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

/**
 * @file   wimdiscr.hpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Mon Aug  3 11:53:19 2015
 */


#ifndef __WIMDISCR_H
#define __WIMDISCR_H 1

#include <iostream>
#include <fstream>
#include <math.h>
#include <algorithm>
// #include <complex>
// #include <cmath>

#include <boost/multi_array.hpp>
#include <boost/array.hpp>
#include <boost/program_options.hpp>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>

#include <boost/any.hpp>
//#include <boost/strong_typedef.hpp>
#include <boost/format.hpp>
#include <iomanip>
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
    void readData(std::string const& filein);
    void writeFile(size_type const& timestp, value_type const& t_out) const;
    void wimInit();
    void wimStep();
    void wimRun();
    void floeScaling(value_type const& dmax, value_type& dave);
    void advAttenSimple(array3_type& Sdir, array2_type& Sfreq,array2_type& taux_omega,array2_type& tauy_omega, array2_type const& ag2d_eff);
    void advAttenIsotropic(array3_type& Sdir, array2_type& Sfreq,array2_type& taux_omega,array2_type& tauy_omega, array2_type const& ag2d_eff);
    void waveAdvWeno(array2_type& h, array2_type const& u, array2_type const& v);
    void weno3pdV2(array2_type const& gin, array2_type const& u, array2_type const& v, array2_type const& scuy,
                   array2_type const& scvx, array2_type const& scp2i, array2_type const& scp2, array2_type& saoout);

    void padVar(array2_type const& u, array2_type& upad);
    //void padVar2d(array2_type const& u, array2_type& upad);
    void calcMWD();

    array2_type getX() const { return X_array; }
    array2_type getY() const { return Y_array; }
    array2_type getSCUY() const { return SCUY_array; }
    array2_type getSCVX() const { return SCVX_array; }
    array2_type getSCP2() const { return SCP2_array; }
    array2_type getSCP2I() const { return SCP2I_array; }
    array2_type getLANDMASK() const { return LANDMASK_array; }


private:

    po::variables_map const& vm;
    int nx, ny, nxext, nyext, nbdy, nbdx;
    array2_type X_array, Y_array, SCUY_array, SCVX_array, SCP2_array, SCP2I_array, LANDMASK_array;

    value_type cfl, dom, guess, Hs_inc, Tp_inc, mwd_inc, Tmin, Tmax, gravity, om;
    value_type xmax, ym, x0, y0, dx, dy, x_edge, unifc, unifh, dfloe_pack_init, amin, amax, dt;
    value_type rhowtr, rhoice, poisson, dmin, xi, fragility, young, visc_rp, kice, kwtr, int_adm, modT, argR, argT, rhoi, rho, rhow;
    value_type fmin, fmax, df, epsc, sigma_c, vbf, vb, flex_rig_coeff;

    int nwavedirn, nwavefreq, advdim, ncs;
    bool ref_Hs_ice, atten, icevel, steady, breaking;
    std::string scatmod, advopt;
    std::vector<value_type> wavedir, wt_simp, wt_om, freq_vec, vec_period, wlng, ag, ap;

    array2_type wave_mask2, wave_mask, ice_mask, wtr_mask, icec, iceh, dfloe, atten_dim, damp_dim, ag2d_eff_temp, tau_x, tau_y, mwd, Hs, Tp;
    array3_type ag_eff, ap_eff, wlng_ice, atten_nond, damping, disp_ratio, sdf3d_dir_temp;
    array4_type sdf_dir, sdf_inc;

    array2_type S_freq, taux_om, tauy_om;
    // array2_type tmp1, mom0, mom2, mom0w, mom2w, var_strain;
    // array2_type uwave, vwave, temp;

    array2_type hp;

    // array2_type ful, fuh, fvl, fvh, gt, sao;

    // array2_type u_pad, v_pad, scp2_pad, scp2i_pad, scuy_pad, scvx_pad, h_pad;

    // variables for calcMWD
    //array2_type cmom0, cmom_dir, CSfreq, cmom_dir0, CF;

    array2_type Fdmax, Ftaux, Ftauy, Fhs, Ftp;

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


    // std::cout<<"x0= "<< x0 <<"\n";
    // std::cout<<"y0= "<< y0 <<"\n";

    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            X_array[i][j] = x0 + i*dx+.5*dx;
            Y_array[i][j] = y0 + j*dy+.5*dy;
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

    dx = vm["dx"].template as<double>();
    dy = vm["dy"].template as<double>();


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
WimDiscr<T>::writeFile (size_type const& timestp, value_type const& t_out) const
{

    std::string str = vm["outparentdir"].template as<std::string>();

    char * senv = ::getenv( "WIM2D_PATH" );
    if ( (str == ".") && (senv != NULL) && (senv[0] != '\0') )
    {
        str = std::string( senv ) + "/CXX";
    }

    fs::path path(str);
    path /= "outputs/binaries/prog";

    if ( !fs::exists(path) )
        fs::create_directories(path);



    std::string timestpstr;

    if (timestp < 10)
        timestpstr = "00"+ std::to_string(timestp);
    else if (timestp < 100)
        timestpstr = "0"+ std::to_string(timestp);
    else
        timestpstr = std::to_string(timestp);

    //std::cout<<"TIMESTEP= "<< timestpstr <<"\n";

    //std::string fileout = (boost::format( "%1%/wim_prog%2%.a" ) % path.string() % timestp).str();
    std::string fileout = (boost::format( "%1%/wim_prog%2%.a" ) % path.string() % timestpstr).str();
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

    // export the txt file for grid field information

    //std::string fileoutb = (boost::format( "%1%/wim_prog%2%.b" ) % path.string() % timestp).str();
    std::string fileoutb = (boost::format( "%1%/wim_prog%2%.b" ) % path.string() % timestpstr).str();
    std::fstream outb(fileoutb, std::ios::out | std::ios::trunc);
    //std::ofstream outb(fileoutb, std::ios::out | std::ios::trunc);

    if (outb.is_open())
    {
        outb << std::setw(15) << std::left << 07    << "    Nrecs    # "<< "Number of records" <<"\n";
        outb << std::setw(15) << std::left << 0     << "    Norder   # "<< "Storage order [column-major (F/matlab) = 1; row-major (C) = 0]" <<"\n";
        outb << std::setw(15) << std::left << nx    << "    nx       # "<< "Record length in x direction (elements)" <<"\n";
        outb << std::setw(15) << std::left << ny    << "    ny       # "<< "Record length in y direction (elements)" <<"\n";
        outb << std::setw(15) << std::left << t_out << "    t_out    # "<< "Model time of output (s)" <<"\n";
        //outb << std::setw(15) << std::left << nwavefreq << "          "<< "Number of wave frequencies" <<"\n";
        //outb << std::setw(15) << std::left << nwavedirn << "          "<< "Number of wave directions" <<"\n";

        outb <<"\n";

        outb << "Record number and name:" <<"\n";
        outb << std::setw(15) << std::left << 01 << "          "<< "icec" <<"\n";
        outb << std::setw(15) << std::left << 02 << "          "<< "iceh" <<"\n";
        outb << std::setw(15) << std::left << 03 << "          "<< "Dmax" <<"\n";
        outb << std::setw(15) << std::left << 04 << "          "<< "tau_x" <<"\n";
        outb << std::setw(15) << std::left << 05 << "          "<< "tau_y" <<"\n";
        outb << std::setw(15) << std::left << 06 << "          "<< "Hs" <<"\n";
        outb << std::setw(15) << std::left << 07 << "          "<< "Tp" <<"\n";
    }
    else
    {
        std::cout << "Cannot open " << fileoutb  << "\n";
        std::cerr << "error: open file " << fileoutb << " for output failed!" <<"\n";
        std::abort();
    }


}

template<typename T>
void
WimDiscr<T>::readData (std::string const& filein)
{

    //array2_type X_array(boost::extents[nx][ny],boost::fortran_storage_order());
    // X_array.resize(boost::extents[nx][ny]);
    // Y_array.resize(boost::extents[nx][ny]);
    // SCUY_array.resize(boost::extents[nx][ny]);
    // SCVX_array.resize(boost::extents[nx][ny]);
    // SCP2_array.resize(boost::extents[nx][ny]);
    // SCP2I_array.resize(boost::extents[nx][ny]);
    // LANDMASK_array.resize(boost::extents[nx][ny]);

    Fdmax.resize(boost::extents[nx][ny]);
    Ftaux.resize(boost::extents[nx][ny]);
    Ftauy.resize(boost::extents[nx][ny]);
    Fhs.resize(boost::extents[nx][ny]);
    Ftp.resize(boost::extents[nx][ny]);

    // dx = vm["dx"].template as<double>();
    // dy = vm["dy"].template as<double>();

    char * senv = ::getenv( "WIM2D_PATH" );

    std::string str = std::string( senv ) + "/fortran/run/out/binaries/prog";


    fs::path path(str);

    std::string _filein = (boost::format("%1%/%2%") % path.string() % filein).str();

    // s::path path(str);
    // path /= "outputs/binaries/prog";

    //std::fstream in(filein, std::ios::binary | std::ios::in);
    std::fstream in(_filein, std::ios::binary | std::ios::in);

    if (in.is_open())
    {
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                in.read((char *)&Fdmax[i][j], sizeof(int));

        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                in.read((char *)&Ftaux[i][j], sizeof(int));

        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                in.read((char *)&Ftauy[i][j], sizeof(int));

        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                in.read((char *)&Fhs[i][j], sizeof(int));

        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                in.read((char *)&Ftp[i][j], sizeof(int));

        in.close();
    }
    else
    {
        std::cout << "Cannot open " << _filein  << "\n";
        std::cerr << "error: open file " << _filein << " for input failed!" <<"\n";
        std::abort();
    }
}



} // namespace WIM2D

#endif
