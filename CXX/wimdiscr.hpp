/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

#ifndef __WIMDISCR_H
#define __WIMDISCR_H 1

#include <iostream>
#include <fstream>
#include <math.h>
//#include <algorithm>

#include <boost/multi_array.hpp>
#include <boost/array.hpp>
#include "boost/program_options.hpp"

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>

#include <boost/any.hpp>
#include <boost/strong_typedef.hpp>
//#include <boost/program_options.hpp>

#define PI M_PI

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

    WimDiscr (po::variables_map const& vmin)
        :
        vm(vmin)
    {}

	void readFile (std::string const& filein, size_type const& nxin, size_type const& nyin);
    void writeFile (std::string const& fileout) const; //, size_type const& nxin, size_type const& nyin) const;
    void wimInit();//(bool steady);
    void wimStep();//(std::string const& scatmod, int const& advdim, std::string const& advopt, bool breaking, bool steady);

	//std::fstream __out( M_filename.c_str(), std::ios::out | std::ios::binary );

    std::vector<std::vector<value_type> > getX() const { return X_array; }
    std::vector<std::vector<value_type> > getY() const { return Y_array; }
    std::vector<std::vector<value_type> > getSCUY() const { return SCUY_array; }
    std::vector<std::vector<value_type> > getSCVX() const { return SCVX_array; }
    std::vector<std::vector<value_type> > getSCP2() const { return SCP2_array; }
    std::vector<std::vector<value_type> > getSCP2I() const { return SCP2I_array; }
    std::vector<std::vector<value_type> > getLANDMASK() const { return LANDMASK_array; }

private:

    po::variables_map const& vm;
    size_type nx, ny; //recno;
    std::vector<std::vector<value_type> > X_array, Y_array, SCUY_array, SCVX_array, SCP2_array, SCP2I_array, LANDMASK_array;

    value_type cfl, dom, guess, Hs_inc, Tp_inc, Tmin, Tmax, gravity, om;
    value_type xm, ym, dx, dy;

    int nwavedirn, nwavefreq, advdim;
    bool refhsice, atten, icevel, steady, breaking;
    std::string scatmod, advopt;
    std::vector<value_type> wavedir, wt_simp, wt_om, freq_vec, vec_period, wlng, ag, ap;

    array2_type wave_mask2, wave_mask, ice_mask, wtr_mask;
    array3_type ag_eff, ap_eff, wlng_ice, atten_nond, damping, disp_ratio;
    array4_type sdf_dir, sdf_inc;

};

template<typename T>
void
WimDiscr<T>::readFile (std::string const& filein, size_type const& nxin, size_type const& nyin)
{

    nx = nxin;
    ny = nyin;

    for (int i = 0; i < nx; i++)
    {
        X_array.push_back(std::vector<value_type>(ny,0));
        Y_array.push_back(std::vector<value_type>(ny,0));
        SCUY_array.push_back(std::vector<value_type>(ny,0));
        SCVX_array.push_back(std::vector<value_type>(ny,0));
        SCP2_array.push_back(std::vector<value_type>(ny,0));
        SCP2I_array.push_back(std::vector<value_type>(ny,0));
        LANDMASK_array.push_back(std::vector<value_type>(ny,0));
    }

    //std::fstream _in( filein, std::ios::binary | std::ios::in | ios::ate );
    std::fstream in(filein, std::ios::binary | std::ios::in);

    if (in.is_open())
    {
        // size_type size = in.tellg();
        // std::cout<<"number of bytes of the input file: "<< size <<"\n";

        //in.seekg(2400, ios::beg);
        // in.seekg(0, ios::beg);
        //in.read((char *)&X_array, sizeof(X_array));


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

        // for (int i = 0; i < nx; i++)
        //     for (int j = 0; j < ny; j++)
        //         std::cout << "Element[" << i << "," << j << "]= " << X_array[i][j] << std::endl;
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
WimDiscr<T>::writeFile (std::string const& fileout) const //, size_type const& nx, size_type const& ny) const
{
    std::fstream out(fileout, std::ios::binary | std::ios::out | std::ios::trunc);

    if (out.is_open())
    {
        size_type size = out.tellg();
        std::cout<<"number of bytes of the input file: "<< size <<"\n";

        for (int i = 0; i < X_array.size(); i++)
            for (int j = 0; j < X_array[i].size(); j++)
                out.write((char *)&X_array[i][j], sizeof(int));


        size = out.tellp();
        std::cout<<"number of bytes of the input file: "<< size <<"\n";


        for (int i = 0; i < Y_array.size(); i++)
            for (int j = 0; j < Y_array[i].size(); j++)
                out.write((char *)&Y_array[i][j], sizeof(int));


        size = out.tellp();
        std::cout<<"number of bytes of the input file: "<< size <<"\n";


        for (int i = 0; i < SCUY_array.size(); i++)
            for (int j = 0; j < SCUY_array[i].size(); j++)
                out.write((char *)&SCUY_array[i][j], sizeof(int));


        size = out.tellp();
        std::cout<<"number of bytes of the input file: "<< size <<"\n";


        for (int i = 0; i < SCVX_array.size(); i++)
            for (int j = 0; j < SCVX_array[i].size(); j++)
                out.write((char *)&SCVX_array[i][j], sizeof(int));


        size = out.tellp();
        std::cout<<"number of bytes of the input file: "<< size <<"\n";


        for (int i = 0; i < SCP2_array.size(); i++)
            for (int j = 0; j < SCP2_array[i].size(); j++)
                out.write((char *)&SCP2_array[i][j], sizeof(int));


        size = out.tellp();
        std::cout<<"number of bytes of the input file: "<< size <<"\n";


        for (int i = 0; i < SCP2I_array.size(); i++)
            for (int j = 0; j < SCP2I_array[i].size(); j++)
                out.write((char *)&SCP2I_array[i][j], sizeof(int));


        size = out.tellp();
        std::cout<<"number of bytes of the input file: "<< size <<"\n";


        for (int i = 0; i < LANDMASK_array.size(); i++)
            for (int j = 0; j < LANDMASK_array[i].size(); j++)
                out.write((char *)&LANDMASK_array[i][j], sizeof(int));


        size = out.tellp();
        std::cout<<"number of bytes of the input file: "<< size <<"\n";

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
WimDiscr<T>::wimInit() //std::string const& scatmod, int const& advdim, std::string const& advopt, bool breaking, bool steady)
{

    nx = vm["nx"].template as<int>();
    ny = vm["ny"].template as<int>();

    nwavedirn = vm["nwavedirn"].template as<int>();
    nwavefreq = vm["nwavefreq"].template as<int>();
    advdim = vm["advdim"].template as<int>();

    refhsice = vm["refhsice"].template as<bool>();
    atten = vm["atten"].template as<bool>();
    icevel = vm["icevel"].template as<bool>();
    steady = vm["steady"].template as<bool>();
    breaking = vm["breaking"].template as<bool>();

    cfl = vm["cfl"].template as<double>();

    scatmod = vm["scatmod"].template as<std::string>();
    advopt = vm["advopt"].template as<std::string>();

    dx = vm["dx"].template as<double>();
    dy = vm["dy"].template as<double>();

    wave_mask2.resize(boost::extents[nx][ny]);

    sdf_dir.resize(boost::extents[nx][ny][nwavedirn][nwavefreq]);
    sdf_inc.resize(boost::extents[nx][ny][nwavedirn][nwavefreq]);

    wt_simp.resize(nwavefreq);
    wt_om.resize(nwavefreq);

    freq_vec.resize(nwavefreq);
    vec_period.resize(nwavefreq);
    wlng.resize(nwavefreq);
    ag.resize(nwavefreq);
    ap.resize(nwavefreq);


    ag_eff.resize(boost::extents[nx][ny][nwavefreq]);
    ap_eff.resize(boost::extents[nx][ny][nwavefreq]);
    wlng_ice.resize(boost::extents[nx][ny][nwavefreq]);

    atten_nond.resize(boost::extents[nx][ny][nwavefreq]);
    damping.resize(boost::extents[nx][ny][nwavefreq]);
    disp_ratio.resize(boost::extents[nx][ny][nwavefreq]);

    // std::cout<<"wave_mask2.size1= "<< wave_mask2.shape()[0] <<"\n";
    // std::cout<<"wave_mask2.size2= "<< wave_mask2.shape()[1] <<"\n";
    //std::fill( wave_mask2.data(), wave_mask2.data() + wave_mask2.num_elements(), 10 );

    // wave_mask2
    if (steady)
        std::fill( &wave_mask2[0][0], &wave_mask2[3][0], 1. );

    // wavedir
    if (nwavedirn==1)
    {
        wavedir.push_back(-90.);
    }
    else
    {
        value_type theta_max = 90.;
        value_type theta_min = -270.;
        value_type dtheta = (theta_min-theta_max)/nwavedirn;

        for (int i = 0; i < nwavedirn; i++)
        {
            wavedir.push_back(theta_max+i*dtheta);
            std::cout<<"wavedir["<< i <<"]= "<< theta_max+i*dtheta <<"\n";
        }
    }


    // freq_vec
    Hs_inc = 2.0;
    Tp_inc = 12.0;
    Tmin = 2.5;
    Tmax = 25.;
    gravity = 9.81;

    value_type fmin,fmax,df;

    if (nwavefreq == 1)
    {
        freq_vec[0] = 1./Tp_inc;
    }
    else
    {
        fmin = 1./Tmax;
        fmax = 1./Tmin;
        df = (fmax-fmin)/(nwavefreq-1.0);

        for (int i = 0; i < nwavefreq; i++)
        {
            freq_vec[i] = fmin+i*df;
        }
    }

    // wt_simp
    if (nwavefreq >1)
    {
        std::fill(wt_simp.begin(), wt_simp.end(), 2.);
        wt_simp[0] = wt_simp[nwavefreq-1] = 1.;

        size_type w = 1;
        while (w < nwavefreq-1)
        {
            wt_simp[w] = 4.;
            w +=2;
        }

        dom   = 2*PI*(freq_vec[nwavefreq-1]-freq_vec[1])/(nwavefreq-1.0);
        wt_om = wt_simp;

        std::for_each(wt_om.begin(), wt_om.end(), [&](value_type& f){ f*=2; });
    }
    else
    {
        wt_om[0] = 1.;
    }


    // vec_period
    vec_period = freq_vec;
    std::for_each(vec_period.begin(), vec_period.end(), [&](value_type& f){ f =1./f; });

    // wlng
    wlng = freq_vec;
    std::for_each(wlng.begin(), wlng.end(), [&](value_type& f){ f = gravity/(2*PI*std::pow(f,2.)); });

    // ap
    ap = wlng;
    std::for_each(ap.begin(), ap.end(), [&](value_type& f){ f = std::sqrt(gravity*f/(2*PI)); });

    // ag
    ag = ap;
    std::for_each(ag.begin(), ag.end(), [&](value_type& f){ f /=2 ; });



    xm = X_array[0][0]+(nx-1)*dx;

    std::cout<<"nx= "<< nx << " and dx= "<< dx << " and X_array[0][0]= "<< X_array[0][0] <<"\n";

    ym = Y_array[0][0]+(ny-1)*dy;

    std::cout<<"XM= "<< xm <<"\n";
    std::cout<<"YM= "<< ym <<"\n";

    om = 2*PI*freq_vec[0];


    // for (size_type k = 0; k < nwavedirn; k++)
    // {
    //     for (size_type i = 0; i < nx; i++)
    //     {
    //         for (size_type j = 0; j < ny; j++)
    //         {

    //         }
    //     }
    // }


    // for (value_type const& it : freq_vec)
    // {
    //     std::cout<<"[F]= "<< it <<"\n";
    // }

    // for (value_type const& it : vec_period)
    // {
    //     std::cout<<"[T]= "<< it <<"\n";
    // }

    // for (value_type const& it : ap)
    // {
    //     std::cout<<"[AP]= "<< it <<"\n";
    // }

    // for (value_type const& it : ag)
    // {
    //     std::cout<<"[AG]= "<< it <<"\n";
    // }

    // for (int i = 0; i < wt_simp.size(); i++)
    // {
    //     std::cout<<"wt_simp["<< i <<"]= "<< wt_simp[i] <<"\n";
    // }
}


template<typename T>
void
WimDiscr<T>::wimStep() //std::string const& scatmod, int const& advdim, std::string const& advopt, bool breaking, bool steady)
{
    // std::vector<std::vector<value_type> > mom0, mom2, var_strain, S_freq, atten_dim, damp_dim, taux_om, tauy_om, tmp1, mom0w, mom2w;
    // value_type adv_dir, E_tot, sig_strain, Pstrain, P_crit, wlng_crest, Dc, om, ommin, ommax, om1, lam1, lam2, dom, F, Dave, c1d, kice, tmp;

    value_type adv_dir;

    if (vm["steady"].template as<bool>())
    {
        for (size_type i = 0; i < nwavefreq; i++)
        {
            for (size_type j = 0; j < nwavedirn; j++)
            {
                adv_dir = (-PI/180)*(wavedir[j]+90.);
                std::cout<<"COS(adv_dir)= "<< std::cos(adv_dir) <<"\n";

                if (std::cos(adv_dir) >= 0)
                {
                    for (size_type k = 0; k < nx; k++)
                    {
                        for (size_type l = 0; l < ny; l++)
                        {
                            if (wave_mask2[k][l] >=0)
                                sdf_dir[k][l][j][i] = sdf_inc[k][l][j][i];
                        }
                    }
                }

            }
        }
    }


}

#endif
