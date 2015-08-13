/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

/**
 * @file   wimdiscr.cpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Mon Aug  3 11:52:35 2015
 */


#include <wimdiscr.hpp>

namespace WIM2D
{
template<typename T>
void
WimDiscr<T>::wimInit()
{
	// wim2d grid generation
	wimGrid();
    //readFile("wim_grid.a");

	nwavedirn = vm["nwavedirn"].template as<int>();
    nwavefreq = vm["nwavefreq"].template as<int>();
    advdim = vm["advdim"].template as<int>();

    ref_Hs_ice = vm["refhsice"].template as<bool>();
    atten = vm["atten"].template as<bool>();
    icevel = vm["icevel"].template as<bool>();
    steady = vm["steady"].template as<bool>();
    breaking = vm["breaking"].template as<bool>();

    scatmod = vm["scatmod"].template as<std::string>();
    advopt = vm["advopt"].template as<std::string>();

    cfl = vm["cfl"].template as<double>();

    nbdx = vm["nbdx"].template as<int>();
    nbdy = vm["nbdy"].template as<int>();

    if (advdim == 1)
        nbdy = 0;

    nxext = nx+2*nbdx;
    nyext = ny+2*nbdy;

    wave_mask2.resize(boost::extents[nx][ny]);
    wave_mask.resize(boost::extents[nx][ny]);

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

    ag2d_eff_temp.resize(boost::extents[nx][ny]);
    sdf3d_dir_temp.resize(boost::extents[nx][ny][nwavedirn]);

    atten_nond.resize(boost::extents[nx][ny][nwavefreq]);
    damping.resize(boost::extents[nx][ny][nwavefreq]);
    disp_ratio.resize(boost::extents[nx][ny][nwavefreq]);


    ice_mask.resize(boost::extents[nx][ny]);
    wtr_mask.resize(boost::extents[nx][ny]);


    icec.resize(boost::extents[nx][ny]);
    iceh.resize(boost::extents[nx][ny]);
    dfloe.resize(boost::extents[nx][ny]);

    atten_dim.resize(boost::extents[nx][ny]);
    damp_dim.resize(boost::extents[nx][ny]);


    tau_x.resize(boost::extents[nx][ny]);
    tau_y.resize(boost::extents[nx][ny]);

    Hs.resize(boost::extents[nx][ny]);
    Tp.resize(boost::extents[nx][ny]);
    mwd.resize(boost::extents[nx][ny]);

    //std::fill( wave_mask2.data(), wave_mask2.data() + wave_mask2.num_elements(), 10 );
    // wave_mask2
    if (steady)
        std::fill( &wave_mask2[0][0], &wave_mask2[3][0], 1. );

    // parameters
    Hs_inc = 2.0;
    Tp_inc = 12.0;
    mwd_inc = -90.;//-135.;//-90.;
    Tmin = 2.5;
    Tmax = 25.;
    gravity = 9.81;

    unifc = 0.7;
    unifh = 2.0;
    dfloe_pack_init = 300.0;

    rhowtr = 1025.;
    rhoice = 922.5;
    poisson = 0.3;
    dmin = 20.;
    xi = 2.;
    fragility = 0.9;

    young = vm["young"].template as<double>();
    visc_rp = vm["viscrp"].template as<double>();

    std::fill( disp_ratio.data(), disp_ratio.data() + disp_ratio.num_elements(), 1. );

    // freq_vec
    if (nwavefreq == 1)
    {
        freq_vec[0] = 1./Tp_inc;
    }
    else
    {
        // multiple frequencies
        fmin = 1./Tmax;
        fmax = 1./Tmin;
        df = (fmax-fmin)/(nwavefreq-1.0);

        for (int i = 0; i < nwavefreq; i++)
        {
            freq_vec[i] = fmin+i*df;
        }
    }

    // set directional grid
    // wavedir
    if (nwavedirn == 1)
    {
        wavedir.push_back(mwd_inc);
    }
    else
    {
        value_type theta_max = 90.;
        value_type theta_min = -270.;
        value_type dtheta = (theta_min-theta_max)/nwavedirn;

        for (int i = 0; i < nwavedirn; i++)
        {
            wavedir.push_back(theta_max+i*dtheta);
        }
    }


    // wt_simp
    if (nwavefreq >1)
    {
        std::fill(wt_simp.begin(), wt_simp.end(), 2.);
        wt_simp[0] = 1.;
        wt_simp[nwavefreq-1] = 1.;

        int w = 1;
        while (w < nwavefreq-1)
        {
            wt_simp[w] = 4.;
            w +=2;
        }

        dom   = 2*PI*(freq_vec[nwavefreq-1]-freq_vec[1])/(nwavefreq-1.0);
        //wt_om = dom*wt_simp/3.0;
        wt_om = wt_simp;
        std::for_each(wt_om.begin(), wt_om.end(), [&](value_type& f){ f = dom*f/3.0; });
    }
    else
    {
        wt_om[0] = 1.;
    }


    // vec_period
    vec_period = freq_vec;
    std::for_each(vec_period.begin(), vec_period.end(), [&](value_type& f){ f = 1./f; });

    // wlng
    wlng = freq_vec;
    std::for_each(wlng.begin(), wlng.end(), [&](value_type& f){ f = gravity/(2*PI*std::pow(f,2.)); });

    // ap
    ap = wlng;
    std::for_each(ap.begin(), ap.end(), [&](value_type& f){ f = std::sqrt(gravity*f/(2*PI)); });

    // ag
    ag = ap;
    std::for_each(ag.begin(), ag.end(), [&](value_type& f){ f = f/2 ; });


    x0 = X_array[0][0];
    xmax =X_array[nx-1][ny-1]; //x0+(nx-1)*dx;



    y0 = Y_array[0][0];
    ym = Y_array[0][0]+(ny-1)*dy;
    value_type ymax = Y_array[nx-1][ny-1];

    x_edge = 0.5*(x0+xmax)-0.8*(0.5*(xmax-x0));

    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            if (X_array[i][j] < x_edge)
            {
                wave_mask[i][j] = 1.;
                Hs[i][j] = Hs_inc;
                Tp[i][j] = Tp_inc;
                mwd[i][j] = mwd_inc;
            }
        }
    }

    x_edge = 0.5*(x0+xmax)-0.7*(0.5*(xmax-x0));

    std::cout<<"x0= "<< x0 <<"\n";
    std::cout<<"xmax= "<< xmax <<"\n";
    std::cout<<"x_edge= "<< x_edge <<"\n";

    //value_type ymax = y0+(ny-1)*dy;
    std::cout<<"y0= "<< y0 <<"\n";
    std::cout<<"ymax= "<< ymax <<"\n";

    // wtr_mask
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            if (X_array[i][j] < x_edge)
                wtr_mask[i][j] = 1.;

            ice_mask[i][j] = (1-wtr_mask[i][j])*(1-LANDMASK_array[i][j]);

            icec[i][j] = unifc*ice_mask[i][j];
            iceh[i][j] = unifh*ice_mask[i][j];
            dfloe[i][j] = dfloe_pack_init*ice_mask[i][j];

            //std::cout<<"ICE_MASK["<< i  << "," << j << "]= " << ice_mask[i][j] <<"\n";
        }
    }



    om = 2*PI*freq_vec[0];
    std::vector<value_type> Sfreq(nwavefreq);
    std::vector<value_type> theta_fac(nwavedirn,0.);
    value_type f1, f2, f3, t_m, om_m, chi;

    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            if (wave_mask[i][j] == 1.)
            {
                if (nwavefreq > 1)
                {
                    for (int fq = 0; fq < nwavefreq; fq++)
                    {
                        t_m = 2*PI/om;
                        om_m  = 2*PI/Tp[i][j];
                        f1 = 5.0/16.0*std::pow(Hs[i][j],2.)*std::pow(om_m,4.);
                        f2 = 1.0/std::pow(om,5.);
                        f3 = std::exp(-1.25*std::pow(t_m/Tp[i][j],4.));
                        Sfreq[fq] = f1*f2*f3;
                    }
                }
                else
                {
                    Sfreq[0] = std::pow(Hs[i][j]/4.0,2.);
                }
            }

            // direcional spreading
            if (nwavedirn == 1)
            {
                std::fill(theta_fac.begin(), theta_fac.end(), 1.);
            }
            else
            {
                for (int dn = 0; dn < nwavedirn; dn++)
                {
                    chi = PI*(wavedir[dn]-mwd[i][j])/180.0;
                    if (std::cos(chi) > 0.)
                        theta_fac[dn] = 2.0*std::pow(cos(chi),2.)/PI;
                    else
                        theta_fac[dn] = 0.;
                }
            }

            // combine freq and dir
            for (int fq = 0; fq < nwavefreq; fq++)
            {
                for (int dn = 0; dn < nwavedirn; dn++)
                {
                    sdf_inc[i][j][dn][fq] = Sfreq[fq]*theta_fac[dn];

                    // if ((i==1) && (j==1))
                    //     adv_dir = (-PI/180.0)*(wavedir[dn]+90.0);
                }
            }

        }
    }


    double params[5];
    params[0] = young;
    params[1] = gravity;
    params[2] = rhowtr;
    params[3] = rhoice;
    params[4] = poisson;

    double outputs[8];

    for (int fq = 0; fq < nwavefreq; fq++)
    {
        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                if (ice_mask[i][j] == 1.)
                {
                    om = 2*PI*freq_vec[fq];

                    if (fq == 0)
                        guess = std::pow(om,2.);
                    else
                        guess = 2*PI*wlng_ice[i][j][fq-1];


                    RTparam_outer(outputs,iceh[i][j],double(om),double(visc_rp),double(guess),params);

                    damping[i][j][fq] = outputs[0];
                    kice = outputs[1];
                    kwtr = outputs[2];
                    int_adm = outputs[3];
                    atten_nond[i][j][fq] = outputs[4];
                    modT = outputs[5];
                    argR = outputs[6];
                    argT = outputs[7];


                    // for (int i=0; i < 8; ++i)
                    //     std::cout<<"params["<< i <<"]= " << outputs[i] << "\n";


                    disp_ratio[i][j][fq] = (kice*modT)/kwtr;
                    //std::cout<<"RATIO= "<< iceh[i][j] <<"\n";
                    wlng_ice[i][j][fq] = 2*PI/kice;


                    ag_eff[i][j][fq] = ag[fq];
                    ap_eff[i][j][fq] = ap[fq];

                    if (icevel)
                    {
                        wlng_ice[i][j][fq] = wlng[fq];
                    }
                }
                else
                {
                    ag_eff[i][j][fq] = ag[fq];
                    ap_eff[i][j][fq] = ap[fq];
                    wlng_ice[i][j][fq] = wlng[fq];
                }
            }
        }

        if (fq == 0)
        {
            amax = *std::max_element(ag_eff.data(),ag_eff.data() + ag_eff.num_elements());
        }

        if (fq == nwavefreq-1)
        {
            amin = *std::min_element(ag_eff.data(),ag_eff.data() + ag_eff.num_elements());
        }

    }


    if (!atten)
    {
        std::fill( atten_nond.data(), atten_nond.data() + atten_nond.num_elements(), 0. );
        std::fill( damping.data(), damping.data() + damping.num_elements(), 0. );
    }


    dt = cfl*dx/amax;

    // warning: in wimStep()
    // tmp1.resize(boost::extents[nx][ny]);
    // mom0.resize(boost::extents[nx][ny]);
    // mom2.resize(boost::extents[nx][ny]);
    // var_strain.resize(boost::extents[nx][ny]);
    // mom0w.resize(boost::extents[nx][ny]);
    // mom2w.resize(boost::extents[nx][ny]);


    // warning: in wimStep()
    // uwave.resize(boost::extents[nx][ny]);
    // vwave.resize(boost::extents[nx][ny]);
    // temp.resize(boost::extents[nx][ny]);

    S_freq.resize(boost::extents[nx][ny]);
    taux_om.resize(boost::extents[nx][ny]);
    tauy_om.resize(boost::extents[nx][ny]);


    hp.resize(boost::extents[nxext][nyext]);


    // ful.resize(boost::extents[nxext][nyext]);
    // fuh.resize(boost::extents[nxext][nyext]);
    // fvl.resize(boost::extents[nxext][nyext]);
    // fvh.resize(boost::extents[nxext][nyext]);
    // gt.resize(boost::extents[nxext][nyext]);

    // sao.resize(boost::extents[nxext][nyext]);

    // u_pad.resize(boost::extents[nxext][nyext]);
    // v_pad.resize(boost::extents[nxext][nyext]);
    // scp2_pad.resize(boost::extents[nxext][nyext]);
    // scp2i_pad.resize(boost::extents[nxext][nyext]);
    // scuy_pad.resize(boost::extents[nxext][nyext]);
    // scvx_pad.resize(boost::extents[nxext][nyext]);
    // h_pad.resize(boost::extents[nxext][nyext]);

    //ncs = std::ceil(nwavedirn/2);
    ncs = std::round(nwavedirn/2);

    // cmom0.resize(boost::extents[nx][ny]);
    // cmom_dir.resize(boost::extents[nx][ny]);
    // CSfreq.resize(boost::extents[nx][ny]);
    // cmom_dir0.resize(boost::extents[nx][ny]);
    // CF.resize(boost::extents[nx][ny]);

    vbf = 0.1;
    vb = vbf;

    sigma_c  = (1.76e+6)*std::exp(-5.88*std::sqrt(vbf));
    epsc = sigma_c/young;
    flex_rig_coeff = young/(12.0*(1-std::pow(poisson,2.)));

}


template<typename T>
void
WimDiscr<T>::wimStep()
{

	// std::fill( tmp1.data(), tmp1.data() + tmp1.num_elements(), 0. );
    // std::fill( mom0.data(), mom0.data() + mom0.num_elements(), 0. );
    // std::fill( mom2.data(), mom2.data() + mom2.num_elements(), 0. );
    // std::fill( var_strain.data(), var_strain.data() + var_strain.num_elements(), 0. );
    // std::fill( mom0w.data(), mom0w.data() + mom0w.num_elements(), 0. );
    // std::fill( mom2w.data(), mom2w.data() + mom2w.num_elements(), 0. );


    std::fill( tau_x.data(), tau_x.data() + tau_x.num_elements(), 0. );
    std::fill( tau_y.data(), tau_y.data() + tau_y.num_elements(), 0. );


    array2_type tmp1, mom0, mom2, var_strain, mom0w, mom2w;
    tmp1.resize(boost::extents[nx][ny]);
    mom0.resize(boost::extents[nx][ny]);
    mom2.resize(boost::extents[nx][ny]);
    var_strain.resize(boost::extents[nx][ny]);
    mom0w.resize(boost::extents[nx][ny]);
    mom2w.resize(boost::extents[nx][ny]);


    value_type E_tot, sig_strain, Pstrain, P_crit, wlng_crest, Dc;
    value_type adv_dir, F, kicel, om, ommin, ommax, om1, lam1, lam2, dom, dave, c1d, tmp;
    int jcrest;
    bool break_criterion;

    dom = 2*PI*(freq_vec[nwavefreq-1]-freq_vec[0])/(nwavefreq-1.0);

    if (vm["steady"].template as<bool>())
    {
        for (int i = 0; i < nwavefreq; i++)
        {
            for (int j = 0; j < nwavedirn; j++)
            {
                adv_dir = (-PI/180)*(wavedir[j]+90.);

                if (std::cos(adv_dir) >= 0.)
                {
                    for (int k = 0; k < nx; k++)
                    {
                        for (int l = 0; l < ny; l++)
                        {
                            if (wave_mask2[k][l] >= 0.)
                            {
                                sdf_dir[k][l][j][i] = sdf_inc[k][l][j][i];
                                //std::cout<<"sdf_dir[k][l][j][i]= "<< sdf_dir[k][l][j][i] <<"\n";
                            }
                        }
                    }
                }

            }
        }
    }

    dave = 0;

    for (int fq = 0; fq < nwavefreq; fq++)
    {
        std::fill( atten_dim.data(), atten_dim.data() + atten_dim.num_elements(), 0. );
        std::fill( damp_dim.data(), damp_dim.data() + damp_dim.num_elements(), 0. );

        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                if ((ice_mask[i][j] == 1.) && (atten))
                {
                    if (dfloe[i][j] <200.)
                        floeScaling(dfloe[i][j],dave);
                    else
                        dave = dfloe[i][j];

                    // floes per unit length
                    c1d = icec[i][j]/dave;

                    // scattering
                    atten_dim[i][j] = atten_nond[i][j][fq]*c1d;

                    // damping
                    damp_dim[i][j] = 2*damping[i][j][fq]*icec[i][j];
                }
            }
        }


        // copy for application of advAttenSimple
        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                ag2d_eff_temp[i][j] = ag_eff[i][j][fq];

                for (int dn = 0; dn < nwavedirn; dn++)
                {
                    sdf3d_dir_temp[i][j][dn] = sdf_dir[i][j][dn][fq];
                }
            }
        }

        if (scatmod == "dissipated")
        {
            advAttenSimple(sdf3d_dir_temp, S_freq, taux_om, tauy_om, ag2d_eff_temp);
        }
        else if (scatmod == "isotropic")
        {
            advAttenIsotropic(sdf3d_dir_temp, S_freq, taux_om, tauy_om, ag2d_eff_temp);
        }

        // update after application of advAttenSimple
        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                ag_eff[i][j][fq] = ag2d_eff_temp[i][j];

                for (int dn = 0; dn < nwavedirn; dn++)
                {
                    sdf_dir[i][j][dn][fq] = sdf3d_dir_temp[i][j][dn];
                    // std::cout<<"AFTER: SDIR["<< i << "," << j << "]= "<< sdf_dir[i][j][dn][fq] <<"\n";
                }

                //std::cout<<"AFTER: taux_om["<< i << "," << j << "]= "<< taux_om[i][j] <<"\n";
            }
        }

        // for (int i = 0; i < nx; i++)
        // {
        //     for (int j = 0; j < ny; j++)
        //     {
        //         std::cout<<"iceh["<< i << "," << j << "]= "<< S_freq[i][j] <<"\n";
        //     }
        // }

        // integrate stress densities over frequency
        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                //std::cout<<"tau_x["<< i << "," << j << "]= "<< taux_om[i][j] <<"\n";
                tmp1[i][j] = rhowtr*gravity*taux_om[i][j]/ap_eff[i][j][fq];
                tau_x[i][j] += wt_om[fq]*tmp1[i][j];

                tmp1[i][j] = rhowtr*gravity*tauy_om[i][j]/ap_eff[i][j][fq];
                tau_y[i][j] += wt_om[fq]*tmp1[i][j];

                //std::cout<<"tau_x["<< i << "," << j << "]= "<< taux_om[i][j] <<"\n";
            }
        }

        // value_type _min = *std::min_element(S_freq.data(),S_freq.data() + S_freq.num_elements());
        // value_type _max = *std::max_element(S_freq.data(),S_freq.data() + S_freq.num_elements());
        // std::cout<<"Min= " << _min <<"\n";
        // std::cout<<"Max= " << _max <<"\n";

        // integrals for breaking program
        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                // convert from water amp's to ice amp's
                F = disp_ratio[i][j][fq];
                kicel = 2*PI/wlng_ice[i][j][fq];

                // spectral moments: take abs as small errors can make S_freq negative
                om = 2*PI*freq_vec[fq];
                tmp = wt_om[fq]*S_freq[i][j];
                // variance of displacement (water)
                mom0w[i][j] += std::abs(tmp);

                // variance of displacement (ice)
                mom0[i][j] += std::abs(tmp*std::pow(F,2.));
                //std::cout<<"mom0["<< i << "," << j << "]= "<< mom0w[i][j] <<"\n";

                tmp = wt_om[fq]*std::pow(om,2.)*S_freq[i][j];

                // variance of speed (water)
                mom2w[i][j] += std::abs(tmp);

                // variance of speed (ice)
                mom2[i][j] += std::abs(tmp*std::pow(F,2.));

                // variance of strain
                if (ice_mask[i][j] == 1.)
                {
                    // strain conversion factor
                    tmp = F*std::pow(kice,2.)*iceh[i][j]/2.0;
                    // strain density
                    tmp = wt_om[fq]*S_freq[i][j]*std::pow(tmp,2.);

                    var_strain[i][j] += std::abs(tmp);
                }
            }
        }
    }

    // for (int i = 0; i < nx; i++)
    //     for (int j = 0; j < ny; j++)
    //         std::cout << "VRT[" << i << "," << j << "]= " << var_strain[i][j] <<"\n";

    std::fill( Tp.data(), Tp.data() + Tp.num_elements(), 0. );

    if (ref_Hs_ice)
    {
        Hs = mom0;
        std::for_each(Hs.data(), Hs.data()+Hs.num_elements(), [&](value_type& f){ f = 4*std::sqrt(f); });
        // std::fill( Tp.data(), Tp.data() + Tp.num_elements(), 0. );

        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                if (mom2[i][j] > 0.)
                    Tp[i][j] = 2*PI*std::sqrt(mom0[i][j]/mom2[i][j]);
            }
        }
    }
    else
    {
        Hs = mom0w;
        std::for_each(Hs.data(), Hs.data()+Hs.num_elements(), [&](value_type& f){ f = 4*std::sqrt(f); });
        // std::fill( Tp.data(), Tp.data() + Tp.num_elements(), 0. );

        //std::fill( Hs.data(), Hs.data() + Hs.num_elements(), 0. );

        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                //Hs[i][j] = 4*std::sqrt(mom0w[i][j]);

                if (mom2w[i][j] > 0.)
                    Tp[i][j] = 2*PI*std::sqrt(mom0w[i][j]/mom2w[i][j]);
            }
        }
    }

    calcMWD();

    if (!steady && !breaking)
    {
        auto temparray = Hs;
        std::for_each(temparray.data(), temparray.data()+temparray.num_elements(), [&](value_type& f){ f *= f; });
        E_tot = std::accumulate(temparray.data(), temparray.data()+temparray.num_elements(),0.);

        // std::fill( var_strain.data(), var_strain.data()+var_strain.num_elements(), 1. );
        // E_tot = std::accumulate(var_strain.data(), var_strain.data()+var_strain.num_elements(),0.);
        // std::cout<<"Sum= "<< E_tot <<"\n";
    }

    // finally do floe breaking

    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            //std::cout << "MASK[" << i << "," << j << "]= " << ice_mask[i][j] << " and "<< mom0[i][j]  <<"\n";

            if ((ice_mask[i][j] == 1.) && (mom0[i][j] >= 0.))
            {
                // significant strain amp
                sig_strain = 2*std::sqrt(var_strain[i][j]);

                // probability of critical strain
                // being exceeded from Rayleigh distribution
                Pstrain = std::exp( -std::pow(epsc,2.)/(2*var_strain[i][j]) );
                P_crit = (1-breaking)+std::exp(-1.0);

                break_criterion = (Pstrain >= P_crit) && breaking;

                if (break_criterion)
                {
                    //std::cout << "TP[" << i << "," << j << "]= " << Tp[i][j] <<"\n";
                    om    = 2*PI/Tp[i][j];
                    ommin = 2*PI*freq_vec[0];
                    ommax = 2*PI*freq_vec[nwavefreq-1];

                    if (om <= ommin)
                        wlng_crest = wlng_ice[i][j][0];
                    else if (om >= ommax)
                        wlng_crest = wlng_ice[i][j][nwavefreq-1];
                    else
                    {
                        jcrest = std::floor((om-ommin+dom)/dom);
                        om1 = 2*PI*freq_vec[jcrest];
                        lam1 = wlng_ice[i][j][jcrest];
                        lam2 = wlng_ice[i][j][jcrest+1];
                        wlng_crest = lam1+(om-om1)*(lam2-lam1)/dom;
                    }

                    Dc = std::max<value_type>(dmin,wlng_crest/2.0);
                    dfloe[i][j] = std::min<value_type>(Dc,dfloe[i][j]);
                    //std::cout<<"DMAX= std::MAX("<< Dc << "," << dfloe[i][j] <<")\n";
                }
            }
            else if (wtr_mask[i][j] == 1.)
                dfloe[i][j] = 0;
        }
    }


    // for (int i = 0; i < nx; i++)
    //     for (int j = 0; j < ny; j++)
    //         std::cout << "Dmax[" << i << "," << j << "]= " << dfloe[i][j] <<"\n";

    // for (int i = 0; i < nx; i++)
    //     for (int j = 0; j < ny; j++)
    //         std::cout<<"Hs["<< i << "," << j << "]= "<< Hs[i][j] <<"\n";

    // std::cout<<"Hs_max= "<< *std::max_element(Hs.data(), Hs.data()+Hs.num_elements()) <<"\n";
    // std::cout<<"Hs_min= "<< *std::min_element(Hs.data(), Hs.data()+Hs.num_elements()) <<"\n";



}


template<typename T>
void
WimDiscr<T>::wimRun()
{
    value_type x_ext, y_ext, u_ref, duration;
    int nt;
    bool critter;

    std::fill( sdf_dir.data(), sdf_dir.data() + sdf_dir.num_elements(), 0. );

    for (int i = 0; i < nwavefreq; i++)
    {
        for (int j = 0; j < nwavedirn; j++)
        {
            for (int k = 0; k < nx; k++)
            {
                for (int l = 0; l < ny; l++)
                {
                    if (wave_mask[k][l] == 1.)
                    {
                        sdf_dir[k][l][j][i] = sdf_inc[k][l][j][i];
                        //std::cout<<"sdf_dir[k][l][j][i]= "<< sdf_dir[k][l][j][i] <<"\n";
                    }
                }
            }
        }
    }

    // value_type _min = *std::min_element(sdf_dir.data(),sdf_dir.data() + sdf_dir.num_elements());
    // value_type _max = *std::max_element(sdf_dir.data(),sdf_dir.data() + sdf_dir.num_elements());
    // std::cout<<"Min= " << _min <<"\n";
    // std::cout<<"Max= " << _max <<"\n";


    x_ext = nx*dx/1.e+03;
    y_ext = ny*dy/1.e+03;
    u_ref = amin + 0.7*(amax-amin);
    //duration = 1.0e3*x_ext/u_ref;
    duration = vm["duration"].template as<double>();

    std::cout<<"x_ext= "<< x_ext <<"\n";
    std::cout<<"y_ext= "<< y_ext <<"\n";
    std::cout<<"u_ref= "<< u_ref <<"\n";
    std::cout<<"duration= "<< duration <<"\n";
    std::cout<<"dt= "<< dt <<"\n";

    //nt = std::floor(duration/dt);
    nt = std::round(duration/dt);

    std::cout<<"nt= "<< nt <<"\n";

    int cpt = 0;
    //nt = 1;

    //readData("wim_prog001.a");
    readData("wim_prog10.a");

    while (cpt < nt)
    {
        std::cout <<  ":[WIM2D TIME STEP]^"<< cpt+1 <<"\n";

        value_type t_out = dt*cpt;
        //std::cout<<"T_OUT= "<< t_out <<"\n";

        critter = !(cpt % vm["reps"].template as<int>()) && (vm["checkprog"].template as<bool>());
        if (critter)
            writeFile(cpt,t_out);

        //wimStep();

        //critter = !(cpt % vm["reps"].template as<int>()) && (vm["checkprog"].template as<bool>());
        // if (critter)
        //     writeFile(cpt);

        array2_type diff;
        diff.resize(boost::extents[nx][ny]);

#if 1
        if (cpt==10)
        {
            for (int i = 0; i < nx; i++)
            {
                for (int j = 0; j < ny; j++)
                {
                    diff[i][j] = Hs[i][j]-Fhs[i][j];
                    //diff[i][j] = dfloe[i][j]-Fdmax[i][j];
                    //diff[i][j] = tau_x[i][j]-Ftaux[i][j];
                    //diff[i][j] = Tp[i][j]-Ftp[i][j];
                }
            }

            value_type error = *std::max_element(diff.data(),diff.data() + diff.num_elements());
            std::cout<<"Diff= " << error <<"\n";
        }

#endif

#if 1
        for (int i = 0; i < nx; i++)
        {
            if (cpt==10 && (i<10 || i==1))
                for (int j = 0; j < ny; j++)
                {
                    //std::cout<<"Hs[" << i <<"][" << j <<"]= " << Hs[i][j]-Fhs[i][j] <<"\n";

                    if (j==0)
                        std::cout<<"Hs[" << i <<"][" << j <<"]= " << Hs[i][j] << " : "
                                 << Fhs[i][j] << " . "<< Hs[i][j]-Fhs[i][j] <<"\n";
                }
        }
#endif
        wimStep();

        ++cpt;
    }
}

template<typename T>
void
WimDiscr<T>::floeScaling(value_type const& dmax, value_type& dave)
{
    int mm = 0;
    value_type n, nsum, nd, ndsum, r, dfac;

    r = dmax/dmin;
    n = nsum = nd = ndsum = dfac = 0;

    while (r > xi)
    {
        r  = r/xi;
        ++mm;
    }

    if (mm > 0)
    {
        for (int m=0; m<mm+1; ++m)
        {
            n = (1.0-fragility)*std::pow((fragility*std::pow(xi,2.)),m);
            nd = n/(std::pow(xi,m));
            nsum += n;
            ndsum += nd;
            dfac = ndsum/nsum;
        }

        dave = dfac*dmax;
    }
    else
    {
        dave = dmin;
    }
}

template<typename T>
void
WimDiscr<T>::advAttenSimple(array3_type& Sdir, array2_type& Sfreq, array2_type& taux_omega, array2_type& tauy_omega, array2_type const& ag2d_eff)
{

	// std::fill( uwave.data(), uwave.data() + uwave.num_elements(), 0. );
	// std::fill( vwave.data(), vwave.data() + vwave.num_elements(), 0. );
    // std::fill( temp.data(), temp.data() + temp.num_elements(), 0. );

    array2_type uwave, vwave, temp;
    uwave.resize(boost::extents[nx][ny]);
    vwave.resize(boost::extents[nx][ny]);
    temp.resize(boost::extents[nx][ny]);

	std::vector<value_type> wt_theta(nwavedirn);
	value_type adv_dir, S_th, tmp, alp_dim, source;

	for (int k = 0; k < nwavedirn; k++)
    {
        adv_dir = -PI*(90.0+wavedir[k])/180.0;
        //std::cout<<"adv_dir= "<< std::cos(adv_dir) <<"\n";
        uwave = ag2d_eff;
        std::for_each(uwave.data(), uwave.data()+uwave.num_elements(), [&](value_type& f){ f *= std::cos(adv_dir); });

        if (advdim == 2)
        {
            vwave = ag2d_eff;
            std::for_each(vwave.data(), vwave.data()+vwave.num_elements(), [&](value_type& f){ f *= std::sin(adv_dir); });
        }

        // copy from 3D input array to 2D temporary array
        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                temp[i][j] = Sdir[i][j][k];
            }
        }

#if 0
        if (k==3)
        {
            value_type _min = *std::min_element(temp.data(),temp.data() + temp.num_elements());
            value_type _max = *std::max_element(temp.data(),temp.data() + temp.num_elements());
            std::cout<<"Min 1= " << _min <<"\n";
            std::cout<<"Max 1= " << _max <<"\n";
        }
#endif
        // advection
        waveAdvWeno(temp,uwave,vwave);

#if 0
        if (k==3)
        {
            value_type _min = *std::min_element(temp.data(),temp.data() + temp.num_elements());
            value_type _max = *std::max_element(temp.data(),temp.data() + temp.num_elements());
            std::cout<<"Min 2= " << _min <<"\n";
            std::cout<<"Max 2= " << _max <<"\n";
        }
#endif

        // copy from 2D temporary array to 3D input array
        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                Sdir[i][j][k] = temp[i][j];
                //std::cout<<"BEFORE: Sdir["<< i << "," << j << "]= "<< Sdir[i][j][k] <<"\n";
            }
        }

    }

    if (nwavedirn == 1)
        std::fill( wt_theta.begin(), wt_theta.end(), 1. );
    else
        std::fill( wt_theta.begin(), wt_theta.end(), 2*PI/(1.0*nwavedirn) );


    std::fill( Sfreq.data(), Sfreq.data() + Sfreq.num_elements(), 0. );
    std::fill( taux_omega.data(), taux_omega.data() + taux_omega.num_elements(), 0. );
    std::fill( tauy_omega.data(), tauy_omega.data() + tauy_omega.num_elements(), 0. );

    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            if (ice_mask[i][j] > 0.)
            {
                for (int wnd = 0; wnd < nwavedirn; wnd++)
                {
                    adv_dir = -PI*(90.0+wavedir[wnd])/180.0;
                    S_th = Sdir[i][j][wnd];
                    alp_dim = atten_dim[i][j]+damp_dim[i][j];

                    // stress calculation
                    source = -alp_dim*ag2d_eff[i][j]*S_th;
                    tmp = -std::cos(adv_dir)*wt_theta[wnd]*source;
                    taux_omega[i][j] += tmp;
                    tmp = -std::sin(adv_dir)*wt_theta[wnd]*source;
                    tauy_omega[i][j] += tmp;

                    // do attenuation
                    Sdir[i][j][wnd] = S_th*std::exp(-alp_dim*ag2d_eff[i][j]*dt);

                    //std::cout<<"tau_x["<< i << "," << j << "]= "<< atten_dim[i][j] <<"\n";
                }
            }

            // integrate spectrum over direction
            for (int wnd = 0; wnd < nwavedirn; wnd++)
            {
                Sfreq[i][j] += wt_theta[wnd]*Sdir[i][j][wnd];
            }

            //std::cout<<"BEFORE: taux_om["<< i << "," << j << "]= "<< taux_om[i][j] <<"\n";
        }
    }

    // for (int k = 0; k < nwavedirn; k++)
    //     for (int i = 0; i < nx; i++)
    //     {
    //         for (int j = 0; j < ny; j++)
    //         {
    //             std::cout<<"MASK= "<< Sdir[i][j][k] <<"\n";
    //         }
    //     }



}


template<typename T>
void
WimDiscr<T>::advAttenIsotropic(array3_type& Sdir, array2_type& Sfreq, array2_type& taux_omega, array2_type& tauy_omega, array2_type const& ag2d_eff)
{
	// std::fill( uwave.data(), uwave.data() + uwave.num_elements(), 0. );
	// std::fill( vwave.data(), vwave.data() + vwave.num_elements(), 0. );

    array2_type uwave, vwave, temp;
    uwave.resize(boost::extents[nx][ny]);
    vwave.resize(boost::extents[nx][ny]);
    temp.resize(boost::extents[nx][ny]);

	std::vector<int> nvec(nwavedirn);
	std::vector<value_type> K_fou(nwavedirn), S_th(nwavedirn), theta_vec(nwavedirn), wt_theta(nwavedirn);
	std::vector<value_type> tmp1(nwavedirn), evals_x(nwavedirn);
	value_type adv_dir, tmp, alp_dim, source;

	std::vector<std::complex<value_type> > S_fou(nwavedirn);
	std::complex<value_type> zi, src_fou_p1, src_fou_m1;

	std::vector<value_type> S_cos(ncs), S_sin(ncs);
	value_type cg, q_scat, q_abs, q_tot, src_cos_1, src_sin_1;
    int n, jp1, jm1;

	zi = std::complex<value_type>(0.,1.);

	std::fill( wt_theta.begin(), wt_theta.end(), 2*PI/(1.0*nwavedirn) );

	for (int k = 0; k < nwavedirn; k++)
    {
        adv_dir = -PI*(90.0+wavedir[k])/180.0;

        uwave = ag2d_eff;

        std::for_each(uwave.data(), uwave.data()+uwave.num_elements(), [&](value_type& f){ f *= std::cos(adv_dir); });

        if (advdim == 2)
        {
	        vwave = ag2d_eff;
	        std::for_each(vwave.data(), vwave.data()+vwave.num_elements(), [&](value_type& f){ f *= std::sin(adv_dir); });
        }

        // copy from 3D input array to 2D temporary array
        for (int i = 0; i < nx; i++)
        {
	        for (int j = 0; j < ny; j++)
	        {
		        temp[i][j] = Sdir[i][j][k];
	        }
        }

        // advection
        waveAdvWeno(temp,uwave,vwave);

        // copy from 2D temporary array to 3D input array
        for (int i = 0; i < nx; i++)
        {
	        for (int j = 0; j < ny; j++)
	        {
		        Sdir[i][j][k] = temp[i][j];
	        }
        }

        theta_vec[k] = adv_dir;
        nvec[k] = k;
    }

    std::fill( Sfreq.data(), Sfreq.data() + Sfreq.num_elements(), 0. );
    std::fill( taux_omega.data(), taux_omega.data() + taux_omega.num_elements(), 0. );
    std::fill( tauy_omega.data(), tauy_omega.data() + tauy_omega.num_elements(), 0. );

	for (int i = 0; i < nx; i++)
    {
	    for (int j = 0; j < ny; j++)
	    {
		    for (int k = 0; k < nwavedirn; k++)
			    S_th[k] = Sdir[i][j][k];

		    std::fill( S_fou.begin(), S_fou.end(), zi );

		    //S_fou[0] = std::complex<value_type>( sum(wt_theta*S_th) );
            S_fou[0] = std::complex<value_type>( std::inner_product(wt_theta.begin(), wt_theta.end(), S_th.begin(), 0.) );

            if (ice_mask[i][j] > 0.)
            {
                if (dfloe[i][j] < dfloe_pack_init)
                {
                    q_scat = atten_dim[i][j];
                    q_abs = damp_dim[i][j];
                }
                else
                {
                    q_scat = 0;
                    q_abs = atten_dim[i][j]+damp_dim[i][j];
                }

                q_tot = q_scat+q_abs;
                cg = ag2d_eff[i][j];

                std::fill(K_fou.begin(), K_fou.end(), 0.);
                K_fou[0] = q_scat;

                std::fill(evals_x.begin(), evals_x.end(), -q_tot);
                evals_x[0] = -q_abs;

                // for (int cps = 0; cps < nwavedirn; cps++)
                //     std::cout<<"K[" << cps <<"]= " << S_th[cps] <<"\n";

                //std::cout<< "NCS= "<<ncs <<"\n";

                //ncs = 1;

                // std::vector<int> pos;
                // std::vector<int> neg;

                for (int nth = 0; nth < ncs; nth++)
                {
                    //std::vector<value_type> prodtmp(S_th.size());
                    std::vector<value_type> prodtmp = theta_vec;
                    std::for_each(prodtmp.begin(), prodtmp.end(), [&](value_type& f){ f = std::cos(nth*f); });
                    std::transform(prodtmp.begin(), prodtmp.end(), S_th.begin(), prodtmp.begin(),
                                   std::multiplies<value_type>());
                    S_cos[nth] = std::inner_product(prodtmp.begin(), prodtmp.end(), wt_theta.begin(), 0.);

                    // for (int cps = 0; cps < nwavedirn; cps++)
                    //     std::cout<<"K[" << cps <<"]= " << prodtmp[cps] << " and "<< wt_theta[cps] <<"\n";
                    // std::cout<<"taux= "<< S_cos[nth] <<"\n";

                    prodtmp = theta_vec;
                    std::for_each(prodtmp.begin(), prodtmp.end(), [&](value_type& f){ f = std::sin(nth*f); });
                    std::transform(prodtmp.begin(), prodtmp.end(), S_th.begin(), prodtmp.begin(),
                                   std::multiplies<value_type>());
                    S_sin[nth] = std::inner_product(prodtmp.begin(), prodtmp.end(), wt_theta.begin(), 0.);

                    // pos.push_back(nth+1);
                    S_fou[nth+1] = std::complex<value_type>(S_cos[nth],S_sin[nth]);
                    //S_fou[nwavedirn-nth] = std::complex<value_type>(S_cos[nth],-S_sin[nth]);

                    if (nth != ncs-1)
                    {
                        S_fou[nwavedirn-(nth+1)] = std::complex<value_type>(S_cos[nth],-S_sin[nth]);
                        // neg.push_back(nwavedirn-(nth+1));
                    }

                    // std::cout<<"nth= " << nth << ": and ncs+nth= "<< ncs+nth <<"\n";

                    // std::cout<<"taux= "<< S_sin[nth] <<"\n";

                }

                // stresses
                jp1 = 1;
                jm1 = nwavedirn-1;

                src_fou_p1 = cg*(-q_tot*S_fou[jp1]+K_fou[jp1]*S_fou[jp1]);
                src_fou_m1 = cg*(-q_tot*S_fou[jm1]+K_fou[jm1]*S_fou[jm1]);
                src_cos_1 = std::real(std::complex<value_type>(0.5)*(src_fou_p1+src_fou_m1));
                src_sin_1 = std::real(-zi*std::complex<value_type>(0.5)*(src_fou_p1-src_fou_m1));

                taux_omega[i][j] = -src_cos_1;
                tauy_omega[i][j] = -src_sin_1;

                std::vector<value_type> prodtmp = evals_x;
                std::for_each(prodtmp.begin(), prodtmp.end(), [&](value_type& f){ f = std::exp(cg*dt*f); });
                std::transform(S_fou.begin(), S_fou.end(), prodtmp.begin(), S_fou.begin(),
                               std::multiplies<std::complex<value_type> >());

                std::vector<std::complex<value_type> > Sfoutempcos = S_fou;
                std::vector<std::complex<value_type> > Sfoutempsin = S_fou;

#if 0
                // for (int cps = 0; cps < S_fou.size(); cps++)
                //     std::cout<<"S_fou["<< cps <<"]= "<< S_fou[cps] <<"\n";

                // for (int cps = 0; cps < S_fou.size(); cps++)
                //     std::cout<<"S_fou["<< cps <<"]= "<< Sfoutempcos[cps] <<"\n";

                prodtmp = theta_vec;
                std::transform(prodtmp.begin(), prodtmp.end(), nvec.begin(), prodtmp.begin(),
                               std::multiplies<value_type>());
                std::for_each(prodtmp.begin(), prodtmp.end(), [&](value_type& f){ f = std::cos(f); });

                std::fill(prodtmp.begin(), prodtmp.end(), 1.);
                std::transform(Sfoutempcos.begin(), Sfoutempcos.end(), prodtmp.begin(), Sfoutempcos.begin(),
                               std::multiplies<std::complex<value_type> >());

                for (int cps = 0; cps < S_fou.size(); cps++)
                {
                    //std::cout<<"S_fou["<< cps <<"]= "<< Sfoutempcos[cps] <<"\n";
                    std::cout<<"S_fou["<< cps <<"]= Real= "
                             << std::real(Sfoutempcos[cps])-std::imag(Sfoutempcos[cps])
                             <<"\n"; // << " and Imag= "<< std::imag(Sfoutempcos[cps]) <<"\n";
                    //tmp1[kin] = std::real(Sfoutempcos[kin])-std::imag(Sfoutempcos[kin]);
                }
#endif
                for (int k = 0; k < nwavedirn; k++)
                {
                    // Sfoutempcos
                    // prodtmp = std::vector<value_type>(nvec.begin(), nvec.end()); // also works

                    prodtmp.clear();
                    prodtmp.insert(prodtmp.end(), nvec.begin(), nvec.end());

                    // for (int cvt = 0; cvt < prodtmp.size(); cvt++)
                    //     std::cout<<"prodtmp["<< cvt <<"]= "<< prodtmp[cvt] <<" : and "<< "nvec["<< cvt <<"]= "<< nvec[cvt] <<"\n";

                    // std::transform(prodtmp.begin(), prodtmp.end(), nvec.begin(), prodtmp.begin(), std::multiplies<value_type>());
                    std::for_each(prodtmp.begin(), prodtmp.end(), [&](value_type& f){ f = std::cos(theta_vec[k]*f); });
                    std::transform(Sfoutempcos.begin(), Sfoutempcos.end(), prodtmp.begin(), Sfoutempcos.begin(),
                                   std::multiplies<std::complex<value_type> >());


                    // Sfoutempsin
                    // prodtmp = std::vector<value_type>(nvec.begin(), nvec.end()); // also works

                    prodtmp.clear();
                    prodtmp.insert(prodtmp.end(), nvec.begin(), nvec.end());

                    // for (int cvt = 0; cvt < prodtmp.size(); cvt++)
                    //     std::cout<<"prodtmp["<< cvt <<"]= "<< prodtmp[cvt] <<" : and "<< "nvec["<< cvt <<"]= "<< nvec[cvt] <<"\n";

                    // std::transform(prodtmp.begin(), prodtmp.end(), nvec.begin(), prodtmp.begin(), std::multiplies<value_type>());
                    std::for_each(prodtmp.begin(), prodtmp.end(), [&](value_type& f){ f = std::sin(theta_vec[k]*f); });
                    std::transform(Sfoutempsin.begin(), Sfoutempsin.end(), prodtmp.begin(), Sfoutempsin.begin(),
                                   std::multiplies<std::complex<value_type> >());

                    for (int kin = 0; kin < Sfoutempcos.size(); kin++)
                    {
                        tmp1[kin] = std::real(Sfoutempcos[kin])-std::imag(Sfoutempsin[kin]);
                    }

                    //Sdir[i][j][k] = std::inner_product(tmp1.begin(), tmp1.end(), tmp1.begin(), 0.)/(2*PI);
                    Sdir[i][j][k] = std::accumulate(tmp1.begin(), tmp1.end(), 0.0)/(2*PI);

                    // std::cout<<"Sdir[i][j][k]= "<< Sdir[i][j][k] <<"\n";
                    // std::cout<<"std::numeric_limits<short>::max()=  "<< std::numeric_limits<float>::min() <<"\n";

                    S_th[k] = Sdir[i][j][k];
                }
            }

            Sfreq[i][j] = std::real(S_fou[0]);
	    }
    }

    // std::cout<<"max_element= " << *std::min_element(Sdir.data(), Sdir.data() + Sdir.num_elements()) <<"\n";

}

template<typename T>
void
WimDiscr<T>::waveAdvWeno(array2_type& h, array2_type const& u, array2_type const& v)
{

	// std::fill( sao.data(), sao.data() + sao.num_elements(), 0. );
    array2_type sao;
    array2_type u_pad, v_pad, scp2_pad, scp2i_pad, scuy_pad, scvx_pad, h_pad;
    sao.resize(boost::extents[nxext][nyext]);

    array2_type hp_temp;
    hp_temp.resize(boost::extents[nx][ny]);

	padVar(u, u_pad);
	padVar(v, v_pad);
	padVar(SCP2_array, scp2_pad);
    padVar(SCP2I_array, scp2i_pad);
    padVar(SCUY_array, scuy_pad);
    padVar(SCVX_array, scvx_pad);
    padVar(h, h_pad);

    //std::fill( u_pad.data(), u_pad.data() + u_pad.num_elements(), 1. );
    //value_type sum_ = std::inner_product(u.data(), u.data() + u.num_elements(), u.data(), 0.);
    //value_type sum_ = std::accumulate(u_pad.data(), u_pad.data() + u_pad.num_elements(), 0.0);
    //std::cout<<"*********sum= "<< sum_ <<"\n";


    // prediction step
    weno3pdV2(h_pad, u_pad, v_pad, scuy_pad, scvx_pad, scp2i_pad, scp2_pad, sao);

#if 0
    for (int i = 0; i < nxext; i++)
    {
        for (int j = 0; j < nyext; j++)
        {
            hp[i][j] = h_pad[i][j]+dt*sao[i][j];
        }
    }
#endif

    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            hp[i+nbdx][j+nbdy] = h_pad[i+nbdx][j+nbdy]+dt*sao[i+nbdx][j+nbdy];
        }
    }

    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            hp_temp[i][j] = hp[i+nbdx][j+nbdy];
        }
    }

    padVar(hp_temp,hp);

    // correction step
    weno3pdV2(hp, u_pad, v_pad, scuy_pad, scvx_pad, scp2i_pad, scp2_pad, sao);

#if 0
    for (int i = 0; i < nxext; i++)
    {
        for (int j = 0; j < nyext; j++)
        {
            h_pad[i][j] = 0.5*(h_pad[i][j]+hp[i][j]+dt*sao[i][j]);
        }
    }

    for (int i = nbdx; i < nx+nbdx; i++)
    {
        for (int j = nbdy; j < ny+nbdy; j++)
        {
            h[i-nbdx][j-nbdy] = h_pad[i][j];

            // mask land (no waves on land)
            h[i-nbdx][j-nbdy] = h[i-nbdx][j-nbdy]*(1-LANDMASK_array[i-nbdx][j-nbdy]);
        }
    }
#endif

#if 1
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            h[i][j] = 0.5*(h_pad[i+nbdx][j+nbdy]+hp[i+nbdx][j+nbdy]+dt*sao[i+nbdx][j+nbdy]);
        }
    }

    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            h[i][j] = h[i][j]*(1-LANDMASK_array[i][j]);
        }
    }
#endif

}

template<typename T>
void
WimDiscr<T>::weno3pdV2(array2_type const& gin, array2_type const& u, array2_type const& v, array2_type const& scuy,
                       array2_type const& scvx, array2_type const& scp2i, array2_type const& scp2, array2_type& saoout)
{

	value_type cq00=-1./2 ,cq01=3./2, cq10=1./2, cq11=1./2, ca0=1./3, ca1=2./3, eps=1e-12;
	value_type q0, q1, a0, a1, q;
	int im1, im2, ip1, jm1, jm2, jp1, ymargin;

    array2_type ful, fuh, fvl, fvh, gt;
    ful.resize(boost::extents[nxext][nyext]);
    fuh.resize(boost::extents[nxext][nyext]);
    fvl.resize(boost::extents[nxext][nyext]);
    fvh.resize(boost::extents[nxext][nyext]);
    gt.resize(boost::extents[nxext][nyext]);

    // std::fill( ful.data(), ful.data() + ful.num_elements(), 0. );
    // std::fill( fuh.data(), fuh.data() + fuh.num_elements(), 0. );
    // std::fill( fvl.data(), fvl.data() + fvl.num_elements(), 0. );
    // std::fill( fvh.data(), fvh.data() + fvh.num_elements(), 0. );
    // std::fill( gt.data(), gt.data() + gt.num_elements(), 0. );

    // value_type sum_ = std::accumulate(u.data(), u.data() + u.num_elements(), 0.0);
    // //value_type sum_ = std::inner_product(u.data(), u.data() + u.num_elements(), u.data(), 0.);
    // std::cout<<"*********sum= "<< sum_ <<"\n";

    if (advdim == 2)
        ymargin = 1;
    else
        ymargin = 0;

    // fluxes in x direction
    for (int i = nbdx-1; i < nx+nbdx+2; i++)
    {
        for (int j = nbdy-ymargin; j < ny+nbdy+ymargin; j++)
        {
            im1 = i-1;

            if (u[i][j] > 0.)
            {
                // coefficents to calc higher-order fluxes
                im2 = im1-1;
                q0 = cq00*gin[im2][j]+cq01*gin[im1][j];
                q1 = cq10*gin[im1][j]+cq11*gin[i][j];
                a0 = ca0;
                a1 = ca1*(std::abs(gin[im2][j]-gin[im1][j])+eps)/(std::abs(gin[im1][j]-gin[i][j])+eps);

                // lower-order fluxes
                ful[i][j] = u[i][j]*gin[im1][j]*scuy[i][j];

            }
            else
            {
                // coefficents to calc higher-order fluxes
                ip1 = i+1;
                q0 = cq11*gin[im1][j]+cq10*gin[i][j];
                q1 = cq01*gin[i][j]+cq00*gin[ip1][j];
                a0 = ca1;
                a1 = ca0*(abs(gin[im1][j]-gin[i][j])+eps)/(abs(gin[i][j]-gin[ip1][j])+eps);

                // lower-order fluxes
                ful[i][j] = u[i][j]*gin[i][j]*scuy[i][j];
            }

            // higher-order fluxes
            fuh[i][j] = (u[i][j]*(a0*q0+a1*q1)*scuy[i][j]/(a0+a1))-ful[i][j];
        }
    }

    // fluxes in y direction
    if (advdim == 2)
    {
        for (int i = nbdx-1; i < nx+nbdx+1; i++)
        {
            for (int j = nbdy-1; j < ny+nbdy+2; j++)
            {
                jm1 = j-1;

                if (v[i][j] > 0.)
                {
                    jm2 = jm1-1;
                    q0 = cq00*gin[i][jm2]+cq01*gin[i][jm1];
                    q1 = cq10*gin[i][jm1]+cq11*gin[i][j];
                    a0 = ca0;
                    a1 = ca1*(std::abs(gin[i][jm2]-gin[i][jm1])+eps)/(std::abs(gin[i][jm1]-gin[i][j])+eps);
                    fvl[i][j] = v[i][j]*gin[i][jm1]*scvx[i][j];
                }
                else
                {
                    jp1 = j+1;
                    q0 = cq11*gin[i][jm1]+cq10*gin[i][j];
                    q1 = cq01*gin[i][j]+cq00*gin[i][jp1];
                    a0 = ca1;
                    a1 = ca0*(abs(gin[i][jm1]-gin[i][j])+eps)/(abs(gin[i][j]-gin[i][jp1])+eps);
                    fvl[i][j] = v[i][j]*gin[i][j]*scvx[i][j];
                }

                fvh[i][j] = (v[i][j]*(a0*q0+a1*q1)*scvx[i][j]/(a0+a1))-fvl[i][j];
            }
        }
    }

    // update field with low order fluxes
    for (int i = nbdx-1; i < nx+nbdx+1; i++)
    {
        for (int j = nbdy-ymargin; j < ny+nbdy+ymargin; j++)
        {
            if (advdim == 2)
            {
                gt[i][j] = gin[i][j]-dt*(ful[i+1][j]-ful[i][j]+fvl[i][j+1]-fvl[i][j])*scp2i[i][j];
            }
            else if (advdim == 1)
            {
                gt[i][j] = gin[i][j]-dt*(ful[i+1][j]-ful[i][j])*scp2i[i][j];
            }
        }
    }

    q = 0.25/dt;

    // // obtain fluxes with limited high order correction fluxes
    for (int i = nbdx; i < nx+nbdx+1; i++)
    {
        for (int j = nbdy; j < ny+nbdy; j++)
        {
            fuh[i][j] = ful[i][j]+std::max(-q*gt[i][j]*scp2[i][j],std::min(q*gt[i-1][j]*scp2[i-1][j],fuh[i][j]));
        }
    }

    // obtain fluxes with limited high order correction fluxes
    if (advdim == 2)
    {
        for (int i = nbdx; i < nx+nbdx; i++)
        {
            for (int j = nbdy; j < ny+nbdy+1; j++)
            {
                fvh[i][j]=fvl[i][j]+std::max(-q*gt[i][j]*scp2[i][j],std::min(q*gt[i][j-1]*scp2[i][j-1],fvh[i][j]));
            }
        }
    }

#if 1
    // compute the spatial advective operator
    for (int i = nbdx; i < nx+nbdx; i++)
    {
        for (int j = nbdy; j < ny+nbdy; j++)
        {
            if (advdim == 2)
            {
                saoout[i][j] = -(fuh[i+1][j]-fuh[i][j]+fvh[i][j+1]-fvh[i][j])*scp2i[i][j];
            }
            else if (advdim == 1)
            {
                saoout[i][j] = -(fuh[i+1][j]-fuh[i][j])*scp2i[i][j];
            }
        }
    }
#endif

}

template<typename T>
void
WimDiscr<T>::padVar(array2_type const& u, array2_type& upad)
{
    upad.resize(boost::extents[nxext][nyext]);


    for (int i = 0; i < nxext; i++)
    {
        for (int j = 0; j < nyext; j++)
        {

            if ((nbdx-1 < i) && (i < nx+nbdx) && (nbdy-1 < j) && (j < ny+nbdy))
            {
                upad[i][j] = u[i-nbdx][j-nbdy];
            }

            if (advdim == 1)
            {
                // make periodic in i
                if ((i < nbdx) && (nbdy-1 < j) && (j < ny+nbdy))
                {
                    upad[i][j] = u[nx-nbdx+i][j-nbdy];
                }

                if ((nx+nbdx-1 < i) && (nbdy-1 < j) && (j < ny+nbdy))
                {
                    upad[i][j] = u[i-nx-nbdx][j-nbdy];
                }
            }
            else if (advdim == 2)
            {
                if (advopt != "notperiodic")
                {
                    // make periodic in j
                    // if ((j < nbdy) && (nbdx-1 < i) && (i < nx+nbdx))
                    //     upad[i][j] = u[i-nbdx][ny-nbdy+j];

                    if ((j < nbdy) && (nbdx-1 < i) && (i < nx+nbdx))
                        upad[i][j] = u[i-nbdx][ny-nbdy+j-1];

                    if ((ny+nbdy-1 < j) && (nbdx-1 < i) && (i < nx+nbdx))
                        upad[i][j] = u[i-nbdx][j-ny-nbdy];
                }

                if (advopt == "xy-periodic")
                {
                    // make periodic in i
                    if ((i < nbdx) && (nbdy-1 < j) && (j < ny+nbdy))
                        upad[i][j] = u[nx-nbdx+i][j-nbdy];


                    if ((nx+nbdx-1 < i) && (nbdy-1 < j) && (j < ny+nbdy))
                        upad[i][j] = u[i-nx-nbdx][j-nbdy];


                    // // make periodic in j
                    // if ((j < nbdy) && (nbdy-1 < i) && (i < nx+nbdy))
                    //     upad[i][j] = u[i-nbdy][ny-nbdy+j];


                    if ((ny+nbdy-1 < j) && (nbdx-1 < i) && (i < nx+nbdx))
                        upad[i][j] = u[i-nbdx][j-ny-nbdy];


                    // BR, TL
                    if ((nx+nbdx-1 < i) && (ny+nbdy-1 < j))
                        upad[i][j] = u[i-nx-nbdx][j-ny-nbdy];


                    if ((i < nbdx) && (j < nbdy))
                        upad[i][j] = u[i+nx-nbdx][j];


                    // BL, TR
                    if ((nx+nbdx-1 < i) && (j < nbdy))
                        upad[i][j] = u[i-nx-nbdx][ny-nbdy+j];


                    if ((i < nbdx) && (ny+nbdy-1 < j))
                        upad[i][j] = u[i+nx-nbdx][j-ny-nbdy];

                }
            }
        }
    }
}

template<typename T>
void
WimDiscr<T>::calcMWD()
{
    value_type adv_dir, wt_theta, om;

    array2_type cmom0,cmom_dir,CSfreq, cmom_dir0, CF;
    cmom0.resize(boost::extents[nx][ny]);
    cmom_dir.resize(boost::extents[nx][ny]);
    CSfreq.resize(boost::extents[nx][ny]);
    cmom_dir0.resize(boost::extents[nx][ny]);
    CF.resize(boost::extents[nx][ny]);

    if (nwavedirn == 1)
        wt_theta = 1.;
    else
        wt_theta = 2.0*PI/(nwavedirn);

    // spectral moments
    std::fill( cmom0.data(), cmom0.data() + cmom0.num_elements(), 0. );
    std::fill( cmom_dir.data(), cmom_dir.data() + cmom_dir.num_elements(), 0. );

    for (int fq = 0; fq < nwavefreq; fq++)
    {
        om = 2*PI*freq_vec[fq];

        std::fill( CSfreq.data(), CSfreq.data() + CSfreq.num_elements(), 0. );
        std::fill( cmom_dir0.data(), cmom_dir0.data() + cmom_dir0.num_elements(), 0. );

        for (int dn = 0; dn < nwavedirn; dn++)
        {
            adv_dir = -PI*(90.0+wavedir[dn])/180.0;

            for (int i = 0; i < nx; i++)
            {
                for (int j = 0; j < ny; j++)
                {
                    CSfreq[i][j] += wt_theta*sdf_dir[i][j][dn][fq];
                    cmom_dir0[i][j] += wt_theta*sdf_dir[i][j][dn][fq]*adv_dir;
                }
            }
        }


        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                if (ref_Hs_ice)
                    CF[i][j] = disp_ratio[i][j][fq];
                else
                    CF[i][j] = 1;

                cmom0[i][j] += std::abs(wt_om[fq]*std::pow(CF[i][j],2.)*CSfreq[i][j]);
                cmom_dir[i][j] += std::abs(wt_om[fq]*std::pow(CF[i][j],2.)*cmom_dir0[i][j]);
            }
        }
    }

    std::fill( mwd.data(), mwd.data() + mwd.num_elements(), 0. );

    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            if (cmom0[i][j] > 0.)
                mwd[i][j] = -90.-180*(cmom_dir[i][j]/cmom0[i][j])/PI;
        }
    }

}

template class WimDiscr<float>;

} // namespace WIM2D
