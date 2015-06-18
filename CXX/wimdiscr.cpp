/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

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

    nbdy = vm["nbdy"].template as<int>();
    nxext = nx+2*nbdy;
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
    mwd_inc = -90.;
    Tmin = 2.5;
    Tmax = 25.;
    gravity = 9.81;

    unifc = 0.75;
    unifh = 2.0;
    dfloe_pack_init = 250.0;

    rho = 0.9;
    rhowtr = 1025.;
    rhow = 1025.;
    rhoi = rho*rhow;
    rhoice = 922.5;
    poisson = 0.3;
    dmin = 20.;
    xi = 2.;
    fragility = 0.9;

    young = vm["young"].template as<double>();
    visc_rp = vm["viscrp"].template as<double>();


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
            //std::cout<<"wavedir["<< i <<"]= "<< theta_max+i*dtheta <<"\n";
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


    x0 = X_array[0][0];
    xmax = x0+(nx-1)*dx;

    y0 = Y_array[0][0];
    ym = Y_array[0][0]+(ny-1)*dy;

    x_edge = 0.5*(x0+xmax)-0.8*(0.5*(xmax-x0));

    for (size_type i = 0; i < nx; i++)
    {
        for (size_type j = 0; j < ny; j++)
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

    // wtr_mask
    for (size_type i = 0; i < nx; i++)
    {
        for (size_type j = 0; j < ny; j++)
        {
            if (X_array[i][j] < x_edge)
                wtr_mask[i][j] = 1.;

            ice_mask[i][j] = (1-wtr_mask[i][j])*(1-LANDMASK_array[i][j]);

            icec[i][j] = unifc*ice_mask[i][j];
            iceh[i][j] = unifh*ice_mask[i][j];
            dfloe[i][j] = dfloe_pack_init*ice_mask[i][j];

            //std::cout<<"ICE_MASK["<< i  << "," << j << "]= " << ice_mask[i][j] <<"\n";
            //std::cout<<"ICEC["<< i  << "," << j << "]= " << icec[i][j] <<"\n";
            // if (j==0)
            //     std::cout<<"ICEH["<< i  << "," << j << "]= " << iceh[i][j] <<"\n";
            //std::cout<<"DFLOE["<< i  << "," << j << "]= " << dfloe[i][j] <<"\n";
            //std::cout<<"LANDMASK_array["<< i  << "," << j << "]= " << LANDMASK_array[i][j] <<"\n";
        }
    }



    om = 2*PI*freq_vec[0];
    std::vector<value_type> Sfreq(nwavefreq);
    std::vector<value_type> theta_fac(nwavedirn,0.);
    value_type f1, f2, f3, t_m, om_m, chi;

    for (size_type i = 0; i < nx; i++)
    {
        for (size_type j = 0; j < ny; j++)
        {
            if (wave_mask[i][j] == 1.)
            {
                if (nwavefreq > 1)
                {
                    for (size_type fq = 0; fq < nwavefreq; fq++)
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
                for (size_type dn = 0; dn < nwavedirn; dn++)
                {
                    chi = PI*(wavedir[dn]-mwd[i][j])/180.0;
                    if (std::cos(chi) > 0.)
                        theta_fac[dn] = 2.0*std::pow(cos(chi),2.)/PI;
                    else
                        theta_fac[dn] = 0.;
                }
            }

            // combine freq and dir
            for (size_type fq = 0; fq < nwavefreq; fq++)
            {
                for (size_type dn = 0; dn < nwavedirn; dn++)
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
    params[2] = rhow;
    params[3] = rhoi;
    params[4] = poisson;

    double outputs[8];

    for (size_type fq = 0; fq < nwavefreq; fq++)
    {
        for (size_type i = 0; i < nx; i++)
        {
            for (size_type j = 0; j < ny; j++)
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

                    disp_ratio[i][j][fq] = (kice*modT)/kwtr;
                    wlng_ice[i][j][fq] = 2*PI/kice;


                    ag_eff[i][j][fq] = ag[fq];
                    ap_eff[i][j][fq] = ap[fq];

                    if (icevel)
                    {
                        wlng_ice[i][j][fq] = wlng[fq];
                    }

                    // std::cout<<"damping[i][j][fq]= "<< damping[i][j][fq] <<"\n";
                    // std::cout<<"kice= "<< kice <<"\n";
                    // std::cout<<"kwtr= "<< kwtr <<"\n";
                    // std::cout<<"int_adm= "<< int_adm <<"\n";
                    // std::cout<<"atten_nond[i][j][fq]= "<< atten_nond[i][j][fq] <<"\n";
                    // std::cout<<"modT= " << modT <<"\n";
                    // std::cout<< "argR= "<< argR <<"\n";
                    // std::cout<<"argT= "<< argT <<"\n";
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
    tmp1.resize(boost::extents[nx][ny]);
    mom0.resize(boost::extents[nx][ny]);
    mom2.resize(boost::extents[nx][ny]);
    var_strain.resize(boost::extents[nx][ny]);
    mom0w.resize(boost::extents[nx][ny]);
    mom2w.resize(boost::extents[nx][ny]);


    // warning: in wimStep()
    uwave.resize(boost::extents[nx][ny]);
    vwave.resize(boost::extents[nx][ny]);
    temp.resize(boost::extents[nx][ny]);

    S_freq.resize(boost::extents[nx][ny]);
    taux_om.resize(boost::extents[nx][ny]);
    tauy_om.resize(boost::extents[nx][ny]);


    hp.resize(boost::extents[nxext][nyext]);


    ful.resize(boost::extents[nxext][nyext]);
    fuh.resize(boost::extents[nxext][nyext]);
    fvl.resize(boost::extents[nxext][nyext]);
    fvh.resize(boost::extents[nxext][nyext]);
    gt.resize(boost::extents[nxext][nyext]);

    sao.resize(boost::extents[nxext][nyext]);

    u_pad.resize(boost::extents[nxext][nyext]);
    v_pad.resize(boost::extents[nxext][nyext]);
    scp2_pad.resize(boost::extents[nxext][nyext]);
    scp2i_pad.resize(boost::extents[nxext][nyext]);
    scuy_pad.resize(boost::extents[nxext][nyext]);
    scvx_pad.resize(boost::extents[nxext][nyext]);
    h_pad.resize(boost::extents[nxext][nyext]);

    ncs = nwavedirn/2;

    cmom0.resize(boost::extents[nx][ny]);
    cmom_dir.resize(boost::extents[nx][ny]);
    CSfreq.resize(boost::extents[nx][ny]);
    cmom_dir0.resize(boost::extents[nx][ny]);
    CF.resize(boost::extents[nx][ny]);

    vbf = 0.1;
    vb = vbf;

    sigma_c  = double(1.76e+6)*std::exp(-5.88*std::sqrt(vbf));
    epsc = sigma_c/young;
    flex_rig_coeff = young/(12.0*(1-std::pow(poisson,2.)));

}


template<typename T>
void
WimDiscr<T>::wimStep()
{

	std::fill( tmp1.data(), tmp1.data() + tmp1.num_elements(), 0. );
    std::fill( mom0.data(), mom0.data() + mom0.num_elements(), 0. );
    std::fill( mom2.data(), mom2.data() + mom2.num_elements(), 0. );
    std::fill( var_strain.data(), var_strain.data() + var_strain.num_elements(), 0. );
    std::fill( mom0w.data(), mom0w.data() + mom0w.num_elements(), 0. );
    std::fill( mom2w.data(), mom2w.data() + mom2w.num_elements(), 0. );


    std::fill( tau_x.data(), tau_x.data() + tau_x.num_elements(), 0. );
    std::fill( tau_y.data(), tau_y.data() + tau_y.num_elements(), 0. );


    value_type E_tot, sig_strain, Pstrain, P_crit, wlng_crest, Dc;
    value_type adv_dir, F, kicel, om, ommin, ommax, om1, lam1, lam2, dom, dave, c1d, tmp;
    size_type jcrest;
    bool break_criterion;

    dom = 2*PI*(freq_vec[nwavefreq-1]-freq_vec[0])/(nwavefreq-1.0);

    if (vm["steady"].template as<bool>())
    {
        for (size_type i = 0; i < nwavefreq; i++)
        {
            for (size_type j = 0; j < nwavedirn; j++)
            {
                adv_dir = (-PI/180)*(wavedir[j]+90.);
                //std::cout<<"COS(adv_dir)= "<< std::cos(adv_dir) <<"\n";

                if (std::cos(adv_dir) >= 0.)
                {
                    for (size_type k = 0; k < nx; k++)
                    {
                        for (size_type l = 0; l < ny; l++)
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

    for (size_type fq = 0; fq < nwavefreq; fq++)
    {
        std::fill( atten_dim.data(), atten_dim.data() + atten_dim.num_elements(), 0. );
        std::fill( damp_dim.data(), damp_dim.data() + damp_dim.num_elements(), 0. );

        for (size_type i = 0; i < nx; i++)
        {
            for (size_type j = 0; j < ny; j++)
            {
                if (ice_mask[i][j] == 1. && atten)
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

        if (scatmod == "dissipated")
        {

            // copy for application of advAttenSimple
            for (size_type i = 0; i < nx; i++)
            {
                for (size_type j = 0; j < ny; j++)
                {
                    ag2d_eff_temp[i][j] = ag_eff[i][j][fq];

                    for (size_type dn = 0; dn < nwavedirn; dn++)
                    {
                        sdf3d_dir_temp[i][j][dn] = sdf_dir[i][j][dn][fq];
                    }
                }
            }

            advAttenSimple(sdf3d_dir_temp, S_freq, taux_om, tauy_om, ag2d_eff_temp);

            // update after application of advAttenSimple
            for (size_type i = 0; i < nx; i++)
            {
                for (size_type j = 0; j < ny; j++)
                {
                    ag_eff[i][j][fq] = ag2d_eff_temp[i][j];
                    //std::cout<<"ag_eff= "<< ag_eff[i][j][fq] <<"\n";

                    for (size_type dn = 0; dn < nwavedirn; dn++)
                    {
                        sdf_dir[i][j][dn][fq] = sdf3d_dir_temp[i][j][dn];
                    }
                }
            }



            // for (size_type k = 0; k < nwavedirn; k++)
            //     for (size_type i = 0; i < nx; i++)
            //     {
            //         for (size_type j = 0; j < ny; j++)
            //         {
            //             std::cout<<"MASK= "<< sdf3d_dir_temp[i][j][k] <<"\n";
            //         }
            //     }


            // for (size_type i = 0; i < nx; i++)
            // {
            //     for (size_type j = 0; j < ny; j++)
            //     {
            //         std::cout<<"iceh["<< i << "," << j << "]= "<< S_freq[i][j] <<"\n";
            //     }
            // }



            // integrate stress densities over frequency
            for (size_type i = 0; i < nx; i++)
            {
                for (size_type j = 0; j < ny; j++)
                {
                    //std::cout<<"tau_x["<< i << "," << j << "]= "<< taux_om[i][j] <<"\n";
                    tmp1[i][j] = rhowtr*gravity*taux_om[i][j]/ap_eff[i][j][fq];
                    tau_x[i][j] += wt_om[fq]*tmp1[i][j];

                    tmp1[i][j] = rhowtr*gravity*tauy_om[i][j]/ap_eff[i][j][fq];
                    tau_y[i][j] += wt_om[fq]*tmp1[i][j];
                }
            }

            // size_type ctr=0;
            // integrals for breaking program
            for (size_type i = 0; i < nx; i++)
            {
                for (size_type j = 0; j < ny; j++)
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

                    tmp = wt_om[fq]*std::pow(om,2.)*S_freq[i][j];

                    // variance of speed (water)
                    mom2w[i][j] += std::abs(tmp);

                    // variance of speed (ice)
                    mom2[i][j] += std::abs(tmp*std::pow(F,2.));

                    // variance of strain
                    if (ice_mask[i][j] == 1.)
                    {
                        //++ ctr;
                        // strain conversion factor
                        tmp = F*std::pow(kice,2.)*iceh[i][j]/2.0;
                        // strain density
                        tmp = wt_om[fq]*S_freq[i][j]*std::pow(tmp,2.);

                        var_strain[i][j] += std::abs(tmp);
                    }
                }
            }
        }
    }


    // for (size_type i = 0; i < nx; i++)
    //     for (size_type j = 0; j < ny; j++)
    //         std::cout << "VRT[" << i << "," << j << "]= " << var_strain[i][j] <<"\n";

    if (ref_Hs_ice)
    {
        Hs = mom0;
        std::for_each(Hs.data(), Hs.data()+Hs.num_elements(), [&](value_type& f){ f = 4*std::sqrt(f); });

        std::fill( Tp.data(), Tp.data() + Tp.num_elements(), 0. );

        for (size_type i = 0; i < nx; i++)
        {
            for (size_type j = 0; j < ny; j++)
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

        std::fill( Tp.data(), Tp.data() + Tp.num_elements(), 0. );

        for (size_type i = 0; i < nx; i++)
        {
            for (size_type j = 0; j < ny; j++)
            {
                if (mom2w[i][j] > 0.)
                    Tp[i][j] = 2*PI*std::sqrt(mom0w[i][j]/mom2w[i][j]);
            }
        }
    }


    calcMWD();

    // for (size_type i = 0; i < nx; i++)
    //     for (size_type j = 0; j < ny; j++)
    //         std::cout << "Hs[" << i << "," << j << "]= " << Hs[i][j] <<"\n";


    if (!steady && !breaking)
    {
        auto temparray = Hs;
        std::for_each(temparray.data(), temparray.data()+temparray.num_elements(), [&](value_type& f){ f *= f; });
        E_tot = std::accumulate(temparray.data(), temparray.data()+temparray.num_elements(),0.);

        // std::fill( var_strain.data(), var_strain.data()+var_strain.num_elements(), 1. );
        // E_tot = std::accumulate(var_strain.data(), var_strain.data()+var_strain.num_elements(),0.);
        // std::cout<<"Sum= "<< E_tot <<"\n";
    }

    // else if (scatmod == "isotropic")
    // {
    //   advAttenIsotropic();
    // }

    // finally do floe breaking

    for (size_type i = 0; i < nx; i++)
    {
        for (size_type j = 0; j < ny; j++)
        {
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
                        std::cout<<"JREST= "<< jcrest <<"\n";
                        om1 = 2*PI*freq_vec[jcrest];
                        lam1 = wlng_ice[i][j][jcrest];
                        lam2 = wlng_ice[i][j][jcrest+1];
                        wlng_crest = lam1+(om-om1)*(lam2-lam1)/dom;
                    }

                    Dc = std::max<value_type>(dmin,wlng_crest/2.0);
                    dfloe[i][j] = std::min<value_type>(Dc,dfloe[i][j]);
                }
            }
            else if (wtr_mask[i][j] == 1.)
                dfloe[i][j] = 0;
        }
    }




#if 0
    for (size_type i = 0; i < nx; i++)
    {
        for (size_type j = 0; j < ny; j++)
        {
            std::cout << "BMWD[" << i << "," << j << "]= " << mwd[i][j] <<"\n";
        }
    }

    calcMWD();

    for (size_type i = 0; i < nx; i++)
    {
        for (size_type j = 0; j < ny; j++)
        {
            std::cout << "AMWD[" << i << "," << j << "]= " << mwd[i][j] <<"\n";
        }
    }
#endif


    for (size_type i = 0; i < nx; i++)
        for (size_type j = 0; j < ny; j++)
            std::cout << "TP[" << i << "," << j << "]= " << Tp[i][j] <<"\n";


}


template<typename T>
void
WimDiscr<T>::wimRun()
{
    value_type x_ext, y_ext, u_ref, duration;
    int nt;
    bool critter;

    x_ext = nx*dx/1.e+03;
    y_ext = ny*dy/1.e+03;
    u_ref = amin + 0.7*(amax-amin);
    duration = 1.0e3*x_ext/u_ref;

    std::cout<<"x_ext= "<< x_ext <<"\n";
    std::cout<<"y_ext= "<< y_ext <<"\n";
    std::cout<<"u_ref= "<< u_ref <<"\n";
    std::cout<<"duration= "<< duration <<"\n";
    std::cout<<"dt= "<< dt <<"\n";

    nt = std::floor(duration/dt);

    //std::cout<<"nt= "<< nt <<"\n";

    size_type cpt = 0;
    nt = 2;

    while (cpt < nt)
    {
        std::cout <<  ":[WIM2D TIME STEP]^"<< cpt+1 <<"\n";

        critter = !(cpt % vm["reps"].template as<int>()) && (vm["checkprog"].template as<bool>());
        if (critter)
            writeFile(cpt);


        wimStep();

        //critter = !(cpt % vm["reps"].template as<int>()) && (vm["checkprog"].template as<bool>());
        // if (critter)
        //     writeFile(cpt);

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

	std::fill( uwave.data(), uwave.data() + uwave.num_elements(), 0. );
	std::fill( vwave.data(), vwave.data() + vwave.num_elements(), 0. );

	std::vector<value_type> wt_theta(nwavedirn);
	value_type adv_dir, S_th, tmp, alp_dim, source;

	for (size_type k = 0; k < nwavedirn; k++)
    {
        adv_dir = -PI*(90.0+wavedir[k])/180.0;

        uwave = ag2d_eff;
        atten_nond.data() + atten_nond.num_elements();
        //std::for_each(uwave.begin(), uwave.end(), [&](value_type& f){ f *= std::cos(adv_dir); });
        std::for_each(uwave.data(), uwave.data()+uwave.num_elements(), [&](value_type& f){ f *= std::cos(adv_dir); });

        if (advdim == 2)
        {
            vwave = ag2d_eff;
            std::for_each(vwave.data(), vwave.data()+vwave.num_elements(), [&](value_type& f){ f *= std::sin(adv_dir); });
        }

        // copy from 3D input array to 2D temporary array
        for (size_type i = 0; i < nx; i++)
        {
            for (size_type j = 0; j < ny; j++)
            {
                temp[i][j] = Sdir[i][j][k];
            }
        }

        // advection
        waveAdvWeno(temp,uwave,vwave);

        // copy from 2D temporary array to 3D input array
        for (size_type i = 0; i < nx; i++)
        {
            for (size_type j = 0; j < ny; j++)
            {
                Sdir[i][j][k] = temp[i][j];
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

    // for (size_type i = 0; i < nx; i++)
    // {
    //     for (size_type j = 0; j < ny; j++)
    //     {
    //         std::cout<<"MASK= "<< ice_mask[i][j] <<"\n";
    //     }
    // }

    for (size_type i = 0; i < nx; i++)
    {
        for (size_type j = 0; j < ny; j++)
        {

            if (ice_mask[i][j] > 0.)
                for (size_type wnd = 0; wnd < nwavedirn; wnd++)
                {
                    adv_dir = -PI*(90.0+wavedir[wnd])/180.0;
                    S_th = Sdir[i][j][wnd];
                    //std::cout<<"S_th= "<< S_th <<"\n";
                    alp_dim = atten_dim[i][j]+damp_dim[i][j];
                    //std::cout<<"alp_dim= "<< atten_dim[i][j] <<"\n";
                    //std::cout<<"alp_dim= "<< damp_dim[i][j] <<"\n";

                    // stress calculation
                    source = -alp_dim*ag2d_eff[i][j]*S_th;
                    tmp = -std::cos(adv_dir)*wt_theta[wnd]*source;
                    taux_omega[i][j] += tmp;
                    tmp = -std::sin(adv_dir)*wt_theta[wnd]*source;
                    tauy_omega[i][j] += tmp;

                    // do attenuation
                    Sdir[i][j][wnd] = S_th*std::exp(-alp_dim*ag2d_eff[i][j]*dt);
                }

            // integrate spectrum over direction
            for (size_type wnd = 0; wnd < nwavedirn; wnd++)
            {
                Sfreq[i][j] += wt_theta[wnd]*Sdir[i][j][wnd];
            }
        }
    }

    // for (size_type k = 0; k < nwavedirn; k++)
    //     for (size_type i = 0; i < nx; i++)
    //     {
    //         for (size_type j = 0; j < ny; j++)
    //         {
    //             std::cout<<"MASK= "<< Sdir[i][j][k] <<"\n";
    //         }
    //     }



}


template<typename T>
void
WimDiscr<T>::advAttenIsotropic(array3_type& Sdir, array2_type& Sfreq, array2_type& taux_omega, array2_type& tauy_omega, array2_type const& ag2d_eff)
{
	std::fill( uwave.data(), uwave.data() + uwave.num_elements(), 0. );
	std::fill( vwave.data(), vwave.data() + vwave.num_elements(), 0. );

	std::vector<int> nvec(nwavedirn);
	std::vector<value_type> K_fou(nwavedirn), S_th(nwavedirn), theta_vec(nwavedirn), wt_theta(nwavedirn);
	std::vector<value_type> tmp1(nwavedirn), evals_x(nwavedirn);
	value_type adv_dir, tmp, alp_dim, source;

	std::vector<std::complex<value_type> > S_fou;
	std::complex<value_type> zi, src_fou_p1, src_fou_m1;

	std::vector<value_type> S_cos, S_sin;
	value_type cg, q_scat, q_abs, q_tot, src_cos_1, src_sin_1;

	zi = std::complex<value_type>(0,1);
	std::cout<<"zi= " << zi <<"\n";

	std::fill( wt_theta.begin(), wt_theta.end(), 2*PI/(1.0*nwavedirn) );

	for (size_type k = 0; k < nwavedirn; k++)
    {
        adv_dir = -PI*(90.0+wavedir[k])/180.0;

        uwave = ag2d_eff;
        atten_nond.data() + atten_nond.num_elements();

        std::for_each(uwave.data(), uwave.data()+uwave.num_elements(), [&](value_type& f){ f *= std::cos(adv_dir); });

        if (advdim == 2)
        {
	        vwave = ag2d_eff;
	        std::for_each(vwave.data(), vwave.data()+vwave.num_elements(), [&](value_type& f){ f *= std::sin(adv_dir); });
        }

        // copy from 3D input array to 2D temporary array
        for (size_type i = 0; i < nx; i++)
        {
	        for (size_type j = 0; j < ny; j++)
	        {
		        temp[i][j] = Sdir[i][j][k];
	        }
        }

        // advection
        waveAdvWeno(temp,uwave,vwave);

        // copy from 2D temporary array to 3D input array
        for (size_type i = 0; i < nx; i++)
        {
	        for (size_type j = 0; j < ny; j++)
	        {
		        Sdir[i][j][k] = temp[i][j];
	        }
        }

        theta_vec[k] = adv_dir;
        nvec[k] = k;
    }

	for (size_type i = 0; i < nx; i++)
    {
	    for (size_type j = 0; j < ny; j++)
	    {
		    for (size_type k = 0; k < nwavedirn; k++)
			    S_th[k] = Sdir[i][j][k];

		    std::fill( S_fou.begin(), S_fou.end(), zi );

		    //S_fou[0] = std::complex<value_type>( sum(wt_theta*S_th) )
	    }
    }


}

template<typename T>
void
WimDiscr<T>::waveAdvWeno(array2_type& h, array2_type const& u, array2_type const& v)
{

	std::fill( sao.data(), sao.data() + sao.num_elements(), 0. );

	padVar(u, u_pad);
	padVar(v, v_pad);
	padVar(SCP2_array, scp2_pad);
    padVar(SCP2I_array, scp2i_pad);
    padVar(SCUY_array, scuy_pad);
    padVar(SCVX_array, scvx_pad);
    padVar(h, h_pad);

    // prediction step
    weno3pdV2(h_pad, u_pad, v_pad, scuy_pad, scvx_pad, scp2i_pad, scp2_pad, sao);



    for (size_type i = 0; i < nxext; i++)
    {
        for (size_type j = 0; j < nyext; j++)
        {
            hp[i][j] = h_pad[i][j]+dt*sao[i][j];
        }
    }

    // correction step
    weno3pdV2(hp, u_pad, v_pad, scuy_pad, scvx_pad, scp2i_pad, scp2_pad, sao);



    for (size_type i = 0; i < nxext; i++)
    {
        for (size_type j = 0; j < nyext; j++)
        {
            h_pad[i][j] = 0.5*(h_pad[i][j]+hp[i][j]+dt*sao[i][j]);
        }
    }


    for (size_type i = nbdy; i < nx+nbdy; i++)
    {
        for (size_type j = nbdy; j < ny+nbdy; j++)
        {
            h[i-nbdy][j-nbdy] = h_pad[i][j];

            // mask land (no waves on land)
            h[i-nbdy][j-nbdy] = h[i-nbdy][j-nbdy]*(1-LANDMASK_array[i-nbdy][j-nbdy]);
        }
    }
}

template<typename T>
void
WimDiscr<T>::weno3pdV2(array2_type const& gin, array2_type const& u, array2_type const& v, array2_type const& scuy,
                       array2_type const& scvx, array2_type const& scp2i, array2_type const& scp2, array2_type& saoout)
{

	value_type cq00=-1./2 ,cq01=3./2, cq10=1./2, cq11=1./2, ca0=1./3, ca1=2./3, eps=1e-12;
	value_type q0, q1, a0, a1, q;
	size_type im1, im2, ip1, jm1, jm2, jp1;

    // fluxes in x direction
    for (size_type i = nbdy-1; i < nx+nbdy+2; i++)
    {
        for (size_type j = nbdy-1; j < ny+nbdy+1; j++)
        {
            im1 = i-1;

            if (u[i][j] > 0.)
            {
                im2 = im1-1;
                q0 = cq00*gin[im2][j]+cq01*gin[im1][j];
                q1 = cq10*gin[im1][j]+cq11*gin[i][j];
                a0 = ca0;
                a1 = ca1*(std::abs(gin[im2][j]-gin[im1][j])+eps)/(std::abs(gin[im1][j]-gin[i][j])+eps);
                ful[i][j] = u[i][j]*gin[im1][j]*scuy[i][j];

            }
            else
            {
                ip1 = i+1;
                q0 = cq11*gin[im1][j]+cq10*gin[i][j];
                q1 = cq01*gin[i][j]+cq00*gin[ip1][j];
                a0 = ca1;
                a1 = ca0*(abs(gin[im1][j]-gin[i][j])+eps)/(abs(gin[i][j]-gin[ip1][j])+eps);
                ful[i][j] = u[i][j]*gin[i][j]*scuy[i][j];
            }

            fuh[i][j] = u[i][j]*(a0*q0+a1*q1)/(a0+a1)*scuy[i][j]-ful[i][j];
        }
    }


    // fluxes in y direction
    for (size_type i = nbdy-1; i < nx+nbdy+1; i++)
    {
        for (size_type j = nbdy-1; j < ny+nbdy+2; j++)
        {
            jm1 = j-1;

            if (v[i][j] > 0.)
            {
                jm2 = jm1-1;
                q0 = cq00*gin[i][jm2]+cq01*gin[i][jm1];
                q1 = cq10*gin[i][jm1]+cq11*gin[i][j];
                a0 = ca0;
                a1 = ca1*(std::abs(gin[i][jm2]-gin[i][jm1])+eps)/(std::abs(gin[i][jm1]-gin[i][j])+eps);
                fvl[i][j] = u[i][j]*gin[i][jm1]*scvx[i][j];
            }
            else
            {
                jp1 = j+1;
                q0 = cq11*gin[i][jm1]+cq10*gin[i][j];
                q1 = cq01*gin[i][j]+cq00*gin[i][jp1];
                a0 = ca1;
                a1 = ca0*(abs(gin[i][jm1]-gin[i][j])+eps)/(abs(gin[i][j]-gin[i][jp1])+eps);
                fvl[i][j] = u[i][j]*gin[i][j]*scvx[i][j];
            }

            fvh[i][j] = u[i][j]*(a0*q0+a1*q1)/(a0+a1)*scvx[i][j]-fvl[i][j];
        }
    }


    // update field with low order fluxes
    for (size_type i = nbdy-1; i < nx+nbdy+1; i++)
    {
        for (size_type j = nbdy-1; j < ny+nbdy+1; j++)
        {
            gt[i][j] = gin[i][j]-dt*(ful[i+1][j]-ful[i][j]+fvl[i][j+1]-fvl[i][j])*scp2i[i][j];
        }
    }

    q = 0.25/dt;

    for (size_type i = nbdy; i < nx+nbdy+1; i++)
    {
        for (size_type j = nbdy; j < ny+nbdy; j++)
        {
            fuh[i][j] = ful[i][j]+std::max(-q*gt[i][j]*scp2[i][j],std::min(q*gt[i-1][j]*scp2[i-1][j],fuh[i][j]));
        }
    }



    for (size_type i = nbdy; i < nx+nbdy; i++)
    {
        for (size_type j = nbdy; j < ny+nbdy+1; j++)
        {
            fvh[i][j]=fvl[i][j]+std::max(-q*gt[i][j]*scp2[i][j],std::min(q*gt[i][j-1]*scp2[i][j-1],fvh[i][j]));
        }
    }



    for (size_type i = nbdy; i < nx+nbdy; i++)
    {
        for (size_type j = nbdy; j < ny+nbdy; j++)
        {
	        saoout[i][j] = -(fuh[i+1][j]-fuh[i][j]+fvh[i][j+1]-fvh[i][j])*scp2i[i][j];
        }
    }

}


template<typename T>
void
WimDiscr<T>::padVar(array2_type const& u, array2_type& upad)
{

	std::fill( upad.data(), upad.data() + upad.num_elements(), 0. );

    for (size_type i = 0; i < nxext; i++)
    {
        for (size_type j = 0; j < nyext; j++)
        {

            if ((nbdy-1 < i) && (i < nx+nbdy) && (nbdy-1 < j) && (j < ny+nbdy))
            {
                upad[i][j] = u[i-nbdy][j-nbdy];
            }


            if (advopt != "notperiodic")
            {
                // make periodic in j
                if ((j < nbdy) && (nbdy-1 < i) && (i < nx+nbdy))
                    upad[i][j] = u[i-nbdy][ny-nbdy+j];
            }

            if (advopt == "xy-periodic")
            {

                // make periodic in i
                if ((i < nbdy) && (nbdy-1 < j) && (j < ny+nbdy))
                    upad[i][j] = u[nx-nbdy+i][j-nbdy];


                if ((nx+nbdy-1 < i) && (nbdy-1 < j) && (j < ny+nbdy))
                    upad[i][j] = u[i-nx-nbdy][j-nbdy];


                // // make periodic in j
                // if ((j < nbdy) && (nbdy-1 < i) && (i < nx+nbdy))
                //     upad[i][j] = u[i-nbdy][ny-nbdy+j];


                if ((ny+nbdy-1 < j) && (nbdy-1 < i) && (i < nx+nbdy))
                    upad[i][j] = u[i-nbdy][j-ny-nbdy];


                // BR, TL
                if ((nx+nbdy-1 < i) && (ny+nbdy-1 < j))
                    upad[i][j] = u[i-nx-nbdy][j-ny-nbdy];


                if ((i < nbdy) && (j < nbdy))
                    upad[i][j] = u[i+nx-nbdy][j];


                // BL, TR
                if ((nx+nbdy-1 < i) && (j < nbdy))
                    upad[i][j] = u[i-nx-nbdy][ny-nbdy+j];


                if ((i < nbdy) && (ny+nbdy-1 < j))
                    upad[i][j] = u[i+nx-nbdy][j-ny-nbdy];

            }
        }
    }

    // for (size_type i = 0; i < nxext; i++)
    // {
    //     for (size_type j = 0; j < nyext; j++)
    //     {
    //         std::cout<<"COPY["<< i << "," << j << "]= "<< upad[i][j] <<"\n";
    //     }
    // }

}


template<typename T>
void
WimDiscr<T>::calcMWD()
{
    value_type adv_dir, wt_theta, om;

    if (nwavedirn == 1)
        wt_theta = 1.;
    else
        wt_theta = 2.0*PI/(nwavedirn);

    // spectral moments
    std::fill( cmom0.data(), cmom0.data() + cmom0.num_elements(), 0. );
    std::fill( cmom_dir.data(), cmom_dir.data() + cmom_dir.num_elements(), 0. );

    for (size_type fq = 0; fq < nwavefreq; fq++)
    {
        om = 2*PI*freq_vec[fq];

        std::fill( CSfreq.data(), CSfreq.data() + CSfreq.num_elements(), 0. );
        std::fill( cmom_dir0.data(), cmom_dir0.data() + cmom_dir0.num_elements(), 0. );

        for (size_type dn = 0; dn < nwavedirn; dn++)
        {
            adv_dir = -PI*(90.0+wavedir[dn])/180.0;

            for (size_type i = 0; i < nx; i++)
            {
                for (size_type j = 0; j < ny; j++)
                {
                    CSfreq[i][j] += wt_theta*sdf_dir[i][j][dn][fq];
                    cmom_dir0[i][j] += wt_theta*sdf_dir[i][j][dn][fq]+adv_dir;
                }
            }
        }


        for (size_type i = 0; i < nx; i++)
        {
            for (size_type j = 0; j < ny; j++)
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

    for (size_type i = 0; i < nx; i++)
    {
        for (size_type j = 0; j < ny; j++)
        {
            if (cmom0[i][j] > 0.)
                mwd[i][j] = -90.-180*(cmom_dir[i][j]/cmom0[i][j])/PI;
        }
    }

}

template class WimDiscr<float>;

} // namespace WIM2D
