/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

#include "boost/program_options.hpp"

namespace po = boost::program_options;

namespace WIMOPT
{
    po::options_description
    descrOptions()
    {
        po::options_description desc("Options");

        desc.add_options()
            ("help", "Print help messages")
            ("config-file", po::value<std::string>(), "specify .cfg file")
            ("wim.nx", po::value<int>()->default_value( 150 ), "Record length in x direction")
            ("wim.ny", po::value<int>()->default_value( 10 ), "Record length in y direction")
            ("wim.dx", po::value<double>()->default_value( 4e+3 ), "Resolution in x direction [m]")
            ("wim.dy", po::value<double>()->default_value( 4e+3 ), "Resolution in y direction [m]")
            ("wim.xmin", po::value<double>()->default_value( -298.e+3 ), "xmin [m]")
            ("wim.ymin", po::value<double>()->default_value( -60e+3 ), "ymin [m]")
            ("wim.nwavefreq", po::value<int>()->default_value( 1 ), "Number of wave frequencies")
            ("wim.nwavedirn", po::value<int>()->default_value( 16 ), "Number of wave directions")

            //'int_prams' in fortran
            ("wim.scatmod", po::value<std::string>()->default_value( "dissipated" ), "Scattered energy is dissipated (=dissipated), distributed isotropically (=isotropic)")
            ("wim.advopt", po::value<std::string>()->default_value( "y-periodic" ), "Not periodic (=notperiodic), periodic in y only (=y-periodic), periodic in both x,y (=xy-periodic)")
            ("wim.advdim", po::value<int>()->default_value( 2 ), "Dimension of advection scheme (1 or 2)")
            ("wim.steady", po::value<bool>()->default_value( true ), "Steady-state (=true), or not steady-state (=false)")
            ("wim.breaking", po::value<bool>()->default_value( true ), "Do breaking (=true), or turn off breaking (=false)")
            ("wim.atten", po::value<bool>()->default_value( true ), "Do attenuation")
            ("wim.checkprog", po::value<bool>()->default_value( true ), "Do dump intermediate states to binary files (=true), or don't' (=false)")

            //other bool param's
            ("wim.refhsice", po::value<bool>()->default_value( false ), "Inside ice, Hs corresponds to water (=false) or ice (=true) displacement")
            ("wim.useicevel", po::value<bool>()->default_value( false ), "Inside ice, use correct group velocity (=true), or water group velocity (=false)")

            // 'real_prams' in fortran code
            ("wim.young", po::value<double>()->default_value( 5.49e+9 ), "Young's modulus [Pa]")
            ("wim.viscrp", po::value<double>()->default_value( 13. ), "Robinson-Palmer viscosity [Pa.s/m]")
            ("wim.duration", po::value<double>()->default_value( 43200.0 ), "length of simulation [s]")
            ("wim.cfl", po::value<double>()->default_value( 0.7 ), "CFL number")

            //initial conditions in idealised simulation
            ("wim.hsinc", po::value<double>()->default_value( 3. ), "Incident significant wave height [m]")
            ("wim.tpinc", po::value<double>()->default_value( 12. ), "Incident peak period [s]")
            ("wim.mwdinc", po::value<double>()->default_value( -90. ), "Incident mean wave-from direction [deg]")
            ("wim.unifc", po::value<double>()->default_value( 0.7 ), "Initial const conc")
            ("wim.unifh", po::value<double>()->default_value( 1. ), "Initial const thickness [m]")
            ("wim.dfloepackinit", po::value<double>()->default_value( 300. ), "Initial value in pack (unbroken) ice [m]")
            ("wim.landon3edges", po::value<bool>()->default_value( false ), "Add land on upper,lower and RH edges")

            //outputs of WIM
            ("wim.dumpfreq", po::value<int>()->default_value( 10 ), "frequency of dumping (# WIM timesteps)")
            ("wim.outparentdir", po::value<std::string>()->default_value( "out_cpp" ), "Parent directory for the output files")

            //numerical parameters
            ("wim.tmin", po::value<double>()->default_value( 2.5 ), "Minimum wave period in a spectrum [s]")
            ("wim.tmax", po::value<double>()->default_value( 25. ), "Maximum wave period in a spectrum [s]")
            ("wim.dfloemin", po::value<double>()->default_value( 20. ), "Minimum floe size [m]")
            ("wim.cicemin", po::value<double>()->default_value( 0.05 ), "Minimum ice conc considered by WIM")
            ("wim.dfloepackthresh", po::value<double>()->default_value( 400. ), "Don't let Dmax grow above this value [m]")
            ("wim.nbdy", po::value<int>()->default_value( 3 ), "Size of the overlap for periodic boundary conditions")
            ("wim.nbdx", po::value<int>()->default_value( 3 ), "Size of the overlap for periodic boundary conditions")

            //coupling to nextsim
            ("wim.exportresults", po::value<bool>()->default_value( true ), "Export results in coupled mode")
            ("wim.nfloesgridtomesh", po::value<bool>()->default_value( true ), "During neXtSIM regridding interpolate from grid-to-mesh or mesh-to-mesh")
            ("wim.couplingfreq", po::value<int>()->default_value( 20 ), "Coupling frequency between neXtSIM and WIM (# neXtSIM time-steps)")
            ;
        return desc;
    }

} // WIMOPT