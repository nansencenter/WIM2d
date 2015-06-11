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
            ("nx", po::value<int>()->default_value( 150 ), "Record length in x direction")
            ("ny", po::value<int>()->default_value( 4 ), "Record length in y direction")
            ("dx", po::value<double>()->default_value( 4e+3 ), "Resolution in x direction")
            ("dy", po::value<double>()->default_value( 4e+4 ), "Resolution in y direction")
            ("nwavefreq", po::value<int>()->default_value( 1 ), "Number of wave frequency")
            ("nwavedirn", po::value<int>()->default_value( 16 ), "Number of wave direction")
            ("cfl", po::value<double>()->default_value( 0.7 ), "CFL number")
            ("atten", po::value<bool>()->default_value( true ), "Do attenuation")
            ("refhsice", po::value<bool>()->default_value( false ), "Inside ice, Hs corresponds to water (=false) or ice (=true) displacement")
            ("icevel", po::value<bool>()->default_value( false ), "Inside ice, use correct group velocity (=true), or water group velocity (=false)")
            ("scatmod", po::value<std::string>()->default_value( "isotropic" ), "Scattered energy is dissipated (=dissipated), distributed isotropically (=isotropic)")
            ("advopt", po::value<std::string>()->default_value( "y-periodic" ), "Not periodic (=notperiodic), periodic in y only (=y-periodic), periodic in both x,y (=xy-periodic)")
            ("advdim", po::value<int>()->default_value( 2 ), "Dimension of advection scheme (1 or 2)")
            ("nbdy", po::value<int>()->default_value( 3 ), "Size of the overlap for periodic boundary conditions")
            ("steady", po::value<bool>()->default_value( true ), "Steady-state (=true), or not steady-state (=false)")
            ("breaking", po::value<bool>()->default_value( true ), "Do breaking (=true), or turn off breaking (=false)")
            ("wim.nz", po::value<int>()->default_value( 1000 ), "Record length in x direction")
            ("young", po::value<double>()->default_value( 5.45e+9 ), "Young's modulus")
            ("viscrp", po::value<double>()->default_value( 13. ), "Robinson-Palmer viscosity")
            ;
        return desc;
    }

} // WIMOPT
