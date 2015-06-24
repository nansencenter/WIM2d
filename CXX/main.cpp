/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

#include <wimdiscr.hpp>

namespace
{
    const size_t ERROR_IN_COMMAND_LINE = 1;
    const size_t SUCCESS = 0;
} // namespace

namespace po = boost::program_options;

namespace WIMOPT
{
    po::options_description descrOptions();
}// WIMOPT

int main(int argc, char** argv )
{

    using namespace WIMOPT;
    using namespace WIM2D;

    po::options_description desc = descrOptions();
    po::variables_map vm;

    try
    {
        po::store(po::parse_command_line(argc, argv, desc),vm);

        if ( vm.count("help")  )
        {
            std::cout << "Basic Command Line Parameter Application" <<"\n"
                      << desc << "\n";
            return SUCCESS;
        }

        if ( vm.count( "config-file" ) )
        {
            if ( fs::exists( vm["config-file"].as<std::string>() ) )
            {
                std::ifstream ifs( vm["config-file"].as<std::string>().c_str() );
                po::store( parse_config_file( ifs, desc, true ), vm );
                po::notify( vm );
            }
            else
            {
                std::cout << "Cannot found " << "config-file" << "\n";
            }
        }

        po::notify(vm);

    }
    catch(po::error& e)
    {
        std::cerr << "ERROR: " << e.what() << std::endl <<"\n";
        std::cerr << desc <<"\n";
        return ERROR_IN_COMMAND_LINE;
    }

    // instantiation of wim2d
    WimDiscr<float> wim2d(vm);

    // initialization of wim2d
    wim2d.wimInit();

    // run the simulation
    wim2d.wimRun();

#if 0
    //std::vector<float> vec1 = {1.,2.,3.,4.,5.};
    std::vector<std::complex<float> > vec1 = {(1.,2.),2.,3.,4.,5.};
    std::vector<std::complex<float> > vec2 = {(1.5,1.),4.,7.,1.,2.};
    // float sum_of_elems = std::accumulate(vec1.begin(),vec1.end(),0.);
    //std::transform(vec1.begin(), vec1.end(), vec2.begin(), vec1.begin(), std::multiplies<std::complex<float> >());
    //std::complex<float> ans = (1.,1.);
    //ans = std::accumulate(vec1.begin(), vec1.end(), std::complex<float>{1.0}, std::multiplies<std::complex<float> >{});

    std::transform(vec1.begin(), vec1.end(), vec1.begin(), vec1.begin(), [](std::complex<float> a, std::complex<float> b) -> std::complex<float> { return a+b; });

    //std::cout<<"ANS= "<< ans <<"\n";

    for (auto const& it : vec1)
    {
        std::cout<<"Real= "<< std::real(it) <<"\n";
        std::cout<<"Imag= "<< it.imag() <<"\n";
    }
    //float sum_of_elems = std::accumulate(vec1.begin(),vec1.end(),0.);

    //std::cout<<"sum= "<< sum_of_elems <<"\n";

    std::complex<float> a; //(1,2);
    std::complex<float> b; //(3,5);

    a.real(1);
    a.imag(2);
    b = a;

    std::cout<<"A*B= "<< b <<"\n";

    //using namespace std::literals;
    //std::complex<double> z1 = 1i * 1i;
#endif

}
