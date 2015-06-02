/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

#include "wimdiscr.hpp"

namespace
{
    const size_t ERROR_IN_COMMAND_LINE = 1;
    const size_t SUCCESS = 0;
} // namespace

// BOOST_STRONG_TYPEDEF(unsigned, Unsigned)

// namespace po = boost::program_options;

namespace WIMOPT
{
    po::options_description descrOptions();
} // WIMOPT

int main(int argc, char** argv )
{

    using namespace WIMOPT;

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

    std::cout<<"nx= "<< vm["nx"]. template as<int>() <<"\n";
    std::cout<<"ny= "<< vm["ny"]. template as<int>() <<"\n";
    std::cout<<"wim.nz= "<< vm["wim.nz"]. template as<int>() <<"\n";

    WimDiscr<> wim2d(vm);
    wim2d.readFile("wim_grid.a",150,4);

    auto X = wim2d.getX();

    for (int i = 0; i < X.size(); i++)
       for (int j = 0; j < X[i].size(); j++)
          std::cout << "X[" << i << "," << j << "]= " << X[i][j] << std::endl;



    auto Y = wim2d.getY();

    for (int i = 0; i < Y.size(); i++)
       for (int j = 0; j < Y[i].size(); j++)
          std::cout << "Y[" << i << "," << j << "]= " << Y[i][j] << std::endl;



    std::cout<<"PI= "<< float(PI) <<"\n";

    wim2d.writeFile("grid.bin");

    // test init
    wim2d.wimInit();

    // test step
    wim2d.wimStep();

}
