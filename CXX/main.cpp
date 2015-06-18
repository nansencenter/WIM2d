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

    //wim2d.readFile("wim_grid.a");

    // generation and saving of the grid
    //wim2d.wimGrid();

    // auto X = wim2d.getX();
    // for (int i = 0; i < X.shape()[0]; i++)
    //     for (int j = 0; j < X.shape()[1]; j++)
    //         std::cout << "X[" << i << "," << j << "]= " << X[i][j] << std::endl;

    // initialization of wim2d
    wim2d.wimInit();

    //wim2d.readFile("wim_prog000.a");

    // auto X = wim2d.getY();
    // for (int i = 0; i < X.shape()[0]; i++)
    //     for (int j = 0; j < X.shape()[1]; j++)
    //         std::cout << "X[" << i << "," << j << "]= " << X[i][j] << std::endl;


    //wim2d.advAttenIsotropic();

    //wim2d.writeFile(0);

    // run the simulation
    wim2d.wimRun();
    //wim2d.wimStep();

}
