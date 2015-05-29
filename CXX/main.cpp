/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

#include "wimdiscr.hpp"

namespace
{
    const size_t ERROR_IN_COMMAND_LINE = 1;
    const size_t SUCCESS = 0;
} // namespace

BOOST_STRONG_TYPEDEF(unsigned, Unsigned)

namespace po = boost::program_options;

int main(int argc, char** argv )
{

    po::options_description desc("Options");

    desc.add_options()
        ("help", "Print help messages")
        ("nx", po::value<int>()->default_value( 150 ), "number of subdomains in X")
        ("ny", po::value<int>()->default_value( 4 ), "number of subdomains in Y")
        ("nz", po::value<int>()->default_value( 1 ), "number of subdomains in Z")
        ;

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
    std::cout<<"nz= "<< vm["nz"]. template as<int>() <<"\n";

   //  // std::cout<<"run starts\n";

    WimDiscr<> wim2d;
    wim2d.readFile("wim_grid.a",150,4);

    auto X = wim2d.getX();

    // for (int i = 0; i < X.size(); i++)
    //    for (int j = 0; j < X[i].size(); j++)
    //       std::cout << "Element[" << i << "," << j << "]= " << X[i][j] << std::endl;

    wim2d.writeFile("grid.bin",150,4);
}
