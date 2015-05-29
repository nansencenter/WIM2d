/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

#ifndef __WIMDISCR_H
#define __WIMDISCR_H 1

#include <iostream>
#include <fstream>

#include <boost/multi_array.hpp>
#include <boost/array.hpp>
#include "boost/program_options.hpp"

#include <boost/any.hpp>
#include <boost/strong_typedef.hpp>
//#include <boost/program_options.hpp>

template<typename T=float> class WimDiscr
{

	typedef T value_type;
    typedef size_t size_type;
	typedef boost::multi_array<float, 2> array_type;

public:

    //WimDiscr () : {}
	void readFile (std::string const& filein, size_type const& nx, size_type const& ny);
    void writeFile (std::string const& fileout, size_type const& nx, size_type const& ny) const;

	//std::fstream __out( M_filename.c_str(), std::ios::out | std::ios::binary );

    std::vector<std::vector<value_type> > getX() const { return X_array; }
    std::vector<std::vector<value_type> > getY() const { return Y_array; }
    std::vector<std::vector<value_type> > getSCUY() const { return SCUY_array; }
    std::vector<std::vector<value_type> > getSCVX() const { return SCVX_array; }
    std::vector<std::vector<value_type> > getSCP2() const { return SCP2_array; }
    std::vector<std::vector<value_type> > getSCP2I() const { return SCP2I_array; }
    std::vector<std::vector<value_type> > getLANDMASK() const { return LANDMASK_array; }

private:

    //size_type nx, ny, recno;
    std::vector<std::vector<value_type> > X_array, Y_array, SCUY_array, SCVX_array, SCP2_array, SCP2I_array, LANDMASK_array;

};

template<typename T>
void
WimDiscr<T>::readFile (std::string const& filein, size_type const& nx, size_type const& ny)
{

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
WimDiscr<T>::writeFile (std::string const& fileout, size_type const& nx, size_type const& ny) const
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

#endif
