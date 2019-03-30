/**
    io.hpp
    I/O operations for deconvms and wassms
    
    @author Antoine Passemiers
    @version 1.0 30/03/2019
*/

#ifndef IO_HPP__
#define IO_HPP__

#include <iostream>
#include <fstream>
#include <sstream>

#include "spectrum.hpp"

namespace wassersteinms {

/**
    Loads a mass spectrum from text file.

    @param filepath Location of the text file.
    @return A mass spectrum
*/
Spectrum loadRecord(const std::string filepath);

} // wassersteinms

#endif // IO_HPP__