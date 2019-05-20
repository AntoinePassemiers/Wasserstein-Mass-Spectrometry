/**
    io.cpp
    I/O operations for deconvms and wassms
    
    @author Antoine Passemiers
    @version 1.0 30/03/2019
*/

#include "io.hpp"

namespace wassersteinms {

// Load a mass spectrum from text file
Spectrum loadRecord(const std::string filepath) {
    std::ifstream recordFile(filepath);
    if (recordFile.fail()) {
        // TODO: throw error
    }
    Spectrum spectrum;
    std::string line;
    while (std::getline(recordFile, line)) {
        double mz, intensity;
        std::istringstream iss(line);

        if (!(iss >> mz >> intensity)) {
            // Do nothing
        } else {
            spectrum[mz] = intensity;
        }
    }
    std::cout << "Size: " << spectrum.size() << std::endl;
    return spectrum;
}

} // namespace wassersteinms
