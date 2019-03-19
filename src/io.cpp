#include "io.hpp"


Spectrum *loadRecord(std::string filepath) {
    std::ifstream recordFile(filepath);
    Spectrum *spectrum = new Spectrum();
    std::string line;
    while (std::getline(recordFile, line)) {
        double mz, intensity;
        std::istringstream iss(line);

        if (!(iss >> mz >> intensity)) {
            // Do nothing
        } else {
            (*spectrum)[mz] = intensity;
        }
    }
    spectrum->normalize();
    return spectrum;
}