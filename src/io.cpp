#include "io.hpp"


Spectrum loadRecord(const std::string filepath) {
    std::ifstream recordFile(filepath);
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
    return spectrum;
}