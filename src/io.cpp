#include "io.hpp"


Spectrum *loadRecord(std::string filepath) {
    std::ifstream recordFile(filepath);
    Spectrum *spectrum = new Spectrum();
    std::string line;
    while (std::getline(recordFile, line)) {
        float mz, intensity, relIntensity;
        std::istringstream iss(line);

        if (!(iss >> mz >> intensity >> relIntensity)) {
            // Do nothing
        } else {
            spectrum->add(mz, intensity);
        }
    }
    spectrum->sort();
    spectrum->normalize();
    return spectrum;
}