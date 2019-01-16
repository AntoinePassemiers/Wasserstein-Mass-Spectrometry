#include "spectrum.hpp"
#include "wasserstein.hpp"

#include <algorithm>
#include <iostream>
#include <cstddef>
#include <fstream>
#include <sstream>



Spectrum *loadRecord(std::string filepath) {
    std::ifstream record_file(filepath);
    Spectrum *spectrum = new Spectrum();
    std::string line;
    while (std::getline(record_file, line)) {
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


int main(int argc, char *argv[]) {

    std::string filepath1, filepath2;
    if (argc > 2) {
        filepath1 = argv[1];
        filepath2 = argv[2];
    } else {
        filepath1 = "FIO00010.txt";
        filepath2 = "MCH00013.txt";
    }

    Spectrum *spectrum1 = loadRecord(filepath1);
    Spectrum *spectrum2 = loadRecord(filepath2);

    std::cout << "Distance: " << wassersteinDistance(spectrum1, spectrum2) << std::endl;

    std::cout << "Finished";
    return 0;
}
