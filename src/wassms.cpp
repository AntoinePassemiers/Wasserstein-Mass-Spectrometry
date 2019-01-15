#include "spectrum.hpp"

#include <algorithm>
#include <iostream>
#include <cstddef>
#include <fstream>
#include <sstream>




void Spectrum::sort() {
    std::vector<int> indices(ratios.size());
    std::size_t n = 0;
    std::generate(std::begin(indices), std::end(indices), [&]{ return n++; });
    std::sort(
            std::begin(indices),
            std::end(indices),
            [&](int idx1, int idx2) { return ratios[idx1] < ratios[idx2]; } );
    std::vector<mz_t> newRatios(ratios.size());
    std::vector<intensity_t> newIntensities(ratios.size());
    std::vector<intensity_t> newRelIntensities(ratios.size());
    for (int i = 0; i < ratios.size(); i++) {
        newRatios[i] = ratios[indices[i]];
        newIntensities[i] = intensities[indices[i]];
        newRelIntensities[i] = relIntensities[indices[i]];
    }
    ratios = newRatios;
    intensities = newIntensities;
    relIntensities = newRelIntensities;
    sorted = true;
}


int main(int argc, char *argv[]) {

    std::string filepath;
    if (argc > 1) {
        filepath = argv[1];
    } else {
        filepath = "FIO00010.txt";
    }
    std::ifstream record_file(filepath);


    Spectrum *spectrum = new Spectrum();
    std::string line;
    while (std::getline(record_file, line)) {
        float mz, intensity, relIntensity;
        std::istringstream iss(line);

        if (!(iss >> mz >> intensity >> relIntensity)) {
            // Do nothing
        } else {
            spectrum->add(mz, intensity, relIntensity);
        }

    }
    spectrum->sort();


    std::cout << "Finished";
    return 0;
}
