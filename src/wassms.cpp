/**
    wassms.cpp
    Compute dissimilarity between two spectra
    
    @author Antoine Passemiers
    @version 1.0 30/03/2019
*/

#include <algorithm>
#include <iostream>
#include <cstddef>
#include <cstring>
#include <memory>

#include "spectrum.hpp"
#include "io.hpp"
#include "similarity.hpp"

using namespace wassersteinms;

typedef struct _params {
    char *filepath1      = nullptr;
    char *filepath2      = nullptr;
    Similarities method  = Similarities::WASSERSTEIN;
    double resolution    = 0.02;
    bool has_parse_error = 0;
} params;


params parseError(params pars) {
    std::cout << "Error. Calls to wassms must be of the form:\n" << std::endl;
    std::cout << "\twassms <record_file_1> <record_file_2> [--m method] [--r resolution]\n" << std::endl;
    std::cout << "Method can be either 'W' (Wasserstein), 'E' (euclidean) ";
    std::cout << "or 'J' (Jaccard score). Resolution is expressed in Daltons." << std::endl;
    pars.has_parse_error = 1;
    return pars;
}


params parseCLA(int argc, char *argv[]) {
    // Initialize parameters
    params pars;

    if (argc < 3) return parseError(pars);

    pars.filepath1 = argv[1];
    pars.filepath2 = argv[2];
    for (int i = 3; i < argc; i++) {
        if (strcmp(argv[i], "--m") == 0) {
            ++i;
            if (strcmp(argv[i], "E") == 0) {
                pars.method = Similarities::EUCLIDEAN;
            } else if (strcmp(argv[i], "J") == 0) {
                pars.method = Similarities::JACCARD_SCORE;
            } else if (strcmp(argv[i], "W") == 0) {
                pars.method = Similarities::WASSERSTEIN;
            } else {
                std::cout << "Error. Unknown method '" << argv[i];
                std::cout << "'." << std::endl;
            }
        } else if (strcmp(argv[i], "--r") == 0) {
            pars.resolution = atof(argv[++i]);
        }
    }
    return pars;
}


int main(int argc, char *argv[]) {

    params pars = parseCLA(argc, argv);
    if (pars.has_parse_error) return 1;

    Spectrum spectrum1 = loadRecord(pars.filepath1);
    Spectrum spectrum2 = loadRecord(pars.filepath2);

    if (pars.method == Similarities::WASSERSTEIN) {
        spectrum1 = spectrum1.normalize();
        spectrum2 = spectrum2.normalize();
        std::cout << "Wasserstein distance: ";
        std::cout << wassersteinDistance(spectrum1, spectrum2) << std::endl;
    } else if (pars.method == Similarities::EUCLIDEAN) {
        std::cout << "Euclidean distance: ";
        std::cout << euclideanDistance(spectrum1, spectrum2, pars.resolution);
        std::cout << std::endl;
    } else if (pars.method == Similarities::JACCARD_SCORE) {
        std::cout << "Jaccard score: ";
        std::cout << jaccardScore(spectrum1, spectrum2, pars.resolution, 0.05);
        std::cout << std::endl;
    }

    std::cout << "Finished";
    return 0;
}
