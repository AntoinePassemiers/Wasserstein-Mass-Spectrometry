#include "spectrum.hpp"
#include "io.hpp"
#include "similarity.hpp"

#include <algorithm>
#include <iostream>
#include <cstddef>
#include <cstring>
#include <memory>


typedef struct _params {
    char *filepath1;
    char *filepath2;
    Similarities method;
    double resolution;
    bool has_parse_error;
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
    memset(&pars, 0x00, sizeof(params));
    pars.method = Similarities::WASSERSTEIN;

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
            } else {
                pars.method = Similarities::WASSERSTEIN;
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
        std::cout << "Distance: " << wassersteinDistance(spectrum1, spectrum2) << std::endl;
    }

    std::cout << "Finished";
    return 0;
}
