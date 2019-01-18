#include "spectrum.hpp"
#include "io.hpp"
#include "wasserstein.hpp"

#include <algorithm>
#include <iostream>
#include <cstddef>
#include <cstring>
#include <memory>


typedef struct _params {
    char *filepath1;
    char *filepath2;
    float threshold;
    bool has_parse_error;
} params;


params parse_error(params pars) {
    std::cout << "Error. Calls to wassms must be of the form:\n" << std::endl;
    std::cout << "\twassms <record_file_1> <record_file_2> [--thresh value]" << std::endl;
    pars.has_parse_error = 1;
    return pars;
}


params parseCLA(int argc, char *argv[]) {
    params pars;
    memset(&pars, 0x00, sizeof(params));
    if (argc < 3) return parse_error(pars);

    pars.filepath1 = argv[1];
    pars.filepath2 = argv[2];

    for (int i = 3 ; i < argc ; i++) {
        if (strcmp(argv[i], "--thresh") == 0) {
            pars.threshold = atof(argv[++i]);
        }
    }
    return pars;
}


int main(int argc, char *argv[]) {

    params pars = parseCLA(argc, argv);
    if (pars.has_parse_error) return 1;

    std::unique_ptr<Spectrum> spectrum1(loadRecord(pars.filepath1));
    std::unique_ptr<Spectrum> spectrum2(loadRecord(pars.filepath2));

    std::cout << "Distance: " << wassersteinDistance(spectrum1, spectrum2) << std::endl;

    std::cout << "Finished";
    return 0;
}
