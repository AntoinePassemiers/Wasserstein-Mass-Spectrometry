#include "spectrum.hpp"
#include "io.hpp"
#include "wasserstein.hpp"

#include <algorithm>
#include <iostream>
#include <cstddef>
#include <cstring>
#include <memory>
#include <sstream>


typedef struct _params {
    char *filepath1;
    char *filepath2;
    char *folder;
    float threshold;
    bool hasParseError;
} params;


params displayParseError(params pars) {
  std::cout << "Error. Calls to deconvms must be of the form:\n" << std::endl;
  std::cout << "\twassms <mixture_record_file> <molecule_list_file> ";
  std::cout << "<molecules_folder> [--thresh value]" << std::endl;
  pars.hasParseError = 1;
  return pars;
}


params parseCLA(int argc, char *argv[]) {
    params pars;
    memset(&pars, 0x00, sizeof(params));
    if (argc < 4) return displayParseError(pars);

    pars.filepath1 = argv[1];
    pars.filepath2 = argv[2];
    pars.folder = argv[3];

    for (int i = 3; i < argc; i++) {
        if (strcmp(argv[i], "--thresh") == 0) {
            pars.threshold = atof(argv[++i]);
        }
    }
    return pars;
}


int main(int argc, char *argv[]) {

    params pars = parseCLA(argc, argv);
    if (pars.hasParseError) return 1;

    std::unique_ptr<Spectrum> mixture(loadRecord(pars.filepath1));
    std::vector<std::unique_ptr<Spectrum>> theoreticalSpectra;

    std::ifstream recordFile(pars.filepath2);
    std::string line;
    while (std::getline(recordFile, line)) {
        std::string filename;
        std::istringstream iss(line);
        if (!(iss >> filename)) {
            // Do nothing
        } else {
            std::stringstream ss;
            ss << pars.folder << "/" << filename;
            std::string filepath = ss.str();
            theoreticalSpectra.push_back(std::unique_ptr<Spectrum>(loadRecord(filepath)));
        }
    }

    std::cout << "Finished";
    return 0;
}
