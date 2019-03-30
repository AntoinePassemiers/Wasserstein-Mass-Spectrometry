/**
    deconvms.cpp
    Deconvolve empirical spectra
    
    @author Antoine Passemiers
    @version 1.0 30/03/2019
*/

#include <algorithm>
#include <iostream>
#include <cassert>
#include <cstddef>
#include <cstring>
#include <memory>
#include <sstream>

#include "spectrum.hpp"
#include "io.hpp"
#include "ipm.hpp"

using namespace wassersteinms;

typedef struct _params {
    char *filepath1;
    char *filepath2;
    char *folder;
    float epsilon;
    size_t nMaxIterations;
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
    pars.nMaxIterations = 10;

    for (int i = 3; i < argc; i++) {
        if (strcmp(argv[i], "--eps") == 0) {
            pars.epsilon = atof(argv[++i]);
        } else if (strcmp(argv[i], "--niter") == 0) {
	    pars.nMaxIterations = atof(argv[++i]);
	}
    }
    return pars;
}


int main(int argc, char *argv[]) {

    params pars = parseCLA(argc, argv);
    if (pars.hasParseError) return 1;

    Spectrum mixture = loadRecord(pars.filepath1).normalize();

    std::vector<Spectrum> theoreticalSpectra;
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
            Spectrum spectrum = loadRecord(filepath).normalize();
            theoreticalSpectra.push_back(spectrum);
            mixture.addKeys(spectrum);
        }
    }

    for (auto it = theoreticalSpectra.begin(); it != theoreticalSpectra.end(); it++) {
        Spectrum &spectrum = *it;
        spectrum.addKeys(mixture);
        std::cout << mixture.size() << ", " << spectrum.size() << std::endl;
        assert(mixture.size() == spectrum.size());
    }
    
    std::unique_ptr<ProblemInstance> problemInstance = formulateProblem(theoreticalSpectra, mixture);
    size_t k = problemInstance->k;

    std::unique_ptr<IpmSolution> sol = longStepPathFollowingMethod(
		    problemInstance, pars.epsilon, pars.nMaxIterations);
    Eigen::VectorXd p = sol->y.tail(k);
    for (size_t i = 0; i < k; i++) {
        std::cout << "Weight of isotopic enveloppe " << i+1 << ": " << p[i] << std::endl;
    }

    std::cout << "Finished";
    return 0;
}
