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
#include "heuristic.hpp"

using namespace wassersteinms;

typedef struct _params {

    // Location of text file containing mass empirical spectrum
    char *filepath1;

    // Location of text file containing the list of
    // filenames of mass theoretical spectra
    char *filepath2;

    // Location of the folder containing the theoretical
    // spectra text files
    char *folder;

    // Convergence threshold of the interior-point method.
    // The algorithm will stop if the norm of primal residuals
    // and the norm of dual residuals both fall below this
    // threshold epsilon.
    float epsilon = 1e-10;

    // Spectra resolution (may affect performance)
    double resolution = 0.01;

    // Maximum number of iterations of the interior-point method
    size_t nMaxIterations = 10;

    // Presence of an error in command arguments.
    // By default, there is no error.
    bool hasParseError = 0;
} params;


/**
    Displays a parse error in standard output.

    @param pars Parameters struct used to store command line arguments
    @return The same struct with hasParseError flag set to one.
*/
params displayParseError(params pars) {
  std::cout << "Error. Calls to deconvms must be of the form:\n" << std::endl;
  std::cout << "\twassms <mixture_record_file> <molecule_list_file> ";
  std::cout << "<molecules_folder> [--thresh value]" << std::endl;
  pars.hasParseError = 1;
  return pars;
}


/**
    Parses command line arguments and stores them in a struct.

    @param argc The number of arguments
    @param argv Parsed arguments
    @return Parameters struct
*/
params parseCLA(int argc, char *argv[]) {
    params pars;
    if (argc < 4) return displayParseError(pars);

    // Mandatory arguments
    pars.filepath1 = argv[1];
    pars.filepath2 = argv[2];
    pars.folder = argv[3];

    // Optional arguments
    for (int i = 3; i < argc; i++) {
        if (strcmp(argv[i], "--eps") == 0) {
            pars.epsilon = atof(argv[++i]);
        } else if (strcmp(argv[i], "--niter") == 0) {
    	    pars.nMaxIterations = atof(argv[++i]);
    	} else if (strcmp(argv[i], "--res") == 0) {
            pars.resolution = atof(argv[++i]);
        }
    }
    return pars;
}


int main(int argc, char *argv[]) {

    // Parse command line arguments
    params pars = parseCLA(argc, argv);
    if (pars.hasParseError) return 1;

    // Load empirical spectrum from text file
    Spectrum mixture = loadRecord(pars.filepath1) \
            .changeResolution(pars.resolution) \
            .normalize();

    // Load theoretical spectra from text files
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
            Spectrum spectrum = loadRecord(filepath) \
                    .changeResolution(pars.resolution) \
                    .normalize();
            theoreticalSpectra.push_back(spectrum);

            // Adds to empirical spectrum bins that are present only 
            // in theoretical spectra. This step aims at setting empirical spectrum's
            // size to the total number of bins n.
            mixture.addKeys(spectrum);
        }
    }

    // Add to theretical spectra bins that are now only present in empirical spectrum.
    // This step aims at setting theoretical spectra's size to the total number
    // of bins n.
    for (auto it = theoreticalSpectra.begin(); it != theoreticalSpectra.end(); it++) {
        Spectrum &spectrum = *it;
        spectrum.addKeys(mixture);
        assert(mixture.size() == spectrum.size());
    }
    
    // Formulate the spectrum deconvolution problem
    std::unique_ptr<ProblemInstance> problemInstance = formulateProblem(theoreticalSpectra, mixture);
    size_t k = problemInstance->k;

    // Solve deconvolution problem
    std::cout << "Running genetic algorithm for " << pars.nMaxIterations << " iteration(s)..." << std::endl;
    Eigen::VectorXd p = geneticAlgorithm(
		    problemInstance, pars.epsilon, pars.nMaxIterations);
    for (size_t i = 0; i < k; i++) {
        std::cout << "Weight of isotopic enveloppe " << i+1 << ": " << p[i] << std::endl;
    }

    std::cout << "Finished";
    return 0;
}
