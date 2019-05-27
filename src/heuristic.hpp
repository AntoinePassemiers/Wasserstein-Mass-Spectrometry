/**
    heuristic.cpp
    Heuristic method for mass spectrum deconvolution
    
    @author Antoine Passemiers
    @version 1.0 27/05/2019
*/

#ifndef HEURISTIC_HPP__
#define HEURISTIC_HPP__

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <memory>
#include <random>
#include <functional>
#include <math.h>
#include <Eigen/Core>
#include <Eigen/Dense>

#include "ipm.hpp"

namespace wassersteinms {

void applyCorrection(Eigen::VectorXd &p);

double wassersteinDistance(Eigen::VectorXd &p,
                           std::unique_ptr<ProblemInstance> &prob);

Eigen::VectorXd geneticAlgorithm(
        std::unique_ptr<ProblemInstance> &prob,
        double epsilon,
        size_t nMaxIterations);

} // namespace wassersteinms

#endif // HEURISTIC_HPP__