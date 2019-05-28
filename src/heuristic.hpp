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
#include <vector>
#include <math.h>
#include <Eigen/Core>
#include <Eigen/Dense>

#include "ipm.hpp"

namespace wassersteinms {

int argmin(const Eigen::VectorXd &v);

int argmin(const Eigen::VectorXd &v, const Eigen::VectorXi &indices);

int argmax(const Eigen::VectorXd &v);

void applyCorrection(Eigen::VectorXd &p);

Eigen::VectorXd randomSolution(size_t k);

void mutationOperator(Eigen::VectorXd &p);

void crossoverOperator(const Eigen::VectorXd &a,
                       const Eigen::VectorXd &b,
                       Eigen::VectorXd &c);

double wassersteinDistance(const Eigen::VectorXd &p,
                           std::unique_ptr<ProblemInstance> &prob);

Eigen::VectorXd geneticAlgorithm(
        std::unique_ptr<ProblemInstance> &prob,
        double epsilon,
        size_t nMaxIterations);

} // namespace wassersteinms

#endif // HEURISTIC_HPP__