#ifndef IPM_HPP__
#define IPM_HPP__

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

#include "spectrum.hpp"
#include "ldu.hpp"


typedef struct _ipmSolution {
    Eigen::VectorXd x;
    Eigen::VectorXd y;
    Eigen::VectorXd z;
    _ipmSolution() = default;
    _ipmSolution(size_t n, size_t k): x(2*n+k), y(n+k-1), z(2*n+k) {}
} IpmSolution;


typedef struct _problemInstance {
    size_t n; 
    size_t k;
    Eigen::MatrixXd F;
    Eigen::VectorXd b;
    Eigen::VectorXd c;
    Eigen::MatrixXd A;
    _problemInstance(size_t n, size_t k): n(n), k(k) {
        F = Eigen::MatrixXd::Zero(k, n);
        b = Eigen::VectorXd::Zero(n+k-1);
        c = Eigen::VectorXd::Zero(2*n+k);
        A = Eigen::MatrixXd::Zero(n+k-1, 2*n+k);
    }
} ProblemInstance;


double findPositivityConstrainedStepLength(
        Eigen::VectorXd &x,
        Eigen::VectorXd &dx,
        double alpha0);


bool isPrimalFeasible(
        std::unique_ptr<IpmSolution> &sol,
        std::unique_ptr<ProblemInstance> &prob,
        double epsilon);


bool isDualFeasible(
        std::unique_ptr<IpmSolution> &sol,
        std::unique_ptr<ProblemInstance> &prob,
        double epsilon);


bool isFeasible(
        std::unique_ptr<IpmSolution> &sol,
        std::unique_ptr<ProblemInstance> &prob,
        double epsilon);


bool satisfiesKKTConditions(
        std::unique_ptr<IpmSolution> &sol,
        std::unique_ptr<ProblemInstance> &prob,
        double epsilon);


std::unique_ptr<ProblemInstance> formulateProblem(
        std::vector<Spectrum> &mu,
        Spectrum &nu);


std::unique_ptr<IpmSolution> createInitialSolution(
        std::unique_ptr<ProblemInstance> &prob);


std::unique_ptr<IpmSolution> mehrotraPredictorCorrectorMethod(
        std::unique_ptr<ProblemInstance> &prob, double epsilon, size_t nMaxIterations);


std::unique_ptr<IpmSolution> longStepPathFollowingMethod(
        std::unique_ptr<ProblemInstance> &prob, double epsilon, size_t nMaxIterations);


#endif // IPM_HPP__
