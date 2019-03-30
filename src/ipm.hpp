/**
    ipm.hpp
    Interior-point methods for mass spectrum deconvolution
    
    @author Antoine Passemiers
    @version 1.0 30/03/2019
*/

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

namespace wassersteinms {

typedef struct _ipmSolution {

    // Solution to the primal problem
    Eigen::VectorXd x;

    // Solution to the dual problem.
    // The last k elements of this vector are the solution
    // to the spectral deconvolution problem itself.
    Eigen::VectorXd y;

    // Slack variables of the dual problem
    Eigen::VectorXd z;

    _ipmSolution() = default;
    _ipmSolution(size_t n, size_t k): x(2*n+k), y(n+k-1), z(2*n+k) {}
} IpmSolution;


typedef struct _problemInstance {

    // Size of theoretical and empirical spectra
    size_t n;

    // Number of theoretical spectra
    size_t k;

    // Matrix of shape (k x n) where row i is the cumulative
    // probability distribution of theoretical spectrum i.
    Eigen::MatrixXd F;

    // Vector of size n+k-1, where n-1 first elements are
    // negative bin sizes, and k last elements are zeroes
    Eigen::VectorXd b;

    // Vector of size 2n+k, where n first elements are the
    // cumulative probability distribution (cdf) of the empirical
    // spectrum, n next elements are the negative cdf,
    // and k last elements are zeroes
    Eigen::VectorXd c;

    // Matrix of shape (n+k-1 x 2n+k) described as a block matrix:
    // [ [ -J, -J, 0 ], [ F, -F, -I ] ], where I is an identity matrix
    // of shape (k x k) and J is an identity matrix of shape (n-1 x n)
    // without the last row
    Eigen::MatrixXd A;
    
    _problemInstance(size_t n, size_t k): n(n), k(k) {
        F = Eigen::MatrixXd::Zero(k, n);
        b = Eigen::VectorXd::Zero(n+k-1);
        c = Eigen::VectorXd::Zero(2*n+k);
        A = Eigen::MatrixXd::Zero(n+k-1, 2*n+k);
    }
} ProblemInstance;

/**
    Computes the largest step length that can be used to update
    current solution without violating non-negativity constraints

    @param x Current solution (either primal or dual)
    @param dx Search direction
    @param alpha0 Maximum value for the step length
    @return Maximum step length that satisfies
        non-negativity constraints
*/
double findPositivityConstrainedStepLength(
        Eigen::VectorXd &x,
        Eigen::VectorXd &dx,
        double alpha0);

/**
    Checks whether primal solution is epsilon-feasible

    @param sol Current solution
    @param prob An instance of spectral deconvolution problem
    @param epsilon Feasibility threshold
    @return Whether solution is epsilon-feasible
*/
bool isPrimalFeasible(
        std::unique_ptr<IpmSolution> &sol,
        std::unique_ptr<ProblemInstance> &prob,
        double epsilon);

/**
    Checks whether dual solution is epsilon-feasible

    @param sol Current solution
    @param prob An instance of spectral deconvolution problem
    @param epsilon Feasibility threshold
    @return Whether solution is epsilon-feasible
*/
bool isDualFeasible(
        std::unique_ptr<IpmSolution> &sol,
        std::unique_ptr<ProblemInstance> &prob,
        double epsilon);

/**
    Checks whether solution to the primal-dual problem
    is epsilon-feasible

    @param sol Current solution
    @param prob An instance of spectral deconvolution problem
    @param epsilon Feasibility threshold
    @return Whether solution is epsilon-feasible
*/
bool isFeasible(
        std::unique_ptr<IpmSolution> &sol,
        std::unique_ptr<ProblemInstance> &prob,
        double epsilon);

/**
    Checks whether solution to the primal-dual problem
    satisfies KKT solution

    @param sol Current solution
    @param prob An instance of spectral deconvolution problem
    @param epsilon Feasibility threshold
    @return Whether solution is epsilon-optimal
*/
bool satisfiesKKTConditions(
        std::unique_ptr<IpmSolution> &sol,
        std::unique_ptr<ProblemInstance> &prob,
        double epsilon);

/**
    Formulates a spectral deconvolution problem

    @param mu Vector of theoretical spectra
    @param nu Empirical spectrum
*/
std::unique_ptr<ProblemInstance> formulateProblem(
        std::vector<Spectrum> &mu,
        Spectrum &nu);

/**
    Creates an initial primal-dual feasible solution
    to spectral deconvolution problem

    @param prob An instance of the spectral deconvolution problem
    @return A primal-dual feasible solution
*/
std::unique_ptr<IpmSolution> createInitialSolution(
        std::unique_ptr<ProblemInstance> &prob);

/**
    Mehrotra predictor-corrector-method

    @param prob An instance of the spectral deconvolution problem
    @param epsilon Convergence threshold
    @param nMaxIterations Number maximum of iterations of the
        interior-point method
    @return An (eventually) epsilon-optimal solution
*/
std::unique_ptr<IpmSolution> mehrotraPredictorCorrectorMethod(
        std::unique_ptr<ProblemInstance> &prob,
        double epsilon,
        size_t nMaxIterations);

/**
    Long-step path-following method

    @param prob An instance of the spectral deconvolution problem
    @param epsilon Convergence threshold
    @param nMaxIterations Number maximum of iterations of the
        interior-point method
    @return An (eventually) epsilon-optimal solution
*/
std::unique_ptr<IpmSolution> longStepPathFollowingMethod(
        std::unique_ptr<ProblemInstance> &prob,
        double epsilon,
        size_t nMaxIterations);

} // namespace wassersteinms

#endif // IPM_HPP__
