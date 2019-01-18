#ifndef IPM_HPP__
#define IPM_HPP__

#include <memory>
#include <Eigen/Core>
#include <Eigen/Dense>

#include "spectrum.hpp"


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
    _problemInstance(size_t n, size_t k): n(n), k(k) {
        F = Eigen::MatrixXd::Zero(k, n);
        b = Eigen::VectorXd::Zero(n + k - 1);
        c = Eigen::VectorXd::Zero(2*n + k);
    }
} ProblemInstance;


std::unique_ptr<ProblemInstance> formulateProblem(
        std::vector<std::unique_ptr<Spectrum>> &mu,
        std::unique_ptr<Spectrum> &nu);


std::unique_ptr<IpmSolution> createInitialSolution(
        std::unique_ptr<ProblemInstance> &prob);


IpmSolution interiorPointMethod(
        std::unique_ptr<ProblemInstance> &prob, double epsilon);


#endif // IPM_HPP__