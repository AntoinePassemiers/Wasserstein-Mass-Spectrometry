/**
    idu.hpp
    IDU matrix factorization for mass spectrum deconvolution
    
    @author Antoine Passemiers
    @version 1.0 30/03/2019
*/

#ifndef LDU_HPP__
#define LDU_HPP__

#include "solver/solver.hpp"

#include <iostream>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Dense>

namespace wassersteinms {

class LDUDecomposition {
private:
    // Size of spectra
    size_t n;

    // Number of theoretical spectra
    size_t k;

    Eigen::VectorXd alpha;
    Eigen::VectorXd kinv;
    Eigen::VectorXd h;
    Eigen::VectorXd h3;
    Eigen::MatrixXd A;
    Eigen::MatrixXd F;
    Eigen::MatrixXd L;
    Eigen::MatrixXd Kinv;
    Eigen::MatrixXd P;
    Solver solver;
public:
    LDUDecomposition(size_t n, size_t k): n(n), k(k), solver(false) {}
    void factorize(const Eigen::MatrixXd &A, const Eigen::VectorXd &h);
    Eigen::VectorXd solve(const Eigen::VectorXd &r);
};


bool checkSolutionCorrectness(const Eigen::MatrixXd &A,
                              const Eigen::VectorXd &h,
                              const Eigen::VectorXd &v,
                              const Eigen::VectorXd &r);

double logConditionNumber(const Eigen::VectorXd &eigenvalues);

double logConditionNumber(const Eigen::MatrixXd &M);

} // namespace wassersteinms

#endif // LDU_HPP__