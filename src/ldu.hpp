/**
    idu.hpp
    IDU matrix factorization for mass spectrum deconvolution
    
    @author Antoine Passemiers
    @version 1.0 30/03/2019
*/

#ifndef LDU_HPP__
#define LDU_HPP__

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
    Eigen::VectorXd h3;
    Eigen::MatrixXd F;
    Eigen::MatrixXd L;
    Eigen::MatrixXd Kinv;
    Eigen::MatrixXd P;
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> PDecomposition;
    Eigen::MatrixXd Sigma; // FIXME
public:
    LDUDecomposition(size_t n, size_t k): n(n), k(k) {}
    void factorize(Eigen::MatrixXd A, Eigen::VectorXd h);
    Eigen::VectorXd solve(Eigen::VectorXd r);
};


double diagonalMatrixConditionNumber(Eigen::MatrixXd &M);

} // namespace wassersteinms

#endif // LDU_HPP__