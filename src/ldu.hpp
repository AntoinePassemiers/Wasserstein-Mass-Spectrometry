#ifndef LDU_HPP__
#define LDU_HPP__

#include <iostream>
#include <cmath>

#include <Eigen/Core>
#include <Eigen/Dense>


class LDUDecomposition {
private:
    size_t n;
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


#endif // LDU_HPP__