#include "ldu.hpp"


void LDUDecomposition::factorize(Eigen::MatrixXd A, Eigen::VectorXd h) {
    
    Eigen::VectorXd h1 = h.head(this->n);
    Eigen::VectorXd h2 = h.segment(this->n, this->n);
    Eigen::VectorXd h3 = h.tail(this->k);

    this->F = A.block(this->n - 1, 0, this->k, this->n);

    // TODO: precompute sub-matrices
}


Eigen::VectorXd LDUDecomposition::solve(Eigen::VectorXd r) {

    Eigen::VectorXd r1 = r.head(this->n - 1);
    Eigen::VectorXd r2 = r.tail(this->k);

    Eigen::VectorXd v1 = Eigen::VectorXd(this->n + this->k - 1);
    v1.head(this->n - 1) = r1;
    v1.tail(this->k) = r2; // TODO

    // TODO

    Eigen::VectorXd v3 = Eigen::VectorXd(2 * this->n + this->k);

    // TODO

    return v3;
}