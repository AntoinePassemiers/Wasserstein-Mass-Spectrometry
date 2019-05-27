/**
    solver.cpp
    Linear solver with preconditioning
    
    @author Antoine Passemiers
    @version 1.0 26/05/2019
*/

#include "solver.hpp"

namespace wassersteinms {

void Solver::prepare(const Eigen::MatrixXd &A) {
    this->n = A.rows();
    this->A = A;
    Eigen::MatrixXd M;
    if (this->precondition) {
        Eigen::MatrixXd P = this->computePreconditioner(A);
        M = A * this->Pinv;
    } else {
        M = A;
    }
    this->ADecomposition = Eigen::ColPivHouseholderQR<Eigen::MatrixXd>(M);
    this->rank = this->ADecomposition.rank();
}

const Eigen::MatrixXd& Solver::computePreconditioner(const Eigen::MatrixXd &A) {
    // Jacobi preconditioner
    this->P = A.diagonal().asDiagonal();
    this->Pinv = P.inverse();
    return this->P;
}

const Eigen::VectorXd& Solver::solve(const Eigen::VectorXd &r) {
    if (this->rank > 1) {
        Eigen::VectorXd y = this->ADecomposition.solve(r);
        if (this->precondition) {
            this->v = this->Pinv * y;
        } else {
            this->v = y;
        }
    } else {
        // Naive approximation
        this->v = r / (this->n * this->A.mean());
    }
    return this->v;
}

} // namespace wassersteinms