/**
    solver.hpp
    Linear solver with preconditioning
    
    @author Antoine Passemiers
    @version 1.0 26/05/2019
*/

#ifndef SOLVER_HPP__
#define SOLVER_HPP__

#include <iostream>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Dense>

namespace wassersteinms {

class Solver {
private:
    // Size of the system
    size_t n;

    // Matrix representing the linear system
    Eigen::MatrixXd A;

    // Whether to precondition the system
    bool precondition;

    // Preconditioner matrix
    Eigen::MatrixXd P;

    // Inverse preconditioner matrix
    Eigen::MatrixXd Pinv;

    // Rank of the matrix
    size_t rank;

    // Matrix decomposition
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> ADecomposition;

    // Solution to the linear system
    Eigen::VectorXd v;
public:
    Solver(): precondition(false) {};
    Solver(bool precondition): precondition(precondition) {};
    void prepare(const Eigen::MatrixXd &A);
    const Eigen::MatrixXd& computePreconditioner(const Eigen::MatrixXd &A);
    const Eigen::VectorXd& solve(const Eigen::VectorXd &r);
};


} // namespace wassersteinms

#endif // SOLVER_HPP__