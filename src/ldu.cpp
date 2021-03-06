/**
    idu.cpp
    IDU matrix factorization for mass spectrum deconvolution
    
    @author Antoine Passemiers
    @version 1.0 30/03/2019
*/

#include "ldu.hpp"

namespace wassersteinms {

void LDUDecomposition::factorize(const Eigen::MatrixXd &A,
                                 const Eigen::VectorXd &h) {
    Eigen::VectorXd h1 = h.head(this->n);
    Eigen::VectorXd h2 = h.segment(this->n, this->n);
    this->h = h;
    this->kinv = (1.0 / (h1 + h2).array()).matrix();
    this->kinv[this->n - 1] = 0;
    this->alpha = (h2 - h1).cwiseProduct(kinv);

    this->h3 = h.tail(this->k);
    this->A = A;
    this->F = A.block(this->n - 1, 0, this->k, this->n);

    Eigen::MatrixXd J = Eigen::MatrixXd::Identity(n, n).block(0, 0, n-1, n);
    Eigen::MatrixXd L = J * (h2 - h1).asDiagonal();
    Eigen::MatrixXd Hplus = (h1 + h2).asDiagonal();
    Eigen::MatrixXd K = J * Hplus * J.transpose();

    this->Kinv = K.inverse();
    this->L = L;

    // Eigen::MatrixXd G = Hplus - L.transpose() * this->Kinv * L;
    Eigen::VectorXd hmax = h1.cwiseMax(h2);
    Eigen::VectorXd hmin = h1.cwiseMin(h2);
    Eigen::VectorXd g = 4.0 * hmax.cwiseQuotient(h1 + h2).cwiseProduct(hmin);
    g[this->n - 1] = (h1 + h2)[this->n - 1];
    Eigen::MatrixXd G = g.asDiagonal();

    Eigen::MatrixXd H3 = this->h3.asDiagonal();
      this->P = this->F * G * this->F.transpose() + H3;
    this->solver.prepare(this->P);
}


Eigen::VectorXd LDUDecomposition::solve(const Eigen::VectorXd &r) {
    Eigen::VectorXd r1 = r.head(this->n - 1);
    Eigen::VectorXd r2 = r.tail(this->k);

    // Solve L . v1 = r
    Eigen::VectorXd v1 = Eigen::VectorXd(this->n + this->k - 1);
    v1.head(this->n - 1) = r1;
    v1.tail(this->k) = r2 - this->F * this->L.transpose() * this->kinv.head(this->n - 1).cwiseProduct(r1);

    // Solve M . v2 = v1
    Eigen::VectorXd v2 = Eigen::VectorXd(this->n + this->k - 1);
    v2.head(this->n - 1) = this->kinv.head(this->n - 1).cwiseProduct(v1.head(this->n - 1));
    // v2.head(this->n - 1) = v1.head(this->n - 1);
    v2.tail(this->k) = this->solver.solve(v1.tail(this->k));

    // Solve R . v3 = v2
    Eigen::VectorXd v3 = Eigen::VectorXd(this->n + this->k - 1);
    v3.head(this->n - 1) = v2.head(this->n - 1) - this->kinv.head(
            this->n - 1).cwiseProduct(this->L * this->F.transpose() * v2.tail(this->k));
    v3.tail(this->k) = v2.tail(this->k);

    std::cout << checkSolutionCorrectness(this->A, this->h, v3, r) << std::endl;
    return v3;
}


bool checkSolutionCorrectness(const Eigen::MatrixXd &A,
                              const Eigen::VectorXd &h,
                              const Eigen::VectorXd &v,
                              const Eigen::VectorXd &r) {
    Eigen::MatrixXd AHA = A * h.asDiagonal() * A.transpose();
    Eigen::VectorXd rd = (AHA * v - r);
    return (rd.norm() < 1e-5);
}


double logConditionNumber(const Eigen::VectorXd &eigenvalues) {
    Eigen::VectorXd absEigenvalues = eigenvalues.cwiseAbs();
    return std::log(absEigenvalues.maxCoeff()) - std::log(absEigenvalues.minCoeff());
}


double logConditionNumber(const Eigen::MatrixXd &M) {
    Eigen::EigenSolver<Eigen::MatrixXd> es(M);
    Eigen::VectorXd eigenvalues = es \
            .eigenvalues() \
            .asDiagonal() \
            .diagonal() \
            .real();
    return logConditionNumber(eigenvalues);
}

} // namespace wassersteinms