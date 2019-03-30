/**
    idu.cpp
    IDU matrix factorization for mass spectrum deconvolution
    
    @author Antoine Passemiers
    @version 1.0 30/03/2019
*/

#include "ldu.hpp"

namespace wassersteinms {

void LDUDecomposition::factorize(Eigen::MatrixXd A, Eigen::VectorXd h) {
    Eigen::VectorXd h1 = h.head(this->n);
    Eigen::VectorXd h2 = h.segment(this->n, this->n);
    this->kinv = (1.0 / (h1 + h2).array()).matrix();
    this->kinv[this->n - 1] = 0;
    this->alpha = (h2 - h1).cwiseProduct(kinv);

    this->h3 = h.tail(this->k);

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
    //g[this->n - 1] = (h1 + h2)[this->n - 1];
    Eigen::MatrixXd G = g.asDiagonal();

    Eigen::MatrixXd GGG = g.tail(this->k).asDiagonal();
    //std::cout << GGG << std::endl;
    std::cout << "Rank of G: " << Eigen::ColPivHouseholderQR<Eigen::MatrixXd>(G).rank() << std::endl;

    Eigen::MatrixXd H3 = this->h3.asDiagonal();
    this->P = this->F * G * this->F.transpose() + H3;
    this->PDecomposition = Eigen::ColPivHouseholderQR<Eigen::MatrixXd>(this->P);

    this->Sigma = A * h.asDiagonal() * A.transpose(); // FIXME
}


Eigen::VectorXd LDUDecomposition::solve(Eigen::VectorXd r) {
    Eigen::VectorXd r1 = r.head(this->n - 1);
    Eigen::VectorXd r2 = r.tail(this->k);

    // Solve L . v1 = r
    Eigen::VectorXd v1 = Eigen::VectorXd(this->n + this->k - 1);
    v1.head(this->n - 1) = r1;
    v1.tail(this->k) = r2 - this->F * this->L.transpose() * this->kinv.head(this->n - 1).cwiseProduct(r1);

    // Solve M . v2 = v1
    Eigen::VectorXd v2 = Eigen::VectorXd(this->n + this->k - 1);
    v2.head(this->n - 1) = this->kinv.head(this->n - 1).cwiseProduct(v1.head(this->n - 1));
    v2.tail(this->k) = this->PDecomposition.solve(v1.tail(this->k));

    //std::cout << this->M2 * v2.tail(this->k) - v1.tail(this->k) << std::endl;
    // std::cout << this->M2 << std::endl;

    // Solve R . v3 = v2
    Eigen::VectorXd v3 = Eigen::VectorXd(this->n + this->k - 1);
    v3.head(this->n - 1) = v2.head(this->n - 1) - this->kinv.head(
            this->n - 1).cwiseProduct(this->L * this->F.transpose() * v2.tail(this->k));
    v3.tail(this->k) = v2.tail(this->k);

    //std::cout << (this->Sigma * v3 - r).tail(this->k) << std::endl; // FIXME

    return v3;
}


double diagonalMatrixLogConditionNumber(Eigen::MatrixXd &M) {
    Eigen::VectorXd absEigenvalues = M.diagonal().cwiseAbs();
    return std::log(absEigenvalues.maxCoeff()) - std::log(absEigenvalues.minCoeff());
}

} // namespace wassersteinms