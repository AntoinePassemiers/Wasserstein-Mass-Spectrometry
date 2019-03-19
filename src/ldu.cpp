#include "ldu.hpp"


void LDUDecomposition::factorize(Eigen::MatrixXd A, Eigen::VectorXd h) {
    /**
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
    Eigen::MatrixXd G = Hplus - L.transpose() * this->Kinv * L;

    Eigen::MatrixXd H3 = this->h3.asDiagonal();
    Eigen::MatrixXd M2 = this->F * G * this->F.transpose() + H3;
    this->M2Decomposition = Eigen::ColPivHouseholderQR<Eigen::MatrixXd>(M2);
    */
    this->SigmaDecomposition = Eigen::ColPivHouseholderQR<Eigen::MatrixXd>(A * h.asDiagonal() * A.transpose());
}


Eigen::VectorXd LDUDecomposition::solve(Eigen::VectorXd r) {
    /**
    Eigen::VectorXd r1 = r.head(this->n - 1);
    Eigen::VectorXd r2 = r.tail(this->k);

    // Solve L . v1 = r
    Eigen::VectorXd v1 = Eigen::VectorXd(this->n + this->k - 1);
    v1.head(this->n - 1) = r1;
    v1.tail(this->k) = r2 - this->F * this->L.transpose() * this->Kinv * r1;

    // Solve M . v2 = v1
    Eigen::VectorXd v2 = Eigen::VectorXd(this->n + this->k - 1);
    v2.head(this->n - 1) = this->Kinv * v1.head(this->n - 1);
    v2.tail(this->k) = this->M2Decomposition.solve(v1.tail(this->k));

    // Solve R . v3 = v2
    Eigen::VectorXd v3 = Eigen::VectorXd(this->n + this->k - 1);
    v3.head(this->n - 1) = v2.head(this->n - 1) - this->Kinv * this->L * this->F.transpose() * v2.tail(this->k);
    v3.tail(this->k) = v2.tail(this->k);
    */
    Eigen::VectorXd v3 = SigmaDecomposition.solve(r);

    return v3;
}