#include "ipm.hpp"


ipmSolution interiorPointMethod(ipmSolution sol0, problemInstance prob, double epsilon) {

    // Extract data from problem instance
    size_t n = prob.n, k = prob.k;
    Eigen::VectorXd b = prob.b, c = prob.c;
    Eigen::MatrixXd A = prob.A;

    // Initialize linear system
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(5*n+3*k-1, 5*n+3*k-1);
    M.block(0, 0, n+k-1, 2*n+k) = A;
    M.block(n+k-1, 2*n+k, 3*n+2*k-1, 3*n+2*k-1) = A.transpose();
    M.block(n+k-1, 3*n+2*k-1, 3*n+2*k-1, 5*n+3*k-1) = Eigen::MatrixXd::Identity(2*n+k, 2*n+k);
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(5*n+3*k-1); // TODO

    ipmSolution sol = sol0;
    for (int t = 0; t < 5; t++) { // TODO
        Eigen::VectorXd x = sol.x, y = sol.y, z = sol.z;
        double mu = x.dot(z) / static_cast<double>(2*n + k);
        Eigen::VectorXd rp = b - A*x;
        Eigen::VectorXd rd = c - A.transpose()*y - z;

        // Choose sigma
        double sigma = 1.0; // TODO

        // Update matrix M an rhs vector
        M.block(3*n+2*k-1, 0, 5*n+3*k-1, 2*n+k) = z.asDiagonal();
        M.block(3*n+2*k-1, 3*n+2*k-1, 5*n+3*k-1, 5*n+3*k-1) = x.asDiagonal();
        rhs.head(n+k-1) = rp;
        rhs.segment(n+k-1, 3*n+2*k-1) = rd;
        rhs.tail(2*n+k) = sigma * mu * Eigen::VectorXd::Ones(2*n+k) - x.asDiagonal() * z;

        // Solve linear system
        Eigen::VectorXd d = M.colPivHouseholderQr().solve(rhs);
        Eigen::VectorXd dx = d.head(n+k-1);
        Eigen::VectorXd dy = d.segment(n+k-1, 3*n+2*k-1);
        Eigen::VectorXd dz = d.tail(2*n+k);

        // Choose step lengths
        double alphaP = 0.5; // TODO
        double alphaD = 0.5; // TODO

        // Update solution
        x += alphaP * dx;
        z += alphaD * dz;
        y += alphaD * dy;
    }

    return sol; // TODO
}