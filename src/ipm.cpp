#include "ipm.hpp"


std::unique_ptr<ProblemInstance> formulateProblem(
        std::vector<std::unique_ptr<Spectrum>> &mu,
        std::unique_ptr<Spectrum> &nu) {

    size_t k = mu.size(), n = nu->length();
    ProblemInstance *prob = new ProblemInstance(n, k);
    Eigen::MatrixXd &F = prob->F;
    Eigen::VectorXd &b = prob->b;
    Eigen::VectorXd &c = prob->c;
    for (size_t j = 0; j < n-1; j++) {
        b[j] = nu->getRatio(j) - nu->getRatio(j + 1);
    }
    for (size_t j = 0; j < n; j++) {
        intensity_t g = nu->getIntensity(j);
        c[j] = g;
        c[j+n] = -g;
        for (size_t i = 0; i < k; i++) {
            F(i, j) = mu[i]->getIntensity(j);
        }
    }
    return std::unique_ptr<ProblemInstance>(prob);
}


std::unique_ptr<IpmSolution> createInitialSolution(std::unique_ptr<ProblemInstance> &prob) {
    std::unique_ptr<IpmSolution> sol(new IpmSolution(prob->n, prob->k));
    // TODO: initialize x
    // TODO: initialize y
    // TODO: initialize z
    return sol;
}


IpmSolution interiorPointMethod(std::unique_ptr<ProblemInstance> &prob, double epsilon) {

    // Find initial solution
    std::unique_ptr<IpmSolution> initialSolution = createInitialSolution(prob);
    IpmSolution sol = *initialSolution;

    // Extract data from problem instance
    size_t n = prob->n, k = prob->k;
    Eigen::VectorXd b = prob->b, c = prob->c;
    Eigen::MatrixXd F = prob->F; // TODO: sparse matrix?

    // Initialize matrix A
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n+k-1, 2*n+k);
    Eigen::MatrixXd J = Eigen::MatrixXd::Identity(n, n).block(0, 0, n-1, n);
    A.block(0, 0, n-1, n) = -J;
    A.block(0, n, n-1, n) = -J;
    A.block(n-1, 0, k, n) = F;
    A.block(n-1, n, k, n) = -F;
    A.block(n-1, 2*n, k, k) = -Eigen::MatrixXd::Identity(k, k);

    for (int t = 0; t < 5; t++) { // TODO
        Eigen::VectorXd x = sol.x, y = sol.y, z = sol.z;

        // Compute duality gap
        double mu = x.dot(z) / static_cast<double>(2*n + k);

        // Compute feasibility gaps (primal residual and dual residual)
        // and check for epsilon-feasibility
        Eigen::VectorXd rp = b - A*x;
        Eigen::VectorXd rd = c - A.transpose()*y - z;
        if ((rp.norm() < epsilon) && (rd.norm() < epsilon)) break;

        // Choose sigma
        double sigma = 1.0; // TODO

        // Solve linear system
        Eigen::MatrixXd X = x.asDiagonal();
        Eigen::MatrixXd Zinv = z.asDiagonal().inverse();
        Eigen::MatrixXd Sigma = A * Zinv * X * A.transpose();
        Eigen::VectorXd r = b + A * Zinv * (X * rd - sigma * mu * Eigen::VectorXd::Ones(rd.size()));
        Eigen::VectorXd dy = Sigma.colPivHouseholderQr().solve(r);
        Eigen::VectorXd dz = rd - A.transpose() * dy;
        Eigen::VectorXd dx = -x + Zinv * (sigma * mu * Eigen::VectorXd::Ones(dz.size()) - X * dz);

        // Choose step lengths
        double alphaP = 0.5; // TODO: choose alphaP such that new x is positive
        double alphaD = 0.5; // TODO: choose alphaD such that new z is positive

        // Update solution
        x += alphaP * dx;
        z += alphaD * dz;
        y += alphaD * dy;
    }

    return sol; // TODO
}