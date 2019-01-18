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
    Eigen::MatrixXd F = prob->F;

    // Initialize matrix A
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n+k-1, 2*n+k);
    Eigen::MatrixXd J = Eigen::MatrixXd::Identity(n, n).block(0, 0, n-1, n);
    A.block(0, 0, n-1, n) = -J;
    A.block(0, n, n-1, n) = -J;
    A.block(n-1, 0, k, n) = F;
    A.block(n-1, n, k, n) = -F;
    A.block(n-1, 2*n, k, k) = -Eigen::MatrixXd::Identity(k, k);

    // Initialize linear system
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(5*n+3*k-1, 5*n+3*k-1);
    M.block(0, 0, n+k-1, 2*n+k) = A;
    M.block(n+k-1, 2*n+k, 2*n+k, n+k-1) = A.transpose();
    M.block(n+k-1, 3*n+2*k-1, 2*n+k, 2*n+k) = Eigen::MatrixXd::Identity(2*n+k, 2*n+k);
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(5*n+3*k-1); // TODO

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

        // Update matrix M and rhs vector
        M.block(3*n+2*k-1, 0, 2*n+k, 2*n+k) = z.asDiagonal();
        M.block(3*n+2*k-1, 3*n+2*k-1, 2*n+k, 2*n+k) = x.asDiagonal();
        std::cout << x.size() << ", " << rp.size() << ", " << rd.size() << std::endl;
        rhs.head(rp.size()) = rp;
        rhs.segment(rp.size(), rd.size()) = rd;
        rhs.tail(z.size()) = sigma * mu * Eigen::VectorXd::Ones(z.size()) - x.asDiagonal() * z;

        // Solve linear system
        Eigen::VectorXd d = M.colPivHouseholderQr().solve(rhs);
        Eigen::VectorXd dx = d.head(x.size());
        Eigen::VectorXd dy = d.segment(x.size(), y.size());
        Eigen::VectorXd dz = d.tail(z.size());

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