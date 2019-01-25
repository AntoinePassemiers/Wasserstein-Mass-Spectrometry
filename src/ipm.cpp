#include "ipm.hpp"


std::unique_ptr<ProblemInstance> formulateProblem(
        std::vector<std::unique_ptr<Spectrum>> &mu,
        std::unique_ptr<Spectrum> &nu) {

    size_t k = mu.size(), n = nu->length();
    ProblemInstance *prob = new ProblemInstance(n, k);
    Eigen::MatrixXd &F = prob->F;
    Eigen::VectorXd &b = prob->b;
    Eigen::VectorXd &c = prob->c;
    Eigen::MatrixXd &A = prob->A;
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
    F.col(n-1) = Eigen::VectorXd::Ones(k); // TODO?

    Eigen::MatrixXd J = Eigen::MatrixXd::Identity(n, n).block(0, 0, n-1, n);
    A.block(0, 0, n-1, n) = -J;
    A.block(0, n, n-1, n) = -J;
    A.block(n-1, 0, k, n) = F;
    A.block(n-1, n, k, n) = -F;
    A.block(n-1, 2*n, k, k) = -Eigen::MatrixXd::Identity(k, k);
    return std::unique_ptr<ProblemInstance>(prob);
}


std::unique_ptr<IpmSolution> createInitialSolution(std::unique_ptr<ProblemInstance> &prob) {
    std::unique_ptr<IpmSolution> sol(new IpmSolution(prob->n, prob->k));
    size_t k = prob->k, n = prob->n;

    // Initialize x
    for (size_t i = 0; i < n-1; i++) sol->x[i] = sol->x[i+n] = prob->b[i] / 2.0;
    sol->x[n] = 2.0 / 3.0;
    sol->x[2*n] = 1.0 / 3.0;
    for (size_t i = 0; i < k; i++) sol->x[2*n+i] = 1.0 / 3.0;

    // Initialize y
    Eigen::VectorXd p = Eigen::VectorXd::Constant(k, 1.0 / static_cast<float>(k));
    Eigen::VectorXd t = Eigen::VectorXd::Ones(n-1); // TODO
    sol->y.head(n-1) = t;
    sol->y.tail(k) = p;

    // Initialize z
    sol->z = prob->c - prob->A.transpose() * sol->y;

    return sol;
}


std::unique_ptr<IpmSolution> interiorPointMethod(std::unique_ptr<ProblemInstance> &prob, double epsilon) {

    // Find initial solution
    std::unique_ptr<IpmSolution> sol = createInitialSolution(prob);

    // Extract data from problem instance and initial solution
    size_t n = prob->n, k = prob->k;
    Eigen::VectorXd b = prob->b, c = prob->c;
    Eigen::MatrixXd A = prob->A;
    Eigen::VectorXd &x = sol->x, &y = sol->y, &z = sol->z;

    for (int t = 0; t < 100; t++) { // TODO

        // Compute centrality (related to duality gap)
        double mu = x.dot(z) / static_cast<double>(2*n + k);

        // Compute feasibility gaps (primal residual and dual residual)
        // and check for epsilon-feasibility
        Eigen::VectorXd rp = b - A*x;
        Eigen::VectorXd rd = c - A.transpose()*y - z;
        std::cout << "Iteration: " << t+1 << " - Centrality: " << mu << " - Residuals: ";
        std::cout << rp.norm() << ", " << rd.norm() << std::endl;
        if ((rp.norm() < epsilon) and (rd.norm() < epsilon) and (std::abs(mu) < epsilon)) break;

        // Choose sigma
        double sigma = 0.1; // TODO

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

    return sol;
}