#include "ipm.hpp"


double findPositivityConstrainedStepLength(Eigen::VectorXd &x, Eigen::VectorXd &dx, double alpha0) {
    // x_i + alpha * dx_i >= 0
    //     alpha >= -x_i / dx_i     if dx_i > 0 and x_i >= 0 -> ok since alpha >= 0
    //     alpha <  -x_i / dx_i     if dx_i < 0 and x_i >= 0
    double alpha = alpha0;
    for (auto i = 0; i < x.size(); i++) {
        double ratio = -x[i] / dx[i];
        if (dx[i] < 0.0) {
            if (ratio < alpha) alpha = ratio;
        }
    }
    if (std::isnan(alpha)) alpha = 0.0;
    return alpha;
}


bool isFeasible(
        std::unique_ptr<IpmSolution> &sol,
        std::unique_ptr<ProblemInstance> &prob,
        double epsilon) {
    auto rp = ((prob->A * sol->x) - prob->b).array().abs();
    auto rd = (prob->A.transpose() * sol->y + sol->z - prob->c).array().abs();
    return (sol->x.array() >= -epsilon).all()
        && (sol->z.array() >= -epsilon).all()
        && ((rp < epsilon).all())
        && ((rd < epsilon).all());
}


bool satisfiesKKTConditions(
        std::unique_ptr<IpmSolution> &sol,
        std::unique_ptr<ProblemInstance> &prob,
        double epsilon) {
    return isFeasible(sol, prob, epsilon) && (sol->x.dot(sol->z) <= epsilon);
}


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
        mz_t dj = nu->getRatio(j+1) - nu->getRatio(j);
        b[j] = -dj;
    }
    intensity_t gj = 0.0;
    for (size_t j = 0; j < n; j++) {
        gj += nu->getIntensity(j);
        c[j] = gj;
        c[j+n] = -gj;
    }

    for (size_t i = 0; i < k; i++) {
    	intensity_t fj = 0.0;
    	for (size_t j = 0; j < n; j++) {
            fj += mu[i]->getIntensity(j);
            F(i, j) = fj;
        }
    }

    Eigen::MatrixXd J = Eigen::MatrixXd::Identity(n, n).block(0, 0, n-1, n);
    A.block(0, 0, n-1, n) = -J;
    A.block(0, n, n-1, n) = -J;
    A.block(n-1, 0, k, n) = F;
    A.block(n-1, n, k, n) = -F;
    A.block(n-1, 2*n, k, k) = -Eigen::MatrixXd::Identity(k, k);

    // Check positive-definiteness of A
    Eigen::FullPivLU<Eigen::MatrixXd> luDecomposition(A);
    // TODO

    return std::unique_ptr<ProblemInstance>(prob);
}


std::unique_ptr<IpmSolution> createInitialSolution(std::unique_ptr<ProblemInstance> &prob) {
    std::unique_ptr<IpmSolution> sol(new IpmSolution(prob->n, prob->k));
    size_t k = prob->k, n = prob->n;

    // Initialize x
    for (size_t i = 0; i < n-1; i++) sol->x[i] = -prob->b[i] / 2.0;
    sol->x[n-1] = 2.0 / 3.0;
    for (size_t i = n; i < 2*n-1; i++) sol->x[i] = -prob->b[i-n] / 2.0;
    sol->x[2*n-1] = 1.0 / 3.0;
    for (size_t i = 2*n; i < 2*n+k; i++) sol->x[i] = 1.0 / 3.0;

    // Initialize y
    Eigen::VectorXd p = Eigen::VectorXd::Constant(k, 1.0 / static_cast<float>(k));
    Eigen::VectorXd t = Eigen::VectorXd::Ones(n-1);
    for (size_t j = 0; j < n-1; j++) {
        t[j] += std::abs(prob->F.col(j).dot(p) - prob->c[j]);
    }
    sol->y.head(n-1) = t;
    sol->y.tail(k) = p;

    // Initialize z
    sol->z = prob->c - prob->A.transpose() * sol->y;

    return sol;
}


std::unique_ptr<IpmSolution> interiorPointMethod(std::unique_ptr<ProblemInstance> &prob, double epsilon, size_t nMaxIterations) {

    std::function<Eigen::VectorXd (Eigen::VectorXd)> nantonumVec = [](Eigen::VectorXd x) {
        return (x.array() != x.array()).select(0.0, x);
    };

    std::function<Eigen::MatrixXd (Eigen::MatrixXd)> nantonumMat = [](Eigen::MatrixXd x) {
        return (x.array() != x.array()).select(0.0, x);
    };

    // Find initial solution
    std::unique_ptr<IpmSolution> sol = createInitialSolution(prob);
    assert(isFeasible(sol, prob, 1e-5)); // TODO

    // Extract data from problem instance and initial solution
    Eigen::VectorXd b = prob->b, c = prob->c;
    Eigen::MatrixXd A = prob->A;
    Eigen::VectorXd &x = sol->x, &y = sol->y, &z = sol->z;

    for (size_t t = 0; t < nMaxIterations; t++) {

        // assert(isFeasible(sol, prob, 1e-5)); // TODO

        Eigen::MatrixXd X = x.asDiagonal();
        Eigen::MatrixXd Z = z.asDiagonal();
        Eigen::MatrixXd Zinv = (1.0 / z.array()).matrix().asDiagonal().inverse();

        // -- PREDICTOR STEP (AFFINE SCALING) --
	
        // Compute feasibility gaps (primal residual and dual residual)
        Eigen::VectorXd rp = b - A * x;
        Eigen::VectorXd rd = c - A.transpose() * y - z;
    	rp = nantonumVec(rp); rd = nantonumVec(rd);

        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> SigmaDecomposition(A * (Zinv * X) * A.transpose());
    	Eigen::VectorXd r = b + A * Zinv * X * rd;
    	Eigen::VectorXd dyAff = SigmaDecomposition.solve(r);
    	dyAff = nantonumVec(dyAff);
    	Eigen::VectorXd dzAff = rd - A.transpose() * dyAff;
    	dzAff = nantonumVec(dzAff);
    	Eigen::VectorXd dxAff = -x - Zinv * X * dzAff;
    	dxAff = nantonumVec(dxAff);

        // Compute centrality (duality measure)
        double mu = x.dot(z) / static_cast<double>(x.size());
        if (std::isnan(mu)) mu = 0.0;

        // Choose step lengths while satisfying positivity constraints:
        // Choose alphaP such that new x is positive
        // Choose alphaD such that new z is positive
        double alphaP = findPositivityConstrainedStepLength(x, dxAff, 1.0);
        double alphaD = findPositivityConstrainedStepLength(z, dzAff, 1.0);

        // Compute centering parameter
        // double muAff = (x + alphaP*dxAff).dot(z + alphaP*dzAff) / static_cast<double>(x.size());
    	double muAff = (x + alphaP*dxAff).dot(z + alphaD*dzAff) / static_cast<double>(x.size());
        double sigma = std::pow(muAff / mu, 3.0);
    	if (std::isnan(sigma)) sigma = 0.0;

    	// -- CENTER-CORRECTOR STEP (AGGREGATED SYSTEM)

    	// Check the epsilon-feasibility of current solution
        std::cout << "Iteration: " << t+1 << " - Centrality: " << mu << " - Residuals: ";
    	std::cout << rp.norm() << ", " << rd.norm() << std::endl;
        std::cout << "\tObjectives: " << y.dot(b) << " <= " << x.dot(c) << std::endl;
        if ((rp.norm() < epsilon) and (rd.norm() < epsilon) and (std::abs(mu) < epsilon)) break;

        Eigen::VectorXd correction = dxAff.asDiagonal() * dzAff;

    	r = b + A * Zinv * (X * rd - sigma * mu * Eigen::VectorXd::Ones(rd.size()) + correction);
    	Eigen::VectorXd dy = SigmaDecomposition.solve(r);
    	dy = nantonumVec(dy);
    	Eigen::VectorXd dz = rd - A.transpose() * dy;
    	dz = nantonumVec(dz);
    	Eigen::VectorXd dx = -x + Zinv * (-X * dz + sigma * mu * Eigen::VectorXd::Ones(dz.size()) - correction);
    	dx = nantonumVec(dx);

    	double eta = 0.9; // TODO: increase eta at each iteration -> 1.0
        alphaP = findPositivityConstrainedStepLength(x, dx, 1.0);
        alphaD = findPositivityConstrainedStepLength(z, dz, 1.0);
        alphaP = std::min(1.0, eta * alphaP);
        alphaD = std::min(1.0, eta * alphaD);

        std::cout << "\talphaP: " << alphaP << " - alphaD: " << alphaD << " - sigma: ";
        std::cout << sigma << "\n" << std::endl;

        // Update solution
        x += alphaP * dx;
        x = nantonumVec((x.array() < 0).select(1e-20, x));
        z += alphaD * dz;
        z = nantonumVec((z.array() < 0).select(1e-20, z));
        y += alphaD * dy;
        y = nantonumVec(y);
    }

    return sol;
}
