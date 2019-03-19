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

    Eigen::VectorXd d = Eigen::VectorXd::Zero(n-1);
    for (size_t j = 0; j < n-1; j++) {
        d[j] = nu->getRatio(j+1) - nu->getRatio(j);
    }
    b.head(n-1) = -d;
    b.tail(k) = Eigen::VectorXd::Zero(k);

    // Compute cdf of empirical spectrum
    Eigen::VectorXd g = Eigen::VectorXd::Zero(n);
    intensity_t gj = 0.0;
    for (size_t j = 0; j < n; j++) {
        gj += nu->getIntensity(j);
        g[j] = gj;
    }
    c.head(n) = g;
    c.segment(n, n) = -g;
    c.tail(k) = Eigen::VectorXd::Zero(k);

    // Compute cdfs of theoretical spectra
    for (size_t i = 0; i < k; i++) {
    	intensity_t fj = 0.0;
    	for (size_t j = 0; j < n; j++) {
            fj += mu[i]->getIntensity(j);
            F(i, j) = fj;
        }
    }

    // Compute constraint matrix A
    Eigen::MatrixXd J = Eigen::MatrixXd::Identity(n, n).block(0, 0, n-1, n);
    A.block(0, 0, n-1, n) = -J;
    A.block(0, n, n-1, n) = -J;
    A.block(n-1, 0, k, n) = F;
    A.block(n-1, n, k, n) = -F;
    A.block(n-1, 2*n, k, k) = -Eigen::MatrixXd::Identity(k, k);

    // Check rank of A
    Eigen::FullPivLU<Eigen::MatrixXd> luDecomposition(A);
    int rank = luDecomposition.rank();
    if (rank != A.rows()) {
        std::cout << "[Warning] A is of shape (" << A.rows() << ", " << A.cols() << ")";
        std::cout << " and rank " << rank << std::endl;
        std::cout << "[Warning] A is not a full-rank matrix" << std::endl;
    }

    return std::unique_ptr<ProblemInstance>(prob);
}


std::unique_ptr<IpmSolution> createInitialSolution(std::unique_ptr<ProblemInstance> &prob) {
    std::unique_ptr<IpmSolution> sol(new IpmSolution(prob->n, prob->k));
    size_t k = prob->k, n = prob->n;

    // Initialize x
    Eigen::ArrayXd d = -prob->b.head(n-1).array();
    sol->x.segment(0, n-1) = (d / 2.0).matrix();
    sol->x[n-1] = 2.0 / 3.0;
    sol->x.segment(n, n-1) = (d / 2.0).matrix();
    sol->x[2*n-1] = 1.0 / 3.0;
    sol->x.tail(k) = Eigen::VectorXd::Ones(k) * (1.0 / 3.0);

    // Initialize y
    Eigen::VectorXd g = prob->c.head(n);
    Eigen::VectorXd p = Eigen::VectorXd::Constant(k, 1.0 / static_cast<float>(k));
    Eigen::VectorXd t = ((prob->F.transpose() * p - g).head(n-1).cwiseAbs().array() + 1.0).matrix();
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

    // Find initial solution
    std::unique_ptr<IpmSolution> sol = createInitialSolution(prob);
    assert(isFeasible(sol, prob, 1e-5)); // TODO

    LDUDecomposition sigmaDecomposition(prob->n, prob->k);

    // Extract data from problem instance and initial solution
    Eigen::VectorXd b = prob->b, c = prob->c;
    Eigen::MatrixXd A = prob->A;
    Eigen::VectorXd &x = sol->x, &y = sol->y, &z = sol->z;

    for (size_t t = 0; t < nMaxIterations; t++) {

        // assert(isFeasible(sol, prob, 1e-5)); // TODO

        Eigen::VectorXd zinv = (1.0 / z.array()).matrix();
        Eigen::VectorXd zinvx = x.cwiseQuotient(z);
        Eigen::VectorXd _1 = Eigen::VectorXd::Ones(z.size());

        // -- PREDICTOR STEP (AFFINE SCALING) --
	
        // Compute feasibility gaps (primal residual and dual residual)
        Eigen::VectorXd rp = b - A * x;
        Eigen::VectorXd rd = c - A.transpose() * y - z;

        // Compute centrality (duality measure)
        double mu = x.dot(z) / static_cast<double>(x.size());

        // Compute affine scaling direction
        sigmaDecomposition.factorize(A, zinvx);
    	Eigen::VectorXd r = b + A * zinvx.cwiseProduct(rd);
    	Eigen::VectorXd dyAff = sigmaDecomposition.solve(r);
    	//dyAff = nantonumVec(dyAff);
    	Eigen::VectorXd dzAff = rd - A.transpose() * dyAff;
    	//dzAff = nantonumVec(dzAff);
    	Eigen::VectorXd dxAff = -x - zinvx.cwiseProduct(dzAff);
    	//dxAff = nantonumVec(dxAff);

        // Choose step lengths while satisfying positivity constraints:
        // Choose alphaP such that new x is positive
        // Choose alphaD such that new z is positive
        double alphaP = findPositivityConstrainedStepLength(x, dxAff, 1.0);
        double alphaD = findPositivityConstrainedStepLength(z, dzAff, 1.0);

        // Compute centering parameter
        // double muAff = (x + alphaP * dxAff).dot(z + alphaP * dzAff) / static_cast<double>(x.size());
    	double muAff = (x + alphaP * dxAff).dot(z + alphaD * dzAff) / static_cast<double>(x.size());
        double sigma = std::pow(muAff / mu, 3.0);
    	if (std::isnan(sigma)) sigma = 0.0;
        sigma = std::max(0.0, sigma);

    	// -- CENTER-CORRECTOR STEP (AGGREGATED SYSTEM)

    	// Check the epsilon-feasibility of current solution
        std::cout << "Iteration: " << t+1 << " - Centrality: " << mu << " - Residuals: ";
    	std::cout << rp.norm() << ", " << rd.norm() << std::endl;
        std::cout << "\tObjectives: " << y.dot(b) << " <= " << x.dot(c) << std::endl;
        if ((rp.norm() < epsilon) and (rd.norm() < epsilon) and (std::abs(mu) < epsilon)) break;

        Eigen::VectorXd correction = dxAff.cwiseProduct(dzAff);

    	r = b + A * (zinvx.cwiseProduct(rd) - zinv.cwiseProduct(sigma * mu * _1 - correction));
    	Eigen::VectorXd dy = sigmaDecomposition.solve(r);
    	Eigen::VectorXd dz = rd - A.transpose() * dy;
    	Eigen::VectorXd dx = -x - zinvx.cwiseProduct(dz) + zinv.cwiseProduct(sigma * mu * _1 - correction);

        alphaP = findPositivityConstrainedStepLength(x, dx, 1.0);
        alphaD = findPositivityConstrainedStepLength(z, dz, 1.0);
        double eta = std::max(0.995, 1.0 - mu);
        alphaP = std::min(1.0, eta * alphaP);
        alphaD = std::min(1.0, eta * alphaD);

        std::cout << "\talphaP: " << alphaP << " - alphaD: " << alphaD << " - sigma: ";
        std::cout << sigma << "\n" << std::endl;

        // Update solution
        x += alphaP * dx;
        x = nantonumVec((x.array() < 1e-20).select(1e-20, x));
        z += alphaD * dz;
        z = nantonumVec((z.array() < 1e-20).select(1e-20, z));
        y += alphaD * dy;
        //y = nantonumVec(y);
    }

    return sol;
}
