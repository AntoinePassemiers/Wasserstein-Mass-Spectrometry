/**
    heuristic.cpp
    Heuristic method for mass spectrum deconvolution
    
    @author Antoine Passemiers
    @version 1.0 27/05/2019
*/

#include "heuristic.hpp"

namespace wassersteinms {

void applyCorrection(Eigen::VectorXd &p) {
    size_t k = p.size();
    if ((p.array() < 0.0).any()) {
        p = p.cwiseAbs();
    }
    p = p.cwiseQuotient(Eigen::VectorXd::Constant(k, p.sum()));
}

double wassersteinDistance(Eigen::VectorXd &p,
                           std::unique_ptr<ProblemInstance> &prob) {
    Eigen::VectorXd g = prob->c.head(prob->n);
    Eigen::VectorXd d = -prob->b.head(prob->n - 1);
    return (prob->F.transpose() * p - g).cwiseAbs().cwiseProduct(d).sum();
}

Eigen::VectorXd geneticAlgorithm(
        std::unique_ptr<ProblemInstance> &prob,
        double epsilon,
        size_t nMaxIterations) {

    Eigen::VectorXd p = Eigen::VectorXd::Ones(prob->k);
    applyCorrection(p);
    Eigen::VectorXd bestP = p;
    double bestDistance = wassersteinDistance(p, prob);

    for (int i = 0; i < 100000; i++) {
        p = Eigen::VectorXd::Random(prob->k);
        applyCorrection(p);
        double distance = wassersteinDistance(p, prob);
        if (distance < bestDistance) {
            bestDistance = distance;
            bestP = p;
        }
    }

    std::cout << bestP << std::endl;
    return bestP;
}

} // namespace wassersteinms