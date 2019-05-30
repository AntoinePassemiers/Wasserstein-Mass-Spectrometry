/**
    heuristic.cpp
    Heuristic method for mass spectrum deconvolution
    
    @author Antoine Passemiers
    @version 1.0 27/05/2019
*/

#include "heuristic.hpp"

namespace wassersteinms {


int argmin(const Eigen::VectorXd &v) {
    Eigen::VectorXd::Index minIndex;
    v.minCoeff(&minIndex);
    return minIndex;
}

int argmin(const Eigen::VectorXd &v, const Eigen::VectorXi &indices) {
    double bestValue = v[indices[0]];
    int bestIndex = indices[0];
    for (int i = 1; i < indices.size(); i++) {
        if (v[indices[i]] < bestValue) {
            bestIndex = indices[i];
            bestValue = v[indices[i]];
        }
    }
    return bestIndex;
}

int argmax(const Eigen::VectorXd &v) {
    Eigen::VectorXd::Index maxIndex;
    v.maxCoeff(&maxIndex);
    return maxIndex;
}

void applyCorrection(Eigen::VectorXd &p) {
    size_t k = p.size();
    if ((p.array() < 0.0).any()) {
        p = p.cwiseAbs();
    }
    p = p.cwiseQuotient(Eigen::VectorXd::Constant(k, p.sum()));
}

Eigen::VectorXd randomSolution(size_t k) {
    Eigen::VectorXd p = Eigen::VectorXd::Random(k);
    // Each component is set to 0 with a 50% chance
    Eigen::VectorXd _0 = Eigen::VectorXd::Zero(k);
    p = (p.array() < 0.0).select(_0, p);
    applyCorrection(p);
    return p;
}

void mutationOperator(Eigen::VectorXd &p) {
    Eigen::VectorXd rand = Eigen::VectorXd::Random(p.size());
    double mutation = 2.0 * rand[0] / p.size();
    p[argmax(rand)] += mutation;
    applyCorrection(p);
}

void crossoverOperator(const Eigen::VectorXd &pa,
                       const Eigen::VectorXd &pb,
                       Eigen::VectorXd &pc) {
    pc = pa + pb;
}

double wassersteinDistance(const Eigen::VectorXd &p,
                           std::unique_ptr<ProblemInstance> &prob) {
    Eigen::VectorXd g = prob->c.head(prob->n);
    Eigen::VectorXd d = -prob->b.head(prob->n - 1);
    return (prob->F.transpose() * p - g).cwiseAbs().cwiseProduct(d).sum();
}

Eigen::VectorXd geneticAlgorithm(
        std::unique_ptr<ProblemInstance> &prob,
        double epsilon,
        size_t nMaxIterations) {

    // Create initial population
    int populationSize = 1000;
    int partitionSize = 50;
    assert(partitionSize * 2 <= populationSize);
    Eigen::VectorXd distances = Eigen::VectorXd::Ones(populationSize);
    std::vector<Eigen::VectorXd> population;
    for (int i = 0; i < populationSize; i++) {
        population.push_back(randomSolution(prob->k));
        distances[i] = wassersteinDistance(population[i], prob);
    }

    // Find best initial solution and initialize indices
    int i = argmin(distances);
    Eigen::VectorXd bestP = population[i];
    double bestDistance = distances[i];
    Eigen::VectorXi indices = Eigen::VectorXi::Ones(populationSize);
    for (size_t i = 0; i < populationSize; i++) indices[i] = i;

    for (size_t j = 0; j < nMaxIterations; j++) {

        // Shuffle the population and find two parents
        std::random_shuffle(
                indices.array().data(),
                indices.array().data() + indices.size());
        Eigen::VectorXi partitionA = indices.head(partitionSize);
        Eigen::VectorXi partitionB = indices.tail(partitionSize);
        int ia = argmin(distances, partitionA);
        int ib = argmin(distances, partitionB);

        // Find worst solution
        i = argmax(distances);

        // Apply crossover operator
        crossoverOperator(population[ia], population[ib], population[i]);

        // Apply mutation operator
        mutationOperator(population[i]);

        // Evaluate new solution
        distances[i] = wassersteinDistance(population[i], prob);
        if (distances[i] < bestDistance) {
            bestDistance = distances[i];
            bestP = population[i];
        }
    }

    return bestP;
}

} // namespace wassersteinms