#ifndef IPM_HPP__
#define IPM_HPP__

#include <Eigen/Core>
#include <Eigen/Dense>


typedef struct _ipmSolution {
    Eigen::VectorXd x;
    Eigen::VectorXd y;
    Eigen::VectorXd z;
} ipmSolution;


typedef struct _problemInstance {
    size_t n;
    size_t k;
    Eigen::MatrixXd A;
    Eigen::VectorXd b;
    Eigen::VectorXd c;
} problemInstance;


ipmSolution interiorPointMethod(ipmSolution sol0, problemInstance prob, double epsilon);


#endif // IPM_HPP__