#ifndef WASSERSTEIN_HPP__
#define WASSERSTEIN_HPP__

#include "spectrum.hpp"

#include <algorithm>
#include <cmath>
#include <memory>


double wassersteinDistance(
        std::unique_ptr<Spectrum> &S1,
        std::unique_ptr<Spectrum> &S2);


#endif // WASSERSTEIN_HPP__