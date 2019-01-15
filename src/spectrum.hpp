#ifndef SPECTRUM_HPP_
#define SPECTRUM_HPP_

#include <vector>


typedef double mz_t;
typedef double intensity_t;


class Spectrum {
private:
    bool sorted;
    std::vector<mz_t> ratios;
    std::vector<intensity_t> intensities;
    std::vector<intensity_t> relIntensities;
public:
    Spectrum() : sorted(true) {}
    void add(mz_t ratio, intensity_t intensity, intensity_t relIntensity);
    void sort();
};

#endif // SPECTRUM_HPP_