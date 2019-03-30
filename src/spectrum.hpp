/**
    spectrum.hpp
    Mass spectra
    
    @author Antoine Passemiers
    @version 1.0 30/03/2019
*/

#ifndef SPECTRUM_HPP_
#define SPECTRUM_HPP_

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <vector>
#include <numeric>
#include <iostream>

namespace wassersteinms {

// Type of m/z (frequency) ratios
typedef double mz_t;

// Type of m/z intensities
typedef double intensity_t;

// Maximum value for m/z intensities
const double MAX_INTENSITY = std::numeric_limits<double>::infinity();


class Spectrum {
    friend std::ostream& operator<<(std::ostream&, const Spectrum&);
private:

    // Peak intensities represented by a map object.
    // Keys refer to m/z ratios, and values refer to the associated intensities.
    std::map<mz_t, intensity_t> peaks;

public:

    explicit Spectrum(): peaks() {}
    explicit Spectrum(Spectrum &other): peaks(other.peaks) {}
    Spectrum(const Spectrum &other): peaks(other.peaks) {}
    Spectrum& operator=(const Spectrum& other) { this->peaks = other.peaks; return *this; }
    ~Spectrum() = default;

    intensity_t& operator[](const mz_t key) { return this->peaks[key]; }
    const intensity_t& operator[](const mz_t key) const { return this->peaks.at(key); }

    /**
        Normalize spectrum by dividing values by the sum over spectrum's values.

        @return A new Spectrum object with normalized values.
    */
    Spectrum normalize();

    /**
        Normalize spectrum by scaling values the following way:
        x_new = (x_old - x_min) / (x_max - x_min).

        @return A new Spectrum object with scaled values.
    */
    Spectrum scale();

    /**
        Change spectrum's resolution and merge m/z ratios that fall
        in the same newly created bins.

        @param resolution Resolution to use to form the new bins.
        @return A new Spectrum object with bins of given resolution.
    */
    Spectrum changeResolution(const double resolution);

    /**
        Remove a bin (inplace operation).

        @param key M/z ratio operation to be removed from the spectrum.
    */
    void remove(mz_t key) { this->peaks.erase(key); }

    /**
        Given another spectrum, add the bins that are not present in
        current spectrum.

        @param Another spectrum
    */
    void addKeys(Spectrum &other);

    std::map<mz_t, intensity_t>::iterator begin() { return this->peaks.begin(); }
    std::map<mz_t, intensity_t>::iterator end() { return this->peaks.end(); }
    std::map<mz_t, intensity_t>::const_iterator begin() const { return this->peaks.begin(); }
    std::map<mz_t, intensity_t>::const_iterator end() const { return this->peaks.end(); }
    const int size() const { return this->peaks.size(); }
};

} // namespace wassersteinms

#endif // SPECTRUM_HPP_