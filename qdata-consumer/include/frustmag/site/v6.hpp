// SVM Order Parameters for Hidden Spin Order
// Copyright (C) 2018-2019  Jonas Greitemann, Ke Liu, and Lode Pollet

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#pragma once

#include <cmath>
#include <random>

#include <Eigen/Dense>

#include <frustmag/concepts.hpp>


namespace tksvm {
namespace frustmag {
namespace site {

struct v6 : Eigen::Matrix<double, 6, 1> {
    template <typename RNG>
    static v6 random(RNG & rng) {

        //client-specific: Choose POVM
        /* SIC-POVM mapping Spin 1
        double sq2 = std::sqrt(2.);
        double sq6 = std::sqrt(6.);
        int rndint = std::uniform_int_distribution<> {0, 8}(rng);
        v6 ret;
        if (rndint == 0)      ret <<     -sq2,     -sq6,       -2,        1,        1,        0;
        else if (rndint == 1) ret <<        0,        0,        0,       -1,        1,        2;
        else if (rndint == 2) ret <<     -sq2,     -sq6,        2,        1,        1,        0;
        else if (rndint == 3) ret <<     -sq2,      sq6,       -2,        1,        1,        0;
        else if (rndint == 4) ret <<        0,        0,        0,       -1,        1,        2;
        else if (rndint == 5) ret <<     -sq2,      sq6,        2,        1,        1,        0;
        else if (rndint == 6) ret <<    2*sq2,        0,       -2,        1,        1,        0;
        else if (rndint == 7) ret <<        0,        0,        0,        2,       -2,        2;
        else if (rndint == 8) ret <<    2*sq2,        0,        2,        1,        1,        0;
        else
            std::cout << "update() mapping: bad value encountered " << rndint << std::endl;
        return ret;
        */

        //client-specific: Choose POVM
        ///* MUB spin 1 map
        double s23 = std::sqrt(2./3.);
        double f23 = std::sqrt(2.)/3.;
        int rndint = std::uniform_int_distribution<> {10, 21}(rng);
        v6 ret;
        if (rndint == 10)      ret <<        0,        0,        4,        0,        0,        2;
        else if (rndint == 11) ret <<        0,        0,        0,        2,        2,       -2;
        else if (rndint == 12) ret <<        0,        0,       -4,        0,        0,        2;
        else if (rndint == 13) ret <<    2*f23,    2*s23,        0,        0,    4./3.,    2./3.;
        else if (rndint == 14) ret <<   -4*f23,        0,        0,        2,   -2./3.,    2./3.;
        else if (rndint == 15) ret <<    2*f23,   -2*s23,        0,        0,    4./3.,    2./3.;
        else if (rndint == 16) ret <<   -4*f23,        0,        0,        2,   -2./3.,    2./3.;
        else if (rndint == 17) ret <<    2*f23,   -2*s23,        0,        0,    4./3.,    2./3.;
        else if (rndint == 18) ret <<    2*f23,    2*s23,        0,        0,    4./3.,    2./3.;
        else if (rndint == 19) ret <<   -4*f23,    4*s23,        0,        0,    4./3.,    2./3.;
        else if (rndint == 20) ret <<   -4*f23,   -4*s23,        0,        0,    4./3.,    2./3.;
        else if (rndint == 21) ret <<    8*f23,        0,        0,        2,   -2./3.,    2./3.;
        else
            std::cout << "update() mapping: bad value encountered " << rndint << std::endl;
        return ret;
        //*/

        //client-specific: Generate uniformly distributed random numbers instead of POVM
        /* Random numbers to ensure that features average to zero
        return v6(Eigen::Matrix<double, 6, 1>::Random());
        */
    }

    v6() = default;

    v6(Eigen::Matrix<double, 6, 1> const& other) : Eigen::Matrix<double, 6, 1>{other} {}

    static const size_t size = 6;

    template <typename OutputIterator>
    OutputIterator & serialize(OutputIterator & it) const {
        *it = (*this)(0);
        *(++it) = (*this)(1);
        *(++it) = (*this)(2);
        *(++it) = (*this)(3);
        *(++it) = (*this)(4);
        *(++it) = (*this)(5);
        return ++it;
    }

    template <typename InputIterator>
    InputIterator & deserialize(InputIterator & it) {
        (*this)(0) = *it;
        (*this)(1) = *(++it);
        (*this)(2) = *(++it);
        (*this)(3) = *(++it);
        (*this)(4) = *(++it);
        (*this)(5) = *(++it);
        return ++it;
    }
};

static_assert(is_serializable<v6>::value, "v6 is not serializable");
static_assert(!is_archivable<v6>::value,
              "v6 is archivable, but shouldn't be");

#ifdef USE_CONCEPTS
static_assert(SiteState<v6>, "v6 is not a SiteState");
#endif

}
}
}
