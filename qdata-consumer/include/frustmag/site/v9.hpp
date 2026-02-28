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

struct v9 : Eigen::Matrix<double, 9, 1> {
    template <typename RNG>
    static v9 random(RNG & rng) {

        ///* SIC-POVM mapping Spin 1 (POVM)
        double tsq3 = 2.*std::sqrt(3.);
        int rndint = std::uniform_int_distribution<> {0, 8}(rng);
        v9 ret;
        if (rndint == 0)      ret <<        0,        0,        0,        0,       -2,    -tsq3,      -1,       1,       1;
        else if (rndint == 1) ret <<        0,        0,       -2,     tsq3,        0,        0,       1,      -1,       1;
        else if (rndint == 2) ret <<       -2,    -tsq3,        0,        0,        0,        0,       1,       1,      -1;
        else if (rndint == 3) ret <<        0,        0,        0,        0,       -2,     tsq3,      -1,       1,       1;
        else if (rndint == 4) ret <<        0,        0,       -2,    -tsq3,        0,        0,       1,      -1,       1;
        else if (rndint == 5) ret <<       -2,     tsq3,        0,        0,        0,        0,       1,       1,      -1;
        else if (rndint == 6) ret <<        0,        0,        0,        0,        4,        0,      -1,       1,       1;
        else if (rndint == 7) ret <<        0,        0,        4,        0,        0,        0,       1,      -1,       1;
        else if (rndint == 8) ret <<        4,        0,        0,        0,        0,        0,       1,       1,      -1;
        else
            std::cout << "update() mapping: bad value encountered " << rndint << std::endl;
        return ret;
        //*/
    }

    v9() = default;

    v9(Eigen::Matrix<double, 9, 1> const& other) : Eigen::Matrix<double, 9, 1>{other} {}

    static const size_t size = 9;

    template <typename OutputIterator>
    OutputIterator & serialize(OutputIterator & it) const {
        *it = (*this)(0);
        *(++it) = (*this)(1);
        *(++it) = (*this)(2);
        *(++it) = (*this)(3);
        *(++it) = (*this)(4);
        *(++it) = (*this)(5);
        *(++it) = (*this)(6);
        *(++it) = (*this)(7);
        *(++it) = (*this)(8);
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
        (*this)(6) = *(++it);
        (*this)(7) = *(++it);
        (*this)(8) = *(++it);
        return ++it;
    }
};

static_assert(is_serializable<v9>::value, "v9 is not serializable");
static_assert(!is_archivable<v9>::value,
              "v9 is archivable, but shouldn't be");

#ifdef USE_CONCEPTS
static_assert(SiteState<v9>, "v9 is not a SiteState");
#endif

}
}
}
