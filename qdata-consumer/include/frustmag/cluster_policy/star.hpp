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

#include <frustmag/element_policy/common.hpp>


namespace tksvm {
namespace frustmag {
namespace cluster_policy {

    // NOTE: Star works for even linear length with obc arbitrary linear length with pbc
    //       Here we assume pbc in both directions, i.e. torus
    // These are overlapping star clusters
    template <typename Lattice>
    struct star {
        using ElementPolicy = element_policy::common<4, Lattice>;
        using bravais_const_iterator = typename Lattice::unitcell_const_iterator;

        star(ElementPolicy, Lattice const& lat)
            : lat{lat}
        {
        }

        struct unitcell;
        struct const_iterator {
            const_iterator & operator++ () { ++bit; return *this; }
            const_iterator operator++ (int) {
                const_iterator old(*this);
                ++(*this);
                return old;
            }
            const_iterator & operator-- () { --bit; return *this; }
            const_iterator operator-- (int) {
                const_iterator old(*this);
                --(*this);
                return old;
            }
            friend bool operator== (const_iterator lhs, const_iterator rhs) { return lhs.bit == rhs.bit; }
            friend bool operator!= (const_iterator lhs, const_iterator rhs) { return lhs.bit != rhs.bit; }

            unitcell operator* () const { return {bit}; }
            std::unique_ptr<unitcell> operator-> () const {
                return std::unique_ptr<unitcell>(new unitcell(bit));
            }
            friend star;
        private:
            const_iterator (bravais_const_iterator it) : bit {it} {}
            bravais_const_iterator bit;
        };

        struct unitcell {
            auto operator[](size_t i) const {
                switch (i) {
                    case 0: return (*bit)[0];
                    case 1: return (*bit)[1];
                    case 2: return (*(bit.down(0)))[0];
                    case 3: return (*(bit.down(1)))[1];
                    default:
                        throw std::runtime_error("In star cluster policy, unitcell operator[](): bad index");
                }
            }
            friend const_iterator;
        private:
            unitcell (bravais_const_iterator it) : bit {it} {}
            bravais_const_iterator bit;
        };

        const_iterator begin() const {
            return {lat.cellsbegin()};
        }

        const_iterator end() const {
            return {lat.cellsend()};
        }

        auto size() const {
            return lat.cells().size();
        }

    private:
        Lattice const& lat;
    };

}
}
}
