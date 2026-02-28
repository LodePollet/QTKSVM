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

//This cluster is a collection of lin_l ^ Lattice::dim unit cells
template <size_t lin_l, typename Lattice>
struct multicell {
    using ElementPolicy = element_policy::common
                            <Lattice::n_basis*static_cast<size_t>(pow(lin_l, Lattice::dim)), Lattice>;
    using unitcell_const_iterator = typename Lattice::unitcell_collection_type::const_iterator;
    using difference_type = typename unitcell_const_iterator::difference_type;

    struct unitcell;

    struct const_iterator {

        const_iterator & operator++ () {
            /* XXX Non-overlapping clusters for any lattice and dimension
            auto uidx_old = uidx;
            auto incr = [&](int i, auto&& incr) {
                if(i == 0) {uidx += lin_l; return;}
                auto upper = plengths[i];
                auto lower = plengths[i-1];
                size_t shift = i-1 > 0 ? lower - lin_l*plengths[i-2] : 0;
                if ((uidx + lin_l*lower - shift) % upper == 0) {
                    uidx += lin_l*lower - shift;
                    uidx = (uidx/upper + lin_l-1) * upper;
                    return;
                }
                incr(i - 1, incr);
            };
            incr(Lattice::dim - 1,incr);
            uit += uidx - uidx_old;
            return *this;
            */

            ///*XXX Overlapping clusters, hardcoded for chain lattice with open boundary conditions
            auto L = plengths[1];
            auto uidx_old = uidx;
            if (uidx >= int(L) - int(lin_l))
                uidx = L;
            else
                uidx++;
            uit += uidx - uidx_old;
            return *this;
            //*/
        }

        const_iterator operator++ (int) {
            const_iterator old(*this);
            ++(*this);
            return old;
        }

        const_iterator & operator-- () {
            /*XXX Non-overlapping clusters for any lattice and dimension
            auto uidx_old = uidx;
            auto decr = [&](int i, auto&& decr) {
                if(i == 0) {uidx -= lin_l; return;}
                auto upper = plengths[i];
                auto lower = plengths[i-1];
                size_t shift = i-1 > 0 ? (lin_l-1)*lower : 0;
                if (uidx % (lin_l* upper) == 0) {
                    uidx += -(lin_l-1)*upper - lin_l - shift;
                    return;
                }
                decr(i - 1, decr);
            };
            decr(Lattice::dim - 1,decr);
            uit += uidx - uidx_old;
            return *this;
            */

            ///*XXX Overlapping clusters, hardcoded for chain lattice with open boundary conditions
            auto L = plengths[1];
            auto uidx_old = uidx;
            if (uidx >= int(L))
                uidx = L - lin_l;
            else
                uidx--;
            uit += uidx - uidx_old;
            return *this;
            //*/
        }

        const_iterator operator-- (int) {
            const_iterator old(*this);
            --(*this);
            return old;
        }

        friend bool operator== (const_iterator lhs, const_iterator rhs) { return lhs.uit == rhs.uit; }
        friend bool operator!= (const_iterator lhs, const_iterator rhs) { return lhs.uit != rhs.uit; }

        unitcell operator* () const { return {uit, plengths}; }

        std::unique_ptr<unitcell> operator-> () const {
            return std::unique_ptr<unitcell>(new unitcell(uit, plengths));
        }

        friend multicell;

    private:
        const_iterator(unitcell_const_iterator uit, difference_type uidx, 
                       size_t const* plengths)
            : uit(uit), uidx(uidx), plengths(plengths) {}

        unitcell_const_iterator uit;
        difference_type uidx;
        size_t const* plengths;
    };

    struct unitcell {
        auto operator[](size_t i) const {
            size_t L = plengths[1];
            size_t Lsq = pow(L,2.);
            size_t lin_lsq = pow(lin_l,2.);
            size_t ci = i/Lattice::n_basis;
            size_t base_index = ci / lin_lsq * Lsq + 
                                ci % lin_lsq / lin_l * L +
                                ci % lin_lsq % lin_l;
            return (*(uit + base_index))[i%Lattice::n_basis];
        }
        friend const_iterator;
    private:
        unitcell (unitcell_const_iterator uit, size_t const* plengths)
            : uit(uit), plengths(plengths) {}
        unitcell_const_iterator uit;
        size_t const* plengths;
    };


    multicell(ElementPolicy, Lattice const& lat)
        : ucells{lat.cells()}, plengths{lat.get_plengths()} {
        /*XXX This condition must only be fulfilled in case of non-verlapping clusters
        if (ucells.size() % (int)pow(lin_l, Lattice::dim) != 0)
                throw std::runtime_error(
                            "multicell cluster: linear length is incompatible with system size");
        */
        if (Lattice::dim < 1)
                throw std::runtime_error("lattice dimension must be non-zero");
    }

    const_iterator begin() const {
        return const_iterator{ucells.begin(), 0, plengths.data()};
    }

    //points to one after last unitcells
    const_iterator end() const {
        return const_iterator{ucells.end(), static_cast<difference_type>(ucells.size()), 
                              plengths.data()};
    }

    //number of clusters
    auto size() const {
        //XXX Non-verlapping clusters
        //return ucells.size() / pow(lin_l, Lattice::dim);
        //XXX Overlapping clusters, hardcoded for chain lattice with open boundary conditions
        return ucells.size() - (lin_l - 1);
    }

private:
    typename Lattice::unitcell_collection_type const& ucells;
    typename Lattice::lengths_type const& plengths;
};

}
}
}
