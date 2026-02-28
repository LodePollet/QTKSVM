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

//Hyperhoneycomb lattice as described in REF arXiv:1907.05423
//Note that crystal axis a1 a2 a3 correspond to d0 d2 d1, respectively in this implementation

#pragma once

#include <array>

#include <tksvm/config/serializer.hpp>

#include <frustmag/concepts.hpp>
#include <frustmag/lattice/bravais.hpp>
#include <frustmag/lattice/serializer.hpp>


namespace tksvm {
namespace frustmag {
namespace lattice {

template <typename Site>
struct hyperhoneycomb : bravais<Site, 3ul, 4ul> {
    using Base = bravais<Site, 3ul, 4ul>;
    using iterator = typename Base::iterator;
    using const_iterator = typename Base::const_iterator;

    static const size_t coordination = 3;

    using Base::Base;

    auto nearest_neighbors(iterator const& it)
        -> std::array<iterator, coordination>
    {   
        switch (it.basis_index()) {
        case 0:
            return {
                iterator{it.cell_it().down(0), 3},      //x- type
                iterator{it.cell_it().down(2), 3},      //y+ type
                iterator{it.cell_it(), 1},              //z  type
            };
        case 1:
        {
            return {
                iterator{it.cell_it(), 2},              //x+ type
                iterator{it.cell_it().down(1), 2},      //y- type
                iterator{it.cell_it(), 0},              //z  type
            };
        }
        case 2:
        {
            return {
                iterator{it.cell_it(), 1},              //x+ type
                iterator{it.cell_it().up(1), 1},        //y- type
                iterator{it.cell_it(), 3},              //z  type
            };
        }
        case 3:
        {
            return {
                iterator{it.cell_it().up(0), 0},        //x- type
                iterator{it.cell_it().up(2), 0},        //y+ type
                iterator{it.cell_it(), 2},              //z  type
            };
        }
        default:
            return {};
        }
    }

    auto nearest_neighbors(const_iterator const& it) const
        -> std::array<const_iterator, coordination>
    {   
        switch (it.basis_index()) {
        case 0:
            return {
                const_iterator{it.cell_it().down(0), 3},      //x- type
                const_iterator{it.cell_it().down(2), 3},      //y+ type
                const_iterator{it.cell_it(), 1},              //z  type
            };
        case 1:
        {
            return {
                const_iterator{it.cell_it(), 2},              //x+ type
                const_iterator{it.cell_it().down(1), 2},      //y- type
                const_iterator{it.cell_it(), 0},              //z  type
            };
        }
        case 2:
        {
            return {
                const_iterator{it.cell_it(), 1},              //x+ type
                const_iterator{it.cell_it().up(1), 1},        //y- type
                const_iterator{it.cell_it(), 3},              //z  type
            };
        }
        case 3:
        {
            return {
                const_iterator{it.cell_it().up(0), 0},        //x- type
                const_iterator{it.cell_it().up(2), 0},        //y+ type
                const_iterator{it.cell_it(), 2},              //z  type
            };
        }
        default:
            return {};
        }
    }
    
};

}
}

namespace config {

template <typename Site>
struct serializer<frustmag::lattice::hyperhoneycomb<Site>>
    : frustmag::lattice::serializer<frustmag::lattice::hyperhoneycomb<Site>> {};

}
}
