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

//Stripyhoneycomb lattice as described in REF arXiv:1907.05423
//Note that crystal axis a1 a2 a3 correspond to d0 d1 d2, respectively in this implementation

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
struct stripyhoneycomb : bravais<Site, 3ul, 8ul> {
    using Base = bravais<Site, 3ul, 8ul>;
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
                iterator{it.cell_it().down(0).down(2), 7},  //x+ type
                iterator{it.cell_it().down(0), 7},          //y- type
                iterator{it.cell_it(), 1},                  //z  type
            };
        case 1:
        {
            return {
                iterator{it.cell_it(), 2},                  //x+ type
                iterator{it.cell_it().down(2), 2},          //y- type
                iterator{it.cell_it(), 0},                  //z  type
            };
        }
        case 2:
        {
            return {
                iterator{it.cell_it(), 1},                  //x+ type
                iterator{it.cell_it().up(2), 1},            //y- type
                iterator{it.cell_it(), 3},                  //z  type
            };
        }
        case 3:
        {
            return {
                iterator{it.cell_it(), 4},                  //x- type
                iterator{it.cell_it().down(1).up(2), 4},    //y+ type
                iterator{it.cell_it(), 2},                  //z  type
            };
        }
        case 4:
        {
            return {
                iterator{it.cell_it(), 3},                  //x- type
                iterator{it.cell_it().up(1).down(2), 3},    //y+ type
                iterator{it.cell_it(), 5},                  //z  type
            };
        }
        case 5:
        {
            return {
                iterator{it.cell_it().up(1).down(2), 6},    //x- type
                iterator{it.cell_it(), 6},                  //y+ type
                iterator{it.cell_it(), 4},                  //z  type
            };
        }
        case 6:
        {
            return {
                iterator{it.cell_it().down(1).up(2), 5},    //x- type
                iterator{it.cell_it(), 5},                  //y+ type
                iterator{it.cell_it(), 7},                  //z  type
            };
        }
        case 7:
        {
            return {
                iterator{it.cell_it().up(0).up(2), 0},      //x+ type
                iterator{it.cell_it().up(0), 0},            //y- type
                iterator{it.cell_it(), 6},                  //z  type
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
                const_iterator{it.cell_it().down(0).down(2), 7},  //x+ type
                const_iterator{it.cell_it().down(0), 7},          //y- type
                const_iterator{it.cell_it(), 1},                  //z  type
            };
        case 1:
        {
            return {
                const_iterator{it.cell_it(), 2},                  //x+ type
                const_iterator{it.cell_it().down(2), 2},          //y- type
                const_iterator{it.cell_it(), 0},                  //z  type
            };
        }
        case 2:
        {
            return {
                const_iterator{it.cell_it(), 1},                  //x+ type
                const_iterator{it.cell_it().up(2), 1},            //y- type
                const_iterator{it.cell_it(), 3},                  //z  type
            };
        }
        case 3:
        {
            return {
                const_iterator{it.cell_it(), 4},                  //x- type
                const_iterator{it.cell_it().down(1).up(2), 4},    //y+ type
                const_iterator{it.cell_it(), 2},                  //z  type
            };
        }
        case 4:
        {
            return {
                const_iterator{it.cell_it(), 3},                  //x- type
                const_iterator{it.cell_it().up(1).down(2), 3},    //y+ type
                const_iterator{it.cell_it(), 5},                  //z  type
            };
        }
        case 5:
        {
            return {
                const_iterator{it.cell_it().up(1).down(2), 6},    //x- type
                const_iterator{it.cell_it(), 6},                  //y+ type
                const_iterator{it.cell_it(), 4},                  //z  type
            };
        }
        case 6:
        {
            return {
                const_iterator{it.cell_it().down(1).up(2), 5},    //x- type
                const_iterator{it.cell_it(), 5},                  //y+ type
                const_iterator{it.cell_it(), 7},                  //z  type
            };
        }
        case 7:
        {
            return {
                const_iterator{it.cell_it().up(0).up(2), 0},      //x+ type
                const_iterator{it.cell_it().up(0), 0},            //y- type
                const_iterator{it.cell_it(), 6},                  //z  type
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
struct serializer<frustmag::lattice::stripyhoneycomb<Site>>
    : frustmag::lattice::serializer<frustmag::lattice::stripyhoneycomb<Site>> {};

}
}
