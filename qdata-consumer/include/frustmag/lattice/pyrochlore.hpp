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

#include <array>

#include <tksvm/config/serializer.hpp>

#include <frustmag/concepts.hpp>
#include <frustmag/lattice/bravais.hpp>
#include <frustmag/lattice/serializer.hpp>


namespace tksvm {
namespace frustmag {
namespace lattice {

template <typename Site>
struct pyrochlore : bravais<Site, 3ul, 4ul> {
    using Base = bravais<Site, 3ul, 4ul>;
    using iterator = typename Base::iterator;
    using const_iterator = typename Base::const_iterator;

    static const size_t coordination = 6;

    using Base::Base;

    auto nearest_neighbors(iterator const& it)
        -> std::array<iterator, coordination>
    {
        switch (it.basis_index()) {
        case 0:
            return {
                iterator{it.cell_it(), 1},
                iterator{it.cell_it(), 2},
                iterator{it.cell_it(), 3},
                iterator{it.cell_it().up(1).down(2), 1},
                iterator{it.cell_it().up(0).down(2), 2},
                iterator{it.cell_it().up(0).up(1), 3}
            };
        case 1:
        {
            return {
                iterator{it.cell_it(), 0},
                iterator{it.cell_it(), 3},
                iterator{it.cell_it(), 2},
                iterator{it.cell_it().down(1).up(2), 0},
                iterator{it.cell_it().up(0).up(2), 3},
                iterator{it.cell_it().up(0).down(1), 2}
            };
        }
        case 2:
        {
            return {
                iterator{it.cell_it(), 3},
                iterator{it.cell_it(), 0},
                iterator{it.cell_it(), 1},
                iterator{it.cell_it().up(1).up(2), 3},
                iterator{it.cell_it().down(0).up(2), 0},
                iterator{it.cell_it().down(0).up(1), 1}
            };
        }
        case 3:
        {
            return {
                iterator{it.cell_it(), 2},
                iterator{it.cell_it(), 1},
                iterator{it.cell_it(), 0},
                iterator{it.cell_it().down(1).down(2), 2},
                iterator{it.cell_it().down(0).down(2), 1},
                iterator{it.cell_it().down(0).down(1), 0}
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
                const_iterator{it.cell_it(), 1},
                const_iterator{it.cell_it(), 2},
                const_iterator{it.cell_it(), 3},
                const_iterator{it.cell_it().up(1).down(2), 1},
                const_iterator{it.cell_it().up(0).down(2), 2},
                const_iterator{it.cell_it().up(0).up(1), 3}
            };
        case 1:
        {
            return {
                const_iterator{it.cell_it(), 0},
                const_iterator{it.cell_it(), 3},
                const_iterator{it.cell_it(), 2},
                const_iterator{it.cell_it().down(1).up(2), 0},
                const_iterator{it.cell_it().up(0).up(2), 3},
                const_iterator{it.cell_it().up(0).down(1), 2}
            };
        }
        case 2:
        {
            return {
                const_iterator{it.cell_it(), 3},
                const_iterator{it.cell_it(), 0},
                const_iterator{it.cell_it(), 1},
                const_iterator{it.cell_it().up(1).up(2), 3},
                const_iterator{it.cell_it().down(0).up(2), 0},
                const_iterator{it.cell_it().down(0).up(1), 1}
            };
        }
        case 3:
        {
            return {
                const_iterator{it.cell_it(), 2},
                const_iterator{it.cell_it(), 1},
                const_iterator{it.cell_it(), 0},
                const_iterator{it.cell_it().down(1).down(2), 2},
                const_iterator{it.cell_it().down(0).down(2), 1},
                const_iterator{it.cell_it().down(0).down(1), 0}
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
struct serializer<frustmag::lattice::pyrochlore<Site>>
    : frustmag::lattice::serializer<frustmag::lattice::pyrochlore<Site>> {};

}
}
