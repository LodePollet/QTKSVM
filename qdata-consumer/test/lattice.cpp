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

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"

#include <algorithm>
#include <iostream>
#include <iterator>
#include <vector>

#include <frustmag/lattice/pyrochlore.hpp>
#include <frustmag/lattice/hyperhoneycomb.hpp>
#include <frustmag/lattice/stripyhoneycomb.hpp>
#include <frustmag/lattice/kagome.hpp>

struct int_site {
    int i;
};

template <typename Lattice>
void test_nn(Lattice && l, int nn_ind[][Lattice::coordination]) {
    int i = 0;
    for (auto it = l.begin(); it != l.end(); ++it, ++i) {
        CHECK(it->i == i);
        auto nn = l.nearest_neighbors(it);
        std::vector<int> nni;
        std::transform(nn.begin(), nn.end(), std::back_inserter(nni),
                       [end=l.end()] (auto it_n) {
                           return it_n == end ? -1 : it_n->i;
                       });
        //std::sort(nni.begin(), nni.end());
        bool nn_equal = std::equal(std::begin(nn_ind[i]), std::end(nn_ind[i]),
                                   nni.begin(), nni.end());
        if (!nn_equal) {
            std::copy(std::begin(nn_ind[i]), std::end(nn_ind[i]),
                      std::ostream_iterator<int>{std::cout, ", "});
            std::cout << std::endl;
            std::copy(nni.begin(), nni.end(),
                      std::ostream_iterator<int>{std::cout, ", "});
            std::cout << std::endl;

        }
        CHECK(nn_equal);
    }
}

auto increment_gen = [] {
    return [i=0]() mutable -> int_site { return {i++}; };
};

TEST_CASE("nn-pyro-periodic-1") {
    int nn_ind[][6] = {
        {1,2,3,1,2,3},
        {0,3,2,0,3,2},
        {3,0,1,3,0,1},
        {2,1,0,2,1,0},
    };
    test_nn(tksvm::frustmag::lattice::pyrochlore<int_site>(1, true, increment_gen()), nn_ind);
}

TEST_CASE("nn-pyro-periodic-2") {
    int nn_ind[][6] = {
        {1, 2, 3, 25, 22, 15},
        {0, 3, 2, 24, 23, 14},
        {3, 0, 1, 27, 20, 13},
        {2, 1, 0, 26, 21, 12},

        {5, 6, 7, 29, 18, 11},
        {4, 7, 6, 28, 19, 10},
        {7, 4, 5, 31, 16, 9 },
        {6, 5, 4, 30, 17, 8 },

        {9, 10, 11, 17, 30, 7},
        {8, 11, 10, 16, 31, 6},
        {11, 8, 9 , 19, 28, 5},
        {10, 9, 8 , 18, 29, 4},

        {13, 14, 15, 21, 26, 3},
        {12, 15, 14, 20, 27, 2},
        {15, 12, 13, 23, 24, 1},
        {14, 13, 12, 22, 25, 0},

        {17, 18, 19, 9 , 6 , 31},
        {16, 19, 18, 8 , 7 , 30},
        {19, 16, 17, 11, 4 , 29},
        {18, 17, 16, 10, 5 , 28},

        {21, 22, 23, 13, 2 , 27},
        {20, 23, 22, 12, 3 , 26},
        {23, 20, 21, 15, 0 , 25},
        {22, 21, 20, 14, 1 , 24},

        {25, 26, 27, 1 , 14, 23},
        {24, 27, 26, 0 , 15, 22},
        {27, 24, 25, 3 , 12, 21},
        {26, 25, 24, 2 , 13, 20},

        {29, 30, 31, 5 , 10, 19},
        {28, 31, 30, 4 , 11, 18},
        {31, 28, 29, 7 , 8 , 17},
        {30, 29, 28, 6 , 9 , 16},

    };
    test_nn(tksvm::frustmag::lattice::pyrochlore<int_site>(2, true, increment_gen()), nn_ind);
}


TEST_CASE("nn-hyperhoney-periodic-1") {
    int nn_ind[][3] = {
        {3,3,1},
        {2,2,0},
        {1,1,3},
        {0,0,2},
    };
    test_nn(tksvm::frustmag::lattice::hyperhoneycomb<int_site>(1, true, increment_gen()), nn_ind);
}


TEST_CASE("nn-hyperhoney-periodic-2") {
    int nn_ind[][3] = {
        {7,19,1}, {2,10,0}, {1,9,3}, {4,16,2},
        {3,23,5}, {6,14,4}, {5,13,7}, {0,20,6},
        {15,27,9}, {10,2,8}, {9,1,11}, {12,24,10},
        {11,31,13}, {14,6,12}, {13,5,15}, {8,28,14},
        {23,3,17}, {18,26,16}, {17,25,19}, {20,0,18},
        {19,7,21}, {22,30,20}, {21,29,23}, {16,4,22},
        {31,11,25}, {26,18,24}, {25,17,27}, {28,8,26},
        {27,15,29}, {30,22,28}, {29,21,31}, {24,12,30},
    };
    test_nn(tksvm::frustmag::lattice::hyperhoneycomb<int_site>(2, true, increment_gen()), nn_ind);
}


TEST_CASE("nn-stripyhoney-periodic-1") {
    int nn_ind[][3] = {
        {7,7,1},
        {2,2,0},
        {1,1,3},
        {4,4,2},
        {3,3,5},
        {6,6,4},
        {5,5,7},
        {0,0,6},
    };
    test_nn(tksvm::frustmag::lattice::stripyhoneycomb<int_site>(1, true, increment_gen()), nn_ind);
}

/* incomplete but present tests succeed
TEST_CASE("nn-stripyhoney-periodic-2") {
    int nn_ind[][3] = {
        {47,15, 1}, { 2,34, 0}, { 1,33, 3}, { 4,52, 2},
        { 3,51, 5}, {54, 6, 4}, {53, 5, 7}, {40, 8, 6},

        {39, 7, 9}, {10,42, 8}, { 9,41,11}, {12,60,10},
        {11,59,13}, {62,14,12}, {61,13,15}, {32, 0,14},

        {63,31,17}, {18,50,16}, {17,49,19}, {20,36,18},
        {19,35,21}, {38,22,20}, {37,21,23}, {56,24,22},

        {55,23,25}, {26,58,24}, {25,57,27}, {28,44,26},
        {27,43,29}, {46,30,28}, {45,29,31}, {48,16,30},

    };
    test_nn(tksvm::frustmag::lattice::stripyhoneycomb<int_site>(2, true, increment_gen()), nn_ind);
}
*/


/* incomplete but present tests succeed
TEST_CASE("nn-hyperhoney-periodic-3") {
    int nn_ind[][3] = {
        {1,75,11}, {0,2,26}, {3,1,13}, {2,36,4},
    };
    test_nn(tksvm::frustmag::lattice::hyperhoneycomb<int_site>(3, true, increment_gen()), nn_ind);
}
*/

/* incomplete but present tests succeed
TEST_CASE("nn-stripyhoney-periodic-3") {
    int nn_ind[][3] = {
        {167,23,1}, {2,146,0}, {1,73,3}, {4,124,2},
        {3,171,5}, {174,6,4}, {125,5,7}, {80,8,6},
    };
    test_nn(tksvm::frustmag::lattice::stripyhoneycomb<int_site>(3, true, increment_gen()), nn_ind);
}
*/



/* incomplete but present tests succeed
TEST_CASE("nn-pyro-periodic-3") {
    int nn_ind[][6] = {
        {1, 2, 3, 85, 78, 19},
        {0, 3, 2, 60, 43, 30},
        {3, 0, 1, 51, 44, 21},
        {2, 1, 0, 98, 81, 32},
    };
    test_nn(tksvm::frustmag::lattice::pyrochlore<int_site>(3, true, increment_gen()), nn_ind);
}
*/
