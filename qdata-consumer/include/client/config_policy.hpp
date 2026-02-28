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

#include <memory>
#include <stdexcept>
#include <string>

#include <alps/params.hpp>

#include <tksvm/config/clustered_policy.hpp>
#include <tksvm/config/policy.hpp>
#include <tksvm/symmetry_policy/none.hpp>
#include <tksvm/symmetry_policy/symmetrized.hpp>

//client-specific: include clusters
#include <frustmag/cluster_policy/lattice.hpp>
#include <frustmag/cluster_policy/multicell.hpp>
#include <frustmag/cluster_policy/plaquette.hpp>
#include <frustmag/cluster_policy/single.hpp>
#include <frustmag/cluster_policy/star.hpp>
#include <frustmag/cluster_policy/fish.hpp>
#include <frustmag/cluster_policy/twocell_overlap.hpp>


namespace client {

template <typename Config, typename Introspector, typename SymmetryPolicy,
          typename ClusterPolicy>
using client_config_policy = tksvm::config::clustered_policy<Config,
                                                Introspector,
                                                SymmetryPolicy,
                                                ClusterPolicy>;

inline void define_config_policy_parameters(alps::params & parameters) {
    parameters
        .define<bool>("symmetrized", true, "use symmetry <S_x S_y> == <S_y S_x>")
        .define<std::string>("cluster", "lattice", "cluster used for SVM config")
        .define<size_t>("rank", "rank of the order parameter tensor");
}


template <typename Config, typename Introspector>
auto config_policy_from_parameters(alps::params const& parameters,
                                   bool unsymmetrize = true)
    -> std::unique_ptr<tksvm::config::policy<Config, Introspector>>
{
#define CONFPOL_CREATE()                                                        \
    return std::unique_ptr<tksvm::config::policy<Config, Introspector>>(        \
        new client_config_policy<Config, Introspector,                          \
            SymmetryPolicy, ClusterPolicy >(                                    \
                rank, typename ClusterPolicy::ElementPolicy{}, unsymmetrize));

#define CONFPOL_BRANCH_SYMM()                                                   \
    if (parameters["symmetrized"].as<bool>()) {                                 \
        using SymmetryPolicy = tksvm::symmetry_policy::symmetrized;             \
        CONFPOL_CREATE()                                                        \
    } else {                                                                    \
        using SymmetryPolicy = tksvm::symmetry_policy::none;                    \
        CONFPOL_CREATE()                                                        \
    }


    // set up SVM configuration policy
    size_t rank = parameters["rank"].as<size_t>();
    std::string clname = parameters["cluster"];
    //client-specific: one case for each cluster
    if (clname == "lattice") {
        using ClusterPolicy = tksvm::frustmag::cluster_policy::lattice<Config>;
        CONFPOL_BRANCH_SYMM()
    } else if (clname == "single"){
        using ClusterPolicy = tksvm::frustmag::cluster_policy::single<Config>;
        CONFPOL_BRANCH_SYMM()
    } else if (clname == "plaquette") {
        using ClusterPolicy = tksvm::frustmag::cluster_policy::plaquette<Config>;
        CONFPOL_BRANCH_SYMM()
    } else if (clname == "star") {
        using ClusterPolicy = tksvm::frustmag::cluster_policy::star<Config>;
        CONFPOL_BRANCH_SYMM()
    } else if (clname == "fish") {
        using ClusterPolicy = tksvm::frustmag::cluster_policy::fish<Config>;
        CONFPOL_BRANCH_SYMM()
    } else if (clname == "twocell_overlap") {
        using ClusterPolicy = tksvm::frustmag::cluster_policy::twocell_overlap<Config>;
        CONFPOL_BRANCH_SYMM()
    } else if (clname == "1cell"){
        using ClusterPolicy = tksvm::frustmag::cluster_policy::multicell<1,Config>;
        CONFPOL_BRANCH_SYMM()
    } else if (clname == "2cell"){
        using ClusterPolicy = tksvm::frustmag::cluster_policy::multicell<2,Config>;
        CONFPOL_BRANCH_SYMM()
    } else if (clname == "3cell"){
        using ClusterPolicy = tksvm::frustmag::cluster_policy::multicell<3,Config>;
        CONFPOL_BRANCH_SYMM()
    } else if (clname == "4cell"){
        using ClusterPolicy = tksvm::frustmag::cluster_policy::multicell<4,Config>;
        CONFPOL_BRANCH_SYMM()
    } else if (clname == "5cell"){
        using ClusterPolicy = tksvm::frustmag::cluster_policy::multicell<5,Config>;
        CONFPOL_BRANCH_SYMM()
    } else if (clname == "6cell"){
        using ClusterPolicy = tksvm::frustmag::cluster_policy::multicell<6,Config>;
        CONFPOL_BRANCH_SYMM()
    } else if (clname == "7cell"){
        using ClusterPolicy = tksvm::frustmag::cluster_policy::multicell<7,Config>;
        CONFPOL_BRANCH_SYMM()
    } else if (clname == "8cell"){
        using ClusterPolicy = tksvm::frustmag::cluster_policy::multicell<8,Config>;
        CONFPOL_BRANCH_SYMM()
    } else if (clname == "9cell"){
        using ClusterPolicy = tksvm::frustmag::cluster_policy::multicell<9,Config>;
        CONFPOL_BRANCH_SYMM()
    } else if (clname == "10cell"){
        using ClusterPolicy = tksvm::frustmag::cluster_policy::multicell<10,Config>;
        CONFPOL_BRANCH_SYMM()
    } else if (clname == "12cell"){
        using ClusterPolicy = tksvm::frustmag::cluster_policy::multicell<12,Config>;
        CONFPOL_BRANCH_SYMM()
    } else if (clname == "20cell"){
        using ClusterPolicy = tksvm::frustmag::cluster_policy::multicell<20,Config>;
        CONFPOL_BRANCH_SYMM()
    } else {
        throw std::runtime_error("Invalid cluster spec: " + clname);
    }
#undef CONFPOL_BRANCH_SYMM
#undef CONFPOL_CREATE
}

}
