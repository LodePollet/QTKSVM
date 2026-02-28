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

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include <tksvm/config/policy.hpp>
#include <tksvm/phase_space/classifier/policy.hpp>
#include <tksvm/phase_space/sweep/policy.hpp>
#include <tksvm/sim_adapters/embarrassing_adapter.hpp>
#include <tksvm/utilities/convenience_params.hpp>

#include <frustmag/concepts.hpp>
//client-specific: include header for desired lattice
#include <frustmag/lattice/ortho.hpp>
//client-specific: include header for desired site-type
//#include <frustmag/site/v6.hpp>
#include <frustmag/site/spin_O3.hpp>


#include <client/config_policy.hpp>
#include <client/phase_point.hpp>

namespace client {

template<template<typename> typename Lattice>
class sim : public tksvm::embarrassing_adapter<tksvm::phase_space::point::phase_point> {

public:
    using phase_point = tksvm::phase_space::point::phase_point;
    using Base = tksvm::embarrassing_adapter<phase_point>;
    using parameters_type = typename Base::parameters_type;
    //client-specific: adapt site type
    //using site_type = tksvm::frustmag::site::v6;
    using site_type = tksvm::frustmag::site::spin_O3;
    using lattice_type = Lattice<site_type>;
    using site_iterator = typename lattice_type::iterator;
    using const_site_iterator = typename lattice_type::const_iterator;


private:
    std::string data_path;
    std::ifstream is;
    size_t sweeps;
    size_t total_sweeps;
    std::mt19937 rng;
    phase_point ppoint;
    std::vector<lattice_type> all_samples;


public:
    static void define_parameters(parameters_type & parameters) {
        if (parameters.is_restored()) {
            return;
        }
        Base::define_parameters(parameters);
        // Adds the convenience parameters (for save/load)
        // followed by simulation control parameters
        tksvm::define_convenience_parameters(parameters)
            .description("data consumer for quantum models")
            .define<std::string>("datapath", ".", "path to the data");

        //phase_point::define_parameters(parameters);
        define_config_policy_parameters(parameters);
    }

    sim(parameters_type & parms, std::size_t seed_offset)
        : Base(parms, seed_offset),
          data_path{parameters["datapath"].as<std::string>()},
          sweeps(0),
          rng(parameters["SEED"].as<std::size_t>() + seed_offset)
    {

        //phase_point pp{parameters};
        //update_phase_point(pp);

        measurements()
            //client-specific: add measurement names
            << alps::accumulators::FullBinningAccumulator<double>("OBS1")
            << alps::accumulators::FullBinningAccumulator<double>("OBS2");
    }

    virtual void update() override {
        //std::copy_n(std::istream_iterator<site_type>{is}, lattice.size(), lattice.begin());

        for (auto& lattice : all_samples) {
            for (auto it = lattice.begin(); it != lattice.end(); ++it) {

                //client-specific: choose POVM
                /* TETRA POVM Spin 1/2 Decker orientation, see ref [Decker03]
                double povm_outcome;
                is >> povm_outcome;
                if (povm_outcome == 0)      //Pt1
                    *it = site_type{(Eigen::Vector3d() <<  sqrt(2./3.), 0,  1./sqrt(3.)).finished()};
                else if (povm_outcome == 1) //Pt2
                    *it = site_type{(Eigen::Vector3d() << -sqrt(2./3.), 0,  1./sqrt(3.)).finished()};
                else if (povm_outcome == 2) //Pt3
                    *it = site_type{(Eigen::Vector3d() <<  0,  sqrt(2./3.), -1./sqrt(3.)).finished()};
                else if (povm_outcome == 3) //Pt4
                    *it = site_type{(Eigen::Vector3d() <<  0, -sqrt(2./3.), -1./sqrt(3.)).finished()};
                else
                    std::cout << "update() mapping: bad value encountered " << povm_outcome << std::endl;
                */

                //client-specific: choose POVM
                ///* Pauli-6 POVM Mapping (Spin 1/2)
                double povm_outcome;
                is >> povm_outcome;
                if (povm_outcome == 0)      //xup
                    *it = site_type{(Eigen::Vector3d() <<  +1,  0,  0).finished()};
                else if (povm_outcome == 1) //xdn
                    *it = site_type{(Eigen::Vector3d() <<  -1,  0,  0).finished()};
                else if (povm_outcome == 2) //yup
                    *it = site_type{(Eigen::Vector3d() <<   0, +1,  0).finished()};
                else if (povm_outcome == 3) //ydn
                    *it = site_type{(Eigen::Vector3d() <<   0, -1,  0).finished()};
                else if (povm_outcome == 4) //zup
                    *it = site_type{(Eigen::Vector3d() <<   0,  0, +1).finished()};
                else if (povm_outcome == 5) //zdn
                    *it = site_type{(Eigen::Vector3d() <<   0,  0, -1).finished()};
                else
                    std::cout << "update() mapping: bad value encountered " << povm_outcome << std::endl;
                //*/


                //client-specific: choose POVM
                /* SIC-POVM Mapping Spin 1
                double sq2 = std::sqrt(2.);
                double sq6 = std::sqrt(6.);
                double povm_outcome;
                is >> povm_outcome;
                if (povm_outcome == 0)
                    *it = site_type{(Eigen::Matrix<double, 6, 1>() <<     -sq2,     -sq6,       -2,        1,        1,        0).finished()};
                else if (povm_outcome == 1)
                    *it = site_type{(Eigen::Matrix<double, 6, 1>() <<        0,        0,        0,       -1,        1,        2).finished()};
                else if (povm_outcome == 2)
                    *it = site_type{(Eigen::Matrix<double, 6, 1>() <<     -sq2,     -sq6,        2,        1,        1,        0).finished()};
                else if (povm_outcome == 3)
                    *it = site_type{(Eigen::Matrix<double, 6, 1>() <<     -sq2,      sq6,       -2,        1,        1,        0).finished()};
                else if (povm_outcome == 4)
                    *it = site_type{(Eigen::Matrix<double, 6, 1>() <<        0,        0,        0,       -1,        1,        2).finished()};
                else if (povm_outcome == 5)
                    *it = site_type{(Eigen::Matrix<double, 6, 1>() <<     -sq2,      sq6,        2,        1,        1,        0).finished()};
                else if (povm_outcome == 6)
                    *it = site_type{(Eigen::Matrix<double, 6, 1>() <<    2*sq2,        0,       -2,        1,        1,        0).finished()};
                else if (povm_outcome == 7)
                    *it = site_type{(Eigen::Matrix<double, 6, 1>() <<        0,        0,        0,        2,       -2,        2).finished()};
                else if (povm_outcome == 8)
                    *it = site_type{(Eigen::Matrix<double, 6, 1>() <<    2*sq2,        0,        2,        1,        1,        0).finished()};
                else
                    std::cout << "update() mapping: bad value encountered " << povm_outcome << std::endl;
                */

                //client-specific: choose POVM
                /* MUB spin 1 map
                double s23 = std::sqrt(2./3.);
                double f23 = std::sqrt(2.)/3.;
                double povm_outcome;
                is >> povm_outcome;
                if (povm_outcome == 10)
                    *it = site_type{(Eigen::Matrix<double, 6, 1>() <<        0,        0,        4,        0,        0,        2).finished()};
                else if (povm_outcome == 11)
                    *it = site_type{(Eigen::Matrix<double, 6, 1>() <<        0,        0,        0,        2,        2,       -2).finished()};
                else if (povm_outcome == 12)
                    *it = site_type{(Eigen::Matrix<double, 6, 1>() <<        0,        0,       -4,        0,        0,        2).finished()};
                else if (povm_outcome == 13)
                    *it = site_type{(Eigen::Matrix<double, 6, 1>() <<    2*f23,    2*s23,        0,        0,    4./3.,    2./3.).finished()};
                else if (povm_outcome == 14)
                    *it = site_type{(Eigen::Matrix<double, 6, 1>() <<   -4*f23,        0,        0,        2,   -2./3.,    2./3.).finished()};
                else if (povm_outcome == 15)
                    *it = site_type{(Eigen::Matrix<double, 6, 1>() <<    2*f23,   -2*s23,        0,        0,    4./3.,    2./3.).finished()};
                else if (povm_outcome == 16)
                    *it = site_type{(Eigen::Matrix<double, 6, 1>() <<   -4*f23,        0,        0,        2,   -2./3.,    2./3.).finished()};
                else if (povm_outcome == 17)
                    *it = site_type{(Eigen::Matrix<double, 6, 1>() <<    2*f23,   -2*s23,        0,        0,    4./3.,    2./3.).finished()};
                else if (povm_outcome == 18)
                    *it = site_type{(Eigen::Matrix<double, 6, 1>() <<    2*f23,    2*s23,        0,        0,    4./3.,    2./3.).finished()};
                else if (povm_outcome == 19)
                    *it = site_type{(Eigen::Matrix<double, 6, 1>() <<   -4*f23,    4*s23,        0,        0,    4./3.,    2./3.).finished()};
                else if (povm_outcome == 20)
                    *it = site_type{(Eigen::Matrix<double, 6, 1>() <<   -4*f23,   -4*s23,        0,        0,    4./3.,    2./3.).finished()};
                else if (povm_outcome == 21)
                    *it = site_type{(Eigen::Matrix<double, 6, 1>() <<    8*f23,        0,        0,        2,   -2./3.,    2./3.).finished()};
                else
                    std::cout << "update() mapping: bad value encountered " << povm_outcome << std::endl;
                */
                
            }
        }
    }

    virtual void measure() override {
        return;
    }

    virtual double fraction_completed() const override {
        return 0;
    }

    using Base::save;
    using Base::load;

    virtual void save(alps::hdf5::archive & ar) const override {
        Base::save(ar);
        {
            std::stringstream engine_ss;
            engine_ss << rng;
            ar["checkpoint/random"] << engine_ss.str();
        }
        ar["checkpoint/total_sweeps"] << total_sweeps;
        //ar["checkpoint/lattice"] << lattice;
    }

    virtual void load(alps::hdf5::archive & ar) override {
        Base::load(ar);
        {
            std::string engine_str;
            ar["checkpoint/random"] >> engine_str;
            std::istringstream engine_ss(engine_str);
            engine_ss >> rng;
        }
        //ar["checkpoint/sweeps"] >> sweeps;
        //ar["checkpoint/lattice"] >> lattice;
    }

    std::vector<std::string> order_param_names() const {
        //client-specific: return names of observables
        return {"OBS1", "OBS2"};
    }

    virtual void reset_sweeps(bool) override {
        sweeps = 0;
        is.seekg(0);
    }

    bool is_thermalized() const {
        return true;
    }

    virtual phase_point phase_space_point() const override {
        return ppoint;
    }


    std::vector<lattice_type> configuration() const {
        return all_samples;
    }

    std::vector<lattice_type> random_configuration() {
        std::vector<lattice_type> random_samples(total_sweeps);
        for (auto& random_lattice : random_samples) {
            random_lattice = all_samples[0];
            using site_t = typename lattice_type::value_type;
            std::generate(random_lattice.begin(), random_lattice.end(),
                [&] { return site_t::random(rng); });
        }
        return random_samples;
    }

    virtual bool update_phase_point(phase_point const& pp) override {
        std::mt19937 rng{};
        bool changed = (pp != ppoint);
        if (changed) {
            ppoint = pp;
            std::string file_name = [&] {
                std::stringstream ss;
                std::string dataset;
                ///* Access datasets through temperature parameter
                dataset = "Run_" + std::to_string(int(ppoint.temperature() + 0.5));
                //*/
                ss << data_path << "/" << dataset << ".txt";
                return ss.str();
            }();
    #pragma omp critical
            std::clog << "opening file '" << file_name << "'\n";
            if (is.is_open())
                is.close();
            is.open(file_name);
            if (!is)
                throw std::runtime_error("could not open file: " + file_name);

            // determine line length
            std::string first_line;
            std::getline(is, first_line);
            std::stringstream ss{first_line};
            size_t n_line = std::distance(std::istream_iterator<double>{ss},
                                          std::istream_iterator<double>{});

            // determine total number of samples
            is.seekg(0, std::ios::end);
            total_sweeps = is.tellg() / (first_line.size() + 1) - 1;

            if (is.tellg() % (first_line.size() + 1) != 0)
    #pragma omp critical
            {
                std::cerr << "warning: file size not a multiple of first line\n"
                          << "file size: " << is.tellg() << "\t1st line: "
                          << (first_line.size() + 1) << '\n';
            }
            is.seekg(0);

            all_samples.resize(total_sweeps);
            for (auto& lattice : all_samples) {

                //client-specific: Infer lattice size from line length

                // Here we randomly initialize the lattice by infering its size from
                // the length of a line from the input data
                // The reading of the values is done in update()

                ///* DIM=1 (chain)
                lattice = {static_cast<size_t>(n_line), true, [&rng] {
                //*/
                    
                /* DIM=2 N_BASIS=2 (squarelink)
                lattice = {static_cast<size_t>(sqrt(n_line/2)), true, [&rng] {
                */
                    return site_type::random(rng);
                    }};
            }
        }

        return changed;
    }


    template <typename Introspector>
    using config_policy_type = tksvm::config::policy<lattice_type, Introspector>;

    template <typename Introspector>
    static auto config_policy_from_parameters(alps::params const& parameters,
                                              bool unsymmetrize = true)
        -> std::unique_ptr<config_policy_type<Introspector>> {
            return client::config_policy_from_parameters<lattice_type, Introspector>(
                parameters, unsymmetrize);
    }

};

}

namespace tksvm {
    //client-specific: specify which lattice to use
    using sim_base = client::sim<frustmag::lattice::chain>;
}

