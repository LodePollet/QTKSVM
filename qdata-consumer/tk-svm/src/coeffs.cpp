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

#include <algorithm>
#include <exception>
#include <fstream>
#include <iostream>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <argh.h>

#include <alps/hdf5.hpp>
#include <alps/params.hpp>

#include <boost/multi_array.hpp>

#include <colormap/colormap.hpp>

#include <Eigen/SVD>

#include <svm/svm.hpp>
#include <svm/serialization/hdf5.hpp>

#include <tksvm/config_sim_base.hpp>
#include <tksvm/config/block_policy.hpp>
#include <tksvm/element_policy/components.hpp>
#include <tksvm/phase_space/classifier.hpp>
#include <tksvm/utilities/contraction.hpp>
#include <tksvm/utilities/filesystem.hpp>
#include <tksvm/utilities/matrix_output.hpp>
#include <tksvm/utilities/results/results.hpp>


using kernel_t = svm::kernel::linear;
using namespace tksvm;


int main(int argc, char** argv) {
    try {
        argh::parser cmdl;
        cmdl.add_params({"block", "t", "transition", "result"});
        cmdl.parse(argc, argv, argh::parser::SINGLE_DASH_IS_MULTIFLAG);
        alps::params parameters(argc, argv);

        sim_base::define_parameters(parameters);

        if (parameters.help_requested(std::cout) ||
            parameters.has_missing(std::cout)) {
            return 1;
        }

        using phase_point = typename sim_base::phase_point;
        auto phase_classifier = phase_space::classifier::from_parameters<phase_point>(
            parameters, "classifier.");

        std::string arname = parameters.get_archive_name();
        bool verbose = cmdl[{"-v", "--verbose"}] || cmdl[{"-c", "--contraction-weights"}];
        auto log_msg = [verbose] (std::string const& msg) {
            if (verbose)
                std::cout << msg << std::endl;
        };

        log_msg("Reading model...");
        using label_t = typename phase_space::classifier::policy<phase_point>::label_type;
        using model_t = svm::model<kernel_t, label_t>;
        using classifier_t = typename model_t::classifier_type;
        using introspec_t = svm::linear_introspector<classifier_t>;
        model_t model;
        {
            alps::hdf5::archive ar(arname, "r");
            svm::serialization::model_serializer<svm::hdf5_tag, model_t> serial(model);
            ar["model"] >> serial;
        }

        auto treat_transition = [&] (classifier_t const& classifier,
                                     std::string const& basename)
        {
            introspec_t coeff = svm::linear_introspect(classifier);

            auto confpol = sim_base::config_policy_from_parameters<introspec_t>(
                parameters,
                cmdl[{"-u", "--unsymmetrize"}]);

            auto contractions = get_contractions(confpol->rank());
            auto block_inds = confpol->all_block_indices();
            using block_ind_t = decltype(block_inds)::value_type;


            log_msg("Allocating coeffs...");
            boost::multi_array<double,1> coeffs(boost::extents[model.dim()]);
            log_msg("Filling coeffs...");
#pragma omp parallel for
            for (size_t i = 0; i < model.dim(); ++i) {
                coeffs[i] = coeff.coefficient(i);
            }
            std::vector<block_ind_t> block_inds_vec(block_inds.begin(),
                                                    block_inds.end());

            log_msg("Normalizing coeffs...");
            normalize_vector(coeffs);
            log_msg("Writing coeffs...");
            write_vector(coeffs, replace_extension(basename, ".coeffs"));

        };

        // Determine requested transitions
        auto transitions = model.classifiers();
        size_t t;
        bool exclusive = bool(cmdl({"-t", "--transition"}) >> t);
        for (size_t k = 0; k < transitions.size(); ++k) {
            if (exclusive && t != k)
                continue;
            auto const& cl = transitions[k];
            std::cout << k << ":   "
                      << phase_classifier->name(cl.labels().first) << " -- "
                      << phase_classifier->name(cl.labels().second)
                      << "\t rho = " << cl.rho() << std::endl;
            std::stringstream ss;
            ss << replace_extension(arname, "")
               << '-' << phase_classifier->name(cl.labels().first)
               << '-' << phase_classifier->name(cl.labels().second);
            if (!cmdl[{"-l", "--list"}])
                treat_transition(cl, ss.str());
        }

        return 0;
    } catch (const std::exception& exc) {
        std::cout << "Exception caught: " << exc.what() << std::endl;
        return 2;
    }
}
