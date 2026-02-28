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
#include <utility>

#include <alps/mc/mcbase.hpp>

#include <svm/svm.hpp>
#include <svm/serialization/hdf5.hpp>

#include <tksvm/phase_space/classifier.hpp>
#include <tksvm/phase_space/sweep.hpp>
#include <tksvm/utilities/void_t.hpp>


namespace tksvm {

namespace detail {

    template <typename T, typename = void, typename = void>
    struct empty_checker {
        T const& c;
        bool empty() const&& {
            return std::distance(std::begin(c), std::end(c)) == 0;
        }
    };

    template <typename T, typename U>
    struct empty_checker<T, void_t<typename std::enable_if<std::is_same<decltype(std::declval<T const&>().size()), size_t>::value>::type>, U> {
        T const& c;
        bool empty() const&& {
            return c.size() == 0;
        }
    };

    template <typename T>
    struct empty_checker<T, void_t<typename std::enable_if<std::is_same<decltype(std::declval<T const&>().empty()), bool>::value>::type>, void> {
        T const& c;
        bool empty() const&& {
            return c.empty();
        }
    };

}

template <class Simulation>
class training_adapter : public Simulation {
public:
    typedef alps::mcbase::parameters_type parameters_type;

    using phase_point = typename Simulation::phase_point;
    using label_t = typename phase_space::classifier::policy<phase_point>::label_type;

    using kernel_t = svm::kernel::linear;
    using problem_t = svm::problem<kernel_t, phase_point>;
    using model_t = svm::model<kernel_t, label_t>;
    using introspec_t = svm::linear_introspector<typename model_t::classifier_type>;

    using config_policy_t = typename Simulation::template config_policy_type<introspec_t>;
    using config_array = typename config_policy_t::config_array;

    static void define_parameters(parameters_type & parameters) {
        // If the parameters are restored, they are already defined
        if (parameters.is_restored()) {
            return;
        }

        // Adds the parameters of the base class
        Simulation::define_parameters(parameters);
        phase_space::classifier::define_parameters<phase_point>(parameters, "classifier.");
        phase_space::sweep::define_parameters<phase_point>(parameters, "sweep.");
        parameters
            .define<size_t>("sweep.samples", 1000,
                            "number of configuration samples taken"
                            " at each phase point")
            .define<size_t>("sweep.Nc", 1, "number of configuration"
                            " used to construct one feature vector (EAGER MODE)")
            ;
    }

    training_adapter(parameters_type & parms,
        std::size_t seed_offset = 0)
        : Simulation(parms, seed_offset)
        , confpol(Simulation::template config_policy_from_parameters<introspec_t>(parms))
        , N_sample(size_t(parameters["sweep.samples"]))
        , Nc(size_t(parameters["sweep.Nc"]))
        , problem(confpol->size())
        , prob_serializer(problem)
    {
    }

    virtual void measure () final override {
        //double frac = Simulation::fraction_completed();
        //Simulation::measure();
        std::vector<typename Simulation::lattice_type> config = Simulation::configuration();
        using detail::empty_checker;
        if (!empty_checker<decltype(config)>{config}.empty()) {
            sample_config(config, Simulation::phase_space_point());
        }
    }

    virtual void save (alps::hdf5::archive & ar) const override {
        Simulation::save(ar);

        // non-overridable parameters
        ar["training/N_sample"] << N_sample;

        // state
        //ar["training/i_sample"] << i_sample;

        if (problem.size() > 0)
            ar["training/problem"] << prob_serializer;
    }

    virtual void load (alps::hdf5::archive & ar) override {
        Simulation::load(ar);

        // non-overridable parameters
        ar["training/N_sample"] >> N_sample;

        // state
        //ar["training/i_sample"] >> i_sample;

        if (ar.is_group("training/problem"))
            ar["training/problem"] >> prob_serializer;
        else
            problem = {confpol->size()};
        if (problem.dim() != confpol->size())
            throw std::runtime_error("invalid problem dimension");
    }

    problem_t surrender_problem () {
        problem_t other_problem(confpol->size());
        std::swap(other_problem, problem);
        return other_problem;
    }

    virtual void sample_config(std::vector<typename Simulation::lattice_type> const& config,
                               phase_point const& ppoint)
    {
        if (Nc == 1) {
            #pragma omp parallel for
            for (auto lit = config.cbegin(); lit != config.cbegin()+N_sample; ++lit) {
                auto mapped_sample = confpol->configuration(*lit);
                #pragma omp critical
                problem.add_sample(mapped_sample, ppoint);
            }
        }
        else if (Nc > 1) {
            // number of constructed feature vectors will be approximaltely sweep.samples/sweep.Nc
            // the actual number depends on number of OMP threads
            // in general there will be at least (sweep.samples//numOMPthreads)//sweep.Nc * sweep.Nc feature vectors
            // '//' means integer division here and '/' means proper division
            #pragma omp parallel
            {
                size_t sample_counter = 0;
                std::vector<double> cumul_sample(confpol->size(), 0.);
                #pragma omp for
                for (auto lit = config.cbegin(); lit != config.cbegin()+N_sample; ++lit) {
                    auto mapped_sample = confpol->configuration(*lit);
                    std::transform(mapped_sample.begin(), mapped_sample.end(),
                                    cumul_sample.begin(), cumul_sample.begin(),
                                    std::plus<double>());
                    ++sample_counter;
                    if (sample_counter % Nc == 0) {
                        std::transform(cumul_sample.begin(), cumul_sample.end(),
                                        cumul_sample.begin(), [&](double const& a){return a/Nc;});

                        #pragma omp critical
                        {
                            problem.add_sample(cumul_sample, ppoint);
                        }
                        cumul_sample.assign(cumul_sample.size(), 0.);
                    }
                }
            }
        }
        else
            throw std::runtime_error("sample_config(): parameter sweep.Nc must be >= 1");
    }

    void reset_sweeps(bool skip_therm = false) override {
        Simulation::reset_sweeps(skip_therm);
        //i_sample = 0;
    }

protected:
    std::unique_ptr<config_policy_t> confpol;

private:
    using Simulation::parameters;
    using Simulation::random;

    size_t N_sample;
    size_t Nc;
    //size_t i_sample = 0;

    problem_t problem;
    svm::serialization::problem_serializer<svm::hdf5_tag, problem_t> prob_serializer;
};

}
