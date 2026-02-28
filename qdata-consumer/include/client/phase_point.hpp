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

#include <string>

#include <alps/params.hpp>

#include <tksvm/phase_space/point/common.hpp>


namespace tksvm {
namespace phase_space {
namespace point {

    struct phase_point {
        //client-specific: dimension of the phase point
        static const size_t label_dim = 1;
        using iterator = double *;
        using const_iterator = double const *;

        static void define_parameters(alps::params & params, std::string prefix="") {
                params
                    //temperature must be defined!
                    .define<double>(prefix + "T", 0., "temperature");
                    //client-specific: name and description of other phase point components 
        }

        static bool supplied(alps::params const& params, std::string prefix="") {
                return 
                    //client-specific: check if all phase point components are supplied
                    params.supplied(prefix + "T");
        }

        phase_point() {
            for(size_t i = 0; i < label_dim; ++i) p[i] = -1.;
        }

        template <class Iterator>
        phase_point(Iterator begin) {
            for(size_t i = 0; i < label_dim; ++i)
                p[i] = *(begin+i);
        }

        phase_point(alps::params const& params, std::string prefix="") {
            //client-specific: construct phase point from alps parameters
                p[0] = params[prefix + "T"].as<double>();
        }

        //accessor for temperature must be defined (it is used for parallel tempering)
        double const& temperature() const {return p[0];}
        double& temperature() {return p[0];}
        //client-specific: OPTIONALLY define other accessors

        const_iterator begin() const { return p; }
        iterator begin() { return p; }
        const_iterator end() const { return p + label_dim; }
        iterator end() { return p + label_dim; }

        double p[label_dim];
    };

}
}
}

