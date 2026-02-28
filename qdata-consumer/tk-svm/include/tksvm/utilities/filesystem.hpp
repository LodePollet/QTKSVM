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
#include <regex>
#include <string>
#include <fstream>


namespace tksvm {

inline std::string replace_extension(std::string const& old_path,
                                     std::string const& new_ext)
{
	std::regex re{"(.+?)(\\.ini|\\.h5|\\.txt|\\.ppm|\\.out|\\.clone|\\.test)*"};
	std::smatch match;
	if (std::regex_match(old_path, match, re)) {
		return match[1].str() + new_ext;
	} else {
		return old_path + new_ext;
	}
}

template <typename Action>
void parse_merge_string(std::istream& merge_is, Action&& action)
{
    for (std::string mline; std::getline(merge_is, mline, '\n');) {
        std::stringstream iss(mline);
        for (std::string name; std::getline(iss, name, ':');) {
            std::smatch sm;
            if ( std::regex_match(name, sm, std::regex{"(\\[)(.*)(\\])"})) {
                std::ifstream merge_file{sm[2]};
                if (merge_file.is_open()) {
                    parse_merge_string(merge_file, action);
                }
                else
                    std::cerr << "unable to open merge-file: " << name << "\tskipping...\n";
            }
            else
                action(name);
        }
    }
}

}
