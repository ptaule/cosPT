/*
   io.hpp

   Created by Petter Taule on 03.09.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#ifndef IO_HPP
#define IO_HPP

#include <string>
#include <vector>

void read_delimited_file(
        const std::string& filename,
        std::vector<std::vector<double>>& data
        );

class Config;
void write_results(
        const Config& cfg,
        const std::vector<double>& tree_level_result,
        const std::vector<double>& loop_result,
        const std::vector<double>& errors
        );

#endif /* ifndef IO_HPP */
