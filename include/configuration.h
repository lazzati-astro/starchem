#pragma once

#include <string>
#include <boost/program_options.hpp>

namespace options = boost::program_options;

struct configuration
{
    double temperature_0;
    double ode_dt_0;
    double ode_time_0;
    double ode_time_n;
    double ode_abs_err;
    double ode_rel_err;

    int io_disk_n_steps;
    int io_screen_n_steps;

    std::string network_file;
    std::string abundance_file;
    std::string environment_file;

    options::options_description  desc;
    options::variables_map        vm;

    configuration();
    void read_config ( const std::string &filename );
};


