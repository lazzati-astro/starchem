#pragma once

#include <bits/stdc++.h>
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

  std::string network_file;
  std::string abundance_file;
  std::string environment_file;
  std::vector<std::string> nucleation_files;

  options::options_description  desc;
  options::variables_map        vm;

  configuration();
  void read_config(const std::string& filename);
};


