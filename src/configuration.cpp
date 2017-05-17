#include <bits/stdc++.h>
#include <plog/Log.h>
#include <boost/program_options.hpp>

#include "configuration.h"



configuration::configuration(): desc("configuration")
{
  desc.add_options() ("temperature_0", options::value<double>(&temperature_0), "standalone temperature (for testing)");
  desc.add_options() ("ode_dt_0", options::value<double>(&ode_dt_0), "initial time step");
  desc.add_options() ("ode_time_0", options::value<double>(&ode_time_0), "integration start time");
  desc.add_options() ("ode_time_n", options::value<double>(&ode_time_n), "integration end time");
  desc.add_options() ("ode_abs_err", options::value<double>(&ode_abs_err), "solver absolute error criteria");
  desc.add_options() ("ode_rel_err", options::value<double>(&ode_rel_err), "solver relative error criteria");
  desc.add_options() ("network_file", options::value<std::string>(&network_file), "file with network");
  desc.add_options() ("abundance_file", options::value<std::string>(&abundance_file), "file with inital abundances");
  desc.add_options() ("environment_file", options::value<std::string>(&environment_file), "file with stellar environment variables");
}

void
configuration::read_config(const std::string& config_filename)
{
  LOGI << "Reading configuration file: " << config_filename;
  std::ifstream config_file(config_filename.c_str());

  vm = options::variables_map();

  options::store(options::parse_config_file(config_file, desc), vm);
  options::notify(vm);

  LOGI << "configuration read";
  LOGD << "configuration: temperature_0 = " << temperature_0;
  LOGD << "configuration: ode_dt_0 = " << ode_dt_0;
  LOGD << "configuration: ode_time_0 = " << ode_time_0;
  LOGD << "configuration: ode_time_n = " << ode_time_n;
  LOGD << "configuration: ode_abs_err = " << ode_abs_err;
  LOGD << "configuration: ode_rel_err = " << ode_rel_err;
  LOGD << "configuration: network_file = " << network_file;
  LOGD << "configuration: abundance_file = " << abundance_file;
  LOGD << "configuration: environment_file = " << environment_file;
}


