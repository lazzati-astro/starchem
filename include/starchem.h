#pragma once

#include "configuration.h"
#include "network.h"
#include "cell.h"

namespace options = boost::program_options;

class star_chem
{

  configuration                 sc_config;
  std::vector<std::string>      initial_elements;
  
  std::map<uint32_t, cell_input> cell_inputs;

  std::vector<cell>             cells;
  network                       net;

  std::vector<std::string>      following_species;

 public:
  star_chem(const std::string& config_filename);
  virtual ~star_chem() {}

  void load_network();
  void load_initial_abundances();
  void load_environment_data();

  void create_simulation_cells();

  void set_following(const spec_v& following_list) { following_species.assign(following_list.begin(), following_list.end()); }

  void run(); 
  void cleanup() {}
};
