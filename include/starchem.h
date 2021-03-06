#pragma once

#include <vector>
#include <map>
#include <string>

#include "configuration.h"
#include "network.h"
#include "cell.h"

namespace options = boost::program_options;

class StarChem
{

    configuration                 sc_config;
    std::vector<std::string>      initial_elements;

    std::map<uint32_t, cell_input> cell_inputs;

    std::vector<cell>             cells;
    network                       net;

    std::vector<std::string>      following_species;

public:
    StarChem ( const std::string &config_filename );
    virtual ~StarChem() {}

    void load_network();
    void load_initial_abundances();
    void load_environment_data();

    void create_simulation_cells();

    void set_following ( const spec_v &following_list )
    {
        following_species.assign ( following_list.begin(), following_list.end() );
    }

    void run();
    void cleanup() {}
};
