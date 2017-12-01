#include <vector>
#include <string>
#include <chrono>
#include <fstream>

#include <plog/Log.h>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

#include "starchem.h"

namespace options = boost::program_options;

StarChem::StarChem ( const std::string &config_file )
{
    sc_config.read_config ( config_file );
}

void
StarChem::load_network()
{
    LOGI << "load network";
    net.read_network ( sc_config.network_file );

    LOGI << "post process network";
    net.post_process();

    LOGI << "network loaded";
}

void
StarChem::load_initial_abundances()
{
    LOGI << "Loading initial abundances";
    std::ifstream abundance_file ( sc_config.abundance_file );

    std::string line_buffer;
    std::vector<std::string> line_tokens;

    if ( abundance_file.is_open() )
    {
        LOGI << "Reading element list from file headers ";
        // read the header line for column->species mapping
        std::getline ( abundance_file, line_buffer );
        boost::split ( line_tokens, line_buffer, boost::is_any_of ( " \t" ), boost::token_compress_on );
        initial_elements.assign ( line_tokens.begin() + 1, line_tokens.end() );

        for ( const auto &e : initial_elements )
            LOGD << "got element " << e;

        // now we read in the abundance values and cell numbers
        while ( std::getline ( abundance_file, line_buffer ) )
        {
            boost::split ( line_tokens, line_buffer, boost::is_any_of ( " \t" ), boost::token_compress_on );

            auto cell_id = boost::lexical_cast<int> ( line_tokens[0] );

            for ( auto it = line_tokens.begin() + 1; it != line_tokens.end(); ++it )
            {
                cell_inputs[cell_id].initial_abundances.push_back ( boost::lexical_cast<double> ( *it ) );
            }

            LOGI << "Read abundances for cell no. " << cell_id;
        }
    }
    else
    {
        LOGE << "Cannot open abundance file " << sc_config.abundance_file;
    }

}

void
StarChem::load_environment_data()
{
    LOGI << "Loading stellar environment";
    std::ifstream env_file ( sc_config.environment_file );

    std::string line_buffer;
    std::vector<std::string> line_tokens;

    double time;
    if ( env_file.is_open() )
    {
        while ( std::getline ( env_file, line_buffer ) )
        {
            boost::split ( line_tokens, line_buffer, boost::is_any_of ( "\t" ) );

            if ( line_tokens.size() == 1 )
            {
                time = boost::lexical_cast<double> ( line_tokens[0] );
            }
            else
            {
                auto cid = boost::lexical_cast<int> ( line_tokens[0] );
                auto temp = boost::lexical_cast<double> ( line_tokens[1] );
                auto vol = boost::lexical_cast<double> ( line_tokens[2] );

                cell_inputs[cid].inp_times.push_back ( time );
                cell_inputs[cid].inp_temperatures.push_back ( temp );
                cell_inputs[cid].inp_volumes.push_back ( vol );
            }
        }
    }
    else
    {
        LOGE << "Cannont open environment file " << sc_config.environment_file;
    }
}

void
StarChem::create_simulation_cells()
{
    LOGI << "Creating cells with input data";

    for ( const auto &ic : cell_inputs )
    {
        cells.emplace_back ( &net, &sc_config, ic.first, initial_elements, ic.second );
    }

    LOGI << cells.size() << " cells created.";
}

void
StarChem::run()
{

    LOGI << "Entering main integration loop";
//    size_t max_now = cells.size();
//#pragma omp parallel for
//    for ( auto c = 0; c < max_now; ++c )
//    {
    cells[0].solve ( following_species );
//    }
    LOGI << "Leaving main integration loop";
//	cells[0].print_abundances ( following_species );
}
