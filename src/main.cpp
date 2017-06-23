#include <vector>
#include <string>

#include <plog/Log.h>
#include "starchem.h"

int main ( int argc, char *argv[] )
{

    plog::init ( plog::verbose, "log.txt", 0, 1 );

    LOGD << "starchem has started";

    std::string cfn ( "config.ini" );

    StarChem sc ( cfn );

    sc.load_network();
    sc.load_initial_abundances();
    sc.load_environment_data();
    sc.create_simulation_cells();

    spec_v following = {"C", "O", "Mg", "Si", "CO", "C2", "O2","c_grain", "si_grain"};
    sc.set_following ( following );

    sc.run();

    return 0;
}
