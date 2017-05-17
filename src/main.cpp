#include <bits/stdc++.h>
#include <plog/Log.h>
#include <boost/algorithm/string.hpp>

#include "starchem.h"
#include "interp2D.h"
#include "cell.h"
#include "network.h"

int main(int argc, char* argv[])
{

  plog::init(plog::debug, "log.txt");

  LOGD << "starchem has started";

  std::string cfn("config.ini");

  star_chem sc(cfn);

  sc.load_network();
  sc.load_initial_abundances();
  sc.load_environment_data();
  sc.create_simulation_cells();

  spec_v following = {"C","O","CO","C2","C3","CO2"};
  sc.set_following(following);

  
  sc.run();
//  c.set_init_abundances(elems, iabuns[0].);
//  c.solve();
//  c.print_abundances(following);

  return 0;
}
