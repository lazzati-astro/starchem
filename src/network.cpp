#include <bits/stdc++.h>
#include <plog/Log.h>

#include <boost/spirit/include/qi.hpp>
#include "network.h"


/* a map for element bookkeeping. I can't find a use for it yet,
 * but worth keeping around */
const std::map<std::string, int> elements =
  {
    {"H",1},
    {"D",2},
    {"He",2},
    {"C",12},
    {"-13-C",13},
    {"N",14},
    {"-15-N",15},
    {"O",16},
    {"-17-O",17},
    {"-18-O",18},
    {"F",19},
    {"Ne",20},
    {"Na",23},
    {"Mg",24},
    {"Si",28},
    {"-29-Si",29},
    {"S",30},
    {"P",31},
    {"Cl",35},
    {"Fe",56},
    {"grain",0},
    {"X",1},
    {"Y",1}
  };

double
reaction::rate(double tgas)
{
  double k = 0.0;
  switch(type)
  {
    case 2:
    case 3:
    case 4:
    case 5:
    case 6:
    case 7:
    case 8:
    case 9:
    case 10:
    case 11:
    case 12:
    case 14:
      k = alpha * pow( tgas / 300.0, beta ) * exp( -gamma / tgas);
      break;
    default:
      LOGE << "UNKNOWN TYPE (" << type <<")";
  }
  return k;
}

network::network() : 
  n_reactions(0), 
  n_species(0)
{
}

network::~network()
{

}

//is this necessary?
network::network(const std::string& chemfile)
{
  read_network(chemfile);
}

/* read the network file and store the chemical network 
 * the details of input parsing are in the header file.
 */
int
network::read_network(const std::string& chemfile)
{
  LOGI << "Reading network from [" << chemfile << "]";

  std::ifstream ifs(chemfile.c_str(), std::ifstream::in);
  ifs.unsetf(std::ios::skipws);
  boost::spirit::istream_iterator f(ifs), l;

  LOGI << "Network file [" << chemfile << "] opened";

  network_parser<boost::spirit::istream_iterator, qi::blank_type> p;

  LOGI << "parser applied";
  
  bool ok = qi::phrase_parse(f, l, p, qi::blank, reactions);

  if (ok)
  {
    LOGI << "parser ok";
  }
  else
  {
    LOGI << "parser fail!!";
  }

  LOGI << "geting species from network...";
  get_species_list();

  LOGI << "mapping species index to reactions...";
  map_species_to_reactions();  

  n_reactions = reactions.size();
  n_species = species.size();

  LOGI << "network loaded from [" << chemfile << "]";
  LOGI << "got " << n_species << " species in " << n_reactions << " reactions";
}

size_t
network::find_spec_idx(const std::string &spec)
{
  auto it = std::find(species.begin(), species.end(), spec);
  auto idx = -1;
  if (it == species.end())
    LOGE << "ERROR: cannot find species(" << spec <<")";
  else
    idx = std::distance(species.begin(), it);
  return idx;
}

int
network::map_species_to_reactions()
{
  
  for(auto &r : reactions)
  {
    for(auto &reactant : r.reacts)
    { 
      r.reacts_idx.emplace_back(find_spec_idx(reactant));
    }

    for(auto &product : r.prods)
    {
      r.prods_idx.emplace_back(find_spec_idx(product));
    }

  }
}

int
network::get_species_list()
{
  std::unordered_set<std::string> spec_set;

  species.clear();  

  for(const auto &r : reactions)
  {
    spec_set.insert(r.reacts.begin(), r.reacts.end());
    spec_set.insert(r.prods.begin(), r.prods.end());
  } 
  
  species.insert(species.begin(), spec_set.begin(), spec_set.end());

  auto s_id = 0;
  for(const auto &s : species)
  {
    LOGI << "added..." << s << " to species list";
  }
}
