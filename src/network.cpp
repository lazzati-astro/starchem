#include <bits/stdc++.h>
#include <plog/Log.h>

#include <boost/spirit/include/qi.hpp>
#include "network.h"

network::network() : 
  n_reactions(0), 
  n_species(0)
{
}

network::~network()
{

}

//is this necessary?
//network::network(const std::string& chemfile)
//{
//  read_network(chemfile);
//}

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
    LOGE << "parser fail!!";
    exit(1);
  }
  
  if(f!=l) LOGI << "remaining unparsed: '" << std::string(f,l) << "'";

  LOGI << "getting species from network...";
  get_species_list();

  n_reactions = reactions.size();
  n_species = species.size();

  
  //LOGI << "mapping species index to reactions...";
  map_species_to_reactions();  


  LOGI << "network loaded from [" << chemfile << "]";
  LOGI << "got " << n_species << " species in " << n_reactions << " reactions";
}

/*
 * do any updates on read network
 * used to:
 * - read nucleation interpolators for dust grains
 *
 */
void
network::post_process()
{
  for(const auto& reaction : reactions)
  {
    if (reaction.type == REACTION_TYPE_NUCLEATE)
    {
      LOGI << "reaction " << reaction.id << " is dust nucleation";
      LOGI << "reading " << reaction.extra << " for nucleation data";

      std::string nucleation_file(reaction.extra);
      nucl_rate_data.insert(std::make_pair(reaction.id, interp2D(nucleation_file)));

      LOGI << "nucleation data loaded";
    }
  }
}

size_t
network::get_species_index(const std::string &spec)
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
  reactants_idx.resize(n_reactions);
  products_idx.resize(n_reactions);
  for(auto i = 0; i < n_reactions; i++)
  {
    for(const auto &reactant : reactions[i].reacts)
      reactants_idx[i].push_back(get_species_index(reactant));

    for(const auto &product : reactions[i].prods)
      products_idx[i].push_back(get_species_index(product));
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

