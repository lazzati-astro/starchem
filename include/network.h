#pragma once

#include <bits/stdc++.h>

#define BOOST_SPIRIT_DEBUG
#include <boost/fusion/adapted.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>

#include "reaction.h"
#include "interp2D.h"

namespace qi = boost::spirit::qi;
namespace phx = boost::phoenix;

typedef std::vector<reaction> reaction_v;

BOOST_FUSION_ADAPT_STRUCT(
          reaction, 
          (spec_v, reacts)
          (spec_v, prods)
          (double, alpha)
          (double, beta)
          (double, gamma)
          (int, type)
          (int, id)
          (std::string, extra)
)


template<typename Iterator, typename Skipper>
struct network_parser : qi::grammar<Iterator, reaction_v(), Skipper>
{
    network_parser() : network_parser::base_type(network_rule)
    {
        using qi::lexeme;
        using qi::lit;
        using qi::graph;
        using qi::double_;
        using qi::int_;
        using qi::char_;
        using qi::omit;
        using qi::space;
        using qi::eoi;

        spec_rule = lexeme [ +graph ] % '+';
        extra_rule = +(char_ - '\n');
        reaction_rule = spec_rule >> lit("->") >> spec_rule >> 
                  double_ >> double_ >> double_ >> 
                  int_ >> int_ >> *(extra_rule);

        network_rule = reaction_rule % +char_("\n") >> omit[*space] > eoi;

        BOOST_SPIRIT_DEBUG_NODES((network_rule)(reaction_rule)(spec_rule)(extra_rule));
    }
    qi::rule<Iterator, reaction_v(), Skipper> network_rule;
    qi::rule<Iterator, reaction(), Skipper> reaction_rule;
    qi::rule<Iterator, spec_v(), Skipper> spec_rule;
    qi::rule<Iterator, std::string(), Skipper> extra_rule;
};



struct network
{

  std::vector<reaction>       reactions;
  std::vector<std::string>    species;

  std::vector<std::vector<int>> reactants_idx;
  std::vector<std::vector<int>> products_idx;

  std::map<int, interp2D> nucl_rate_data;

  int         n_species;
  int         n_reactions;
  
  int         get_species_list();
  int         map_species_to_reactions();

  size_t      get_species_index(const std::string &spec);

  int         read_network(const std::string& chemfile);
  void        post_process();

  network();
  virtual ~network();

};
