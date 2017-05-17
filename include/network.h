#pragma once

#include <bits/stdc++.h>

#define BOOST_SPIRIT_DEBUG
#include <boost/fusion/adapted.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>

namespace qi = boost::spirit::qi;
namespace phx = boost::phoenix;

typedef std::vector<std::string> spec_v;
typedef std::vector<size_t> specidx_v;

struct reaction
{
  spec_v reacts;
  spec_v prods;

  specidx_v reacts_idx;
  specidx_v prods_idx;

  double alpha;
  double beta;
  double gamma;

  int type;
  int num;

  double rate(double tgas);
};

typedef std::vector<reaction> reaction_v;

BOOST_FUSION_ADAPT_STRUCT(
          reaction, 
          (spec_v, reacts)
          (spec_v, prods)
          (double, alpha)
          (double, beta)
          (double, gamma)
          (int, type)
          (int, num)
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
        reaction_rule = spec_rule >> lit("->") >> spec_rule >> 
                  double_ >> double_ >> double_ >> 
                  int_ >> int_;

        network_rule = reaction_rule % +char_("\n") >> omit[*space] > eoi;

        BOOST_SPIRIT_DEBUG_NODES((network_rule)(reaction_rule)(spec_rule));
    }
    qi::rule<Iterator, reaction_v(), Skipper> network_rule;
    qi::rule<Iterator, reaction(), Skipper> reaction_rule;
    qi::rule<Iterator, spec_v(), Skipper> spec_rule;
};

struct network
{

  static const std::map<std::string, int> elements;
  
  reaction_v reactions;
  spec_v     species;

  int n_species;
  int n_reactions;
  
  int get_species_list();
  int map_species_to_reactions();
  size_t find_spec_idx(const std::string &spec);

  network();
  network(const std::string& chemfile);
  virtual ~network();

  int read_network(const std::string& chemfile);
};
