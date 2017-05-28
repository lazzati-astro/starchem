#pragma once

#include <vector>
#include <string>

typedef std::vector<std::string> spec_v;

const int REACTION_TYPE_NUCLEATE = 101;

struct reaction {
    spec_v reacts;
    spec_v prods;

    double alpha;
    double beta;
    double gamma;

    int type;
    int id;

    std::string extra;



    reaction() : id(-1), type(-1) {}
    ~reaction() {}

    double rate(double tgas);

};


