#pragma once

#include "cell.h"
#include "network.h"
#include <vector>

class CellObserver
{
private:
    std::vector<abundance_v> &m_states;
    std::vector<double> &m_times;
    std::vector<cell_state> &m_vars;

    uint32_t m_nstore, m_ndump;
    uint32_t n_called;

    std::vector<int> following_idx;
    std::vector<std::string> following;
    std::string   ofname;
    std::ofstream ofs;

public:
    CellObserver ( std::vector <abundance_v> &states, std::vector<double> &times, std::vector<cell_state> &vars, uint32_t store_every_n, uint32_t dump_every_n, size_t cid , const spec_v& foll, const network* net);
    virtual ~CellObserver() {}

    void operator() ( const abundance_v &x, const cell_state &s, const double t, const double dt );
    void dump_abundances();

};

