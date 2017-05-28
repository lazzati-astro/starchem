#pragma once

#include <vector>
#include <string>
#include <map>
#include <chrono>

#include <boost/serialization/array_wrapper.hpp>
#include <boost/numeric/odeint.hpp>

#include "configuration.h"
#include "spline.h"
#include "network.h"

const double CELL_MINIMUM_ABUNDANCE = 1.0E-10;
//typedef std::vector<double> abundance_v;
typedef std::vector<double> abundance_v;

//typedef boost::numeric::ublas::vector< double > abundance_v;
//typedef boost::numeric::ublas::matrix< double > jacobi_m;

struct cell_observer {
    std::vector<abundance_v> &m_states;
    std::vector<double> &m_times;

    uint32_t m_nstore, m_ndump;
    uint32_t n_called;

    cell_observer( std::vector <abundance_v> &states, std::vector<double> &times, uint32_t store_every_n, uint32_t dump_every_n)
        : m_states(states), m_times(times), m_nstore(store_every_n), m_ndump(dump_every_n), n_called(0) {}

    void operator() (const abundance_v &x, const double t, const double dt) {
        n_called++;
        if (n_called % m_nstore == 0) {
            m_states.push_back(x);
            m_times.push_back(t);
        }

        if (n_called % m_ndump == 0) {
            std::cout << "t = " << t << " dt = " << dt << std::endl;
        }
    }
};

struct cell_input {
    std::vector<double> initial_abundances;
    std::vector<double> inp_times;
    std::vector<double> inp_temperatures;
    std::vector<double> inp_volumes;
};

struct cell_state {
  double temperature;
  double volume;
  std::vector<double> pressure;
  std::vector<double> saturation;

  cell_state() : 
    temperature(0),
    pressure(0),
    saturation(0),
    volume(0)
  {}

  cell_state(double t, double p, double s, double v):
      temperature(t),
      pressure(p),
      saturation(s),
      volume(v)
      {}
};

class cell {

    network *net;

    abundance_v initial_abundances;
    abundance_v current_abundances;

    std::vector<bool> reaction_switch;

    std::vector<abundance_v> solution_abundances;
    std::vector<double> solution_times;

    std::vector<double> env_times;
    std::vector<double> env_temperatures;
    std::vector<double> env_volumes;

    tk::spline env_temperature_spline;
    tk::spline env_volume_spline;

    cell_state state_vars;

    uint32_t cid;

    void set_init_abundances(const spec_v &init_s, const cell_input &init_data);
    void set_env_data(const cell_input &input_data);

  public:

    cell(network *n, uint32_t id, const spec_v &init_s, const cell_input &input_data);

    virtual ~cell() { };

    void solve(const configuration &c);

    void print_abundances(const spec_v &following);

    void operator() (const abundance_v &x, abundance_v &dxdt, const double t);

    int get_id() {
        return cid;
    }


//  void jacobian(const abundance_v &x, jacobi_m &J, const double &t, abundance_v &dfdt);

};
