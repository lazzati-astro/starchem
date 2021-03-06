#pragma once

#include <vector>
#include <string>
#include <map>
#include <chrono>
#include <fstream>

#include <boost/serialization/array_wrapper.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/format.hpp>

#include "configuration.h"
#include "spline.h"
#include "network.h"

const size_t CELL_MAX_STEPS = 10000000000;
const size_t CELL_MAXIMUM_STEPPER_RESETS    = 16;
const double CELL_MINIMUM_ABUNDANCE = 1.0E-7;
typedef std::vector<double> abundance_v;
typedef std::vector<double> abundance_v;

//typedef boost::numeric::ublas::vector< double > abundance_v;
//typedef boost::numeric::ublas::matrix< double > jacobi_m;

struct cell_input
{
    std::vector<double> initial_abundances;
    std::vector<double> inp_times;
    std::vector<double> inp_temperatures;
    std::vector<double> inp_volumes;
};

struct cell_partial
{
    double reaction_idx;

    double pressure;
    double equilibrium_pressure;
    double saturation;

    bool  is_nucleating;
    double critical_size;
    double log_nucleation_rate;
    double nucleation_rate;

    double grains_nucleating;

};

struct cell_state
{
    double temperature;
    double volume;
    double dvolume;
    double drho;

    std::vector<cell_partial> parts;

};

class cell
{

    /*
     * these should be const but I'm lazy
     * and haven't implemented a "read-only"
     * interface yet.
     */
    network*        net;
    configuration*  config;

    abundance_v initial_abundances;
    abundance_v current_abundances;

    std::vector<abundance_v> solution_abundances;
    std::vector<double> solution_times;
    std::vector<cell_state> solution_vars;
    std::vector<bool> reaction_switch;

    std::vector<double> env_times;
    std::vector<double> env_temperatures;
    std::vector<double> env_volumes;

    tk::spline env_temperature_spline;
    tk::spline env_volume_spline;

    cell_state state_vars;

    bool integration_abandoned;

    uint32_t cid;

    void check_reactions ( const abundance_v &x );
    bool check_solution ( const abundance_v &x );
    void calc_state_vars ( const abundance_v &x, const double time );

    void set_init_abundances ( const spec_v &init_s, const cell_input &init_data );
    void set_env_data ( const cell_input &input_data );

public:

    cell ( network* n, configuration* con, uint32_t id, const spec_v &init_s, const cell_input &input_data );

    virtual ~cell() { };

    void solve ( const spec_v& foll );

    void operator() ( const abundance_v &x, abundance_v &dxdt, const double t );

    int get_id()
    {
        return cid;
    }

    //void jacobian(const abundance_v &x, jacobi_m &J, const double &t, abundance_v &dfdt);

};
