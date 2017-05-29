#include <vector>
#include <string>
#include <fstream>
#include <chrono>

#include <boost/numeric/odeint.hpp>
#include <boost/format.hpp>
#include <plog/Log.h>

#include "constants.h"
#include "spline.h"
#include "configuration.h"
#include "network.h"
#include "cell.h"
#include "cellobserver.h"

using namespace boost::numeric::odeint;
using namespace std::chrono;

cell::cell ( network *n, uint32_t id, const spec_v &init_s,
             const cell_input &input_data )
    : net ( n )
    , cid ( id )
{
    set_init_abundances ( init_s, input_data );
    set_env_data ( input_data );

    state_vars.parts.resize ( net->n_nucleation_reactions );

    LOGD << "Created cell " << id;
}

void
cell::set_init_abundances ( const spec_v &init_s, const cell_input &init_data )
{
    initial_abundances.resize ( net->n_species );
    std::fill ( initial_abundances.begin(), initial_abundances.end(), 0.0E0 );

    LOGI << "Setting " << net->n_species << " initial species";

    for ( auto i = 0; i < init_s.size(); ++i )
    {
        auto idx = net->get_species_index ( init_s[i] );
        if ( idx != -1 )
        {
            initial_abundances[idx] = init_data.initial_abundances[i];
            LOGD << "x[" << idx << "] (" << init_s[i]
                 << ") = " << init_data.initial_abundances[i];
        }
        else
        {
            LOGD << "(" << init_s[i] << ") not found in network! it will be ignored";
        }
    }
}

void
cell::set_env_data ( const cell_input &input_data )
{
    LOGI << "Setting environment data (temperatures, volumes)";
    env_times.assign ( input_data.inp_times.begin(), input_data.inp_times.end() );
    env_temperatures.assign ( input_data.inp_temperatures.begin(),
                              input_data.inp_temperatures.end() );
    env_volumes.assign ( input_data.inp_volumes.begin(),
                         input_data.inp_volumes.end() );

    LOGI << "Creating splines for interpolation of environment variables";
    env_temperature_spline.set_points ( env_times, env_temperatures );
    env_volume_spline.set_points ( env_times, env_volumes );
}

void
cell::solve ( const configuration &config, const spec_v& foll )
{
    auto abs_err = config.ode_abs_err, rel_err = config.ode_rel_err;
    auto time_start = env_times.front(), time_end = env_times.back();
    auto dt0 = config.ode_dt_0;

    auto store_every = config.io_disk_n_steps,
         dump_every = config.io_screen_n_steps;

    LOGI << "Solver settings are ABS_ERR = " << abs_err
         << ", REL_ERR = " << rel_err;
    LOGI << "Beginning integration TIME_0 = " << time_start
         << " to TIME_N = " << time_end;

    auto rkd = runge_kutta_dopri5<abundance_v> {};
    auto stepper = make_dense_output ( abs_err, rel_err, rkd );

    auto solve_steps = 0;

    current_abundances.assign ( initial_abundances.begin(),
                                initial_abundances.end() );

    stepper.initialize ( current_abundances, time_start, dt0 );
    calc_state_vars( current_abundances, time_start );

    CellObserver observer ( solution_abundances, solution_times, solution_vars, store_every,
                            dump_every , cid, foll, net );

    while ( ( stepper.current_time() < time_end ) )
    {

        auto step_start_time = high_resolution_clock::now();
        stepper.do_step ( std::ref ( * ( this ) ) );
        auto step_end_time = high_resolution_clock::now();



        observer ( stepper.current_state(), state_vars, stepper.current_time(),
                   stepper.current_time_step() );

        //    for(auto &interpolators : net->nucl_rate_data)
        //    {
        //
        //    }

        ++solve_steps;
    }

    LOGI << "Integration finished with " << solve_steps << " steps";
}

void cell::calc_state_vars ( const abundance_v &x, const double time )
{
    using constants::kBeta;
    using constants::equPres;
    using constants::shapeFactor;

    // get temperature and volume from the splines
    state_vars.temperature = env_temperature_spline ( time );
    env_volume_spline.val_and_deriv ( time, state_vars.volume, state_vars.dvolume );

    state_vars.drho = -state_vars.dvolume / state_vars.volume;

    for ( auto i = 0; i < net->n_nucleation_reactions; ++i )
    {
        auto reaction_idx = net->nucleation_reactions_idx[i];
        auto key_spec_idx = net->reactants_idx[reaction_idx][0];

        state_vars.parts[i].pressure = x[key_spec_idx] * kBeta ( state_vars.temperature );
        state_vars.parts[i].equilibrium_pressure
            = equPres ( net->reactions[reaction_idx].alpha, net->reactions[reaction_idx].beta, state_vars.temperature );

        state_vars.parts[i].saturation = state_vars.parts[i].pressure / state_vars.parts[i].equilibrium_pressure;

        state_vars.parts[i].critical_size =
            net->nucl_rate_data[net->reactions[reaction_idx].id].interpolate ( 1, state_vars.temperature,
                    state_vars.parts[i].saturation );
        if ( state_vars.parts[i].critical_size > 0.0 )
        {
            state_vars.parts[i].log_nucleation_rate =
                net->nucl_rate_data[net->reactions[reaction_idx].id].interpolate ( 0, state_vars.temperature,
                        state_vars.parts[i].saturation );
            state_vars.parts[i].nucleation_rate = pow ( 10.0, state_vars.parts[i].log_nucleation_rate );
            state_vars.parts[i].grains_nucleating = state_vars.parts[i].nucleation_rate * state_vars.parts[i].critical_size;

            state_vars.parts[i].is_nucleating = true;
        }
        else
        {
            state_vars.parts[i].log_nucleation_rate = 0.0;
            state_vars.parts[i].nucleation_rate = 0.0;
            state_vars.parts[i].grains_nucleating = 0.0;

            state_vars.parts[i].is_nucleating = false;
        }

    }
}

void cell::operator() ( const abundance_v &x, abundance_v &dxdt, const double t )
{
    //  LOGD << "Beginning integration step (t = " << t << ")";
    //  LOGI << "using " << net->reactions.size() << " reactions";

    double fi;

    std::fill ( dxdt.begin(), dxdt.end(), 0.0 );

    calc_state_vars ( x, t );

    for ( auto i = 0; i < dxdt.size(); ++i )
    {
        dxdt[i] = state_vars.drho * x[i];
    }

    for ( auto i = 0; i < net->n_nucleation_reactions; ++i )
    {
        auto reaction_idx = net->nucleation_reactions_idx[i];
        for ( const auto &r : net->reactants_idx[reaction_idx] )
            dxdt[r] -= state_vars.parts[i].grains_nucleating;
        for ( const auto &p : net->products_idx[reaction_idx] )
            dxdt[p] += state_vars.parts[i].grains_nucleating;

    }

    for ( auto i = 0; i < net->n_chemical_reactions; ++i )
    {
        auto reaction_idx = net->chemical_reactions_idx[i];
        fi = 1.0;
        for ( const auto &r : net->reactants_idx[reaction_idx] )
            fi *= x[r];

        fi *= net->reactions[reaction_idx].rate ( state_vars.temperature );

        for ( const auto &r : net->reactants_idx[reaction_idx] )
            dxdt[r] -= fi;

        for ( const auto &p : net->products_idx[reaction_idx] )
            dxdt[p] += fi;

    }

    //  LOGD << "Integration step complete";
}

/*void
cell::jacobian(const abundance_v &x, jacobi_m &J, const double &t, abundance_v
&dfdt)
{

}
*/
