#include <bits/stdc++.h>

#include <boost/numeric/odeint.hpp>
#include <boost/format.hpp>
#include <plog/Log.h>

#include "constants.h"
#include "spline.h"
#include "configuration.h"
#include "network.h"
#include "cell.h"

using namespace boost::numeric::odeint;
//using namespace boost::format;

cell::cell(network *n, uint32_t id, const spec_v& init_s, const cell_input& input_data) : net(n), cid(id)
{
  set_init_abundances(init_s, input_data);
  set_env_data(input_data);
  LOGD << "Created cell " << id;
}

void
cell::set_init_abundances(const spec_v& init_s, const cell_input& init_data)
{
  initial_abundances.resize(net->n_species);
  std::fill(initial_abundances.begin(), initial_abundances.end(), 0.0E0);

  LOGI << "Setting " << net->n_species << " initial species";
  
  for(auto i = 0; i < init_s.size(); i++)
  {
    auto idx = net->find_spec_idx(init_s[i]);
    if (idx != -1)
    {
      initial_abundances[idx] = init_data.initial_abundances[i];
      LOGD << "x[" << idx << "] (" << init_s[i] << ") = " << init_data.initial_abundances[i];
    }
    else
    {
      LOGD << "(" << init_s[i] << ") not found in network! it will be ignored";
    }
  }

}

void
cell::set_env_data(const cell_input& input_data)
{
  LOGI << "Setting environment data (temperatures, volumes)";
  env_times.assign(input_data.inp_times.begin(), input_data.inp_times.end());
  env_temperatures.assign(input_data.inp_temperatures.begin(), input_data.inp_temperatures.end());
  env_volumes.assign(input_data.inp_volumes.begin(), input_data.inp_volumes.end());

  LOGI << "Creating splines for interpolation of environment variables";
  env_temperature_spline.set_points(env_times, env_temperatures);
  env_volume_spline.set_points(env_times, env_volumes);
}

void 
cell::solve(const configuration& config)
{
  auto abs_err = config.ode_abs_err, rel_err = config.ode_rel_err;
  auto time_start = env_times.front(), time_end = env_times.back();
  auto dt0 = config.ode_dt_0;
 
  auto store_every=config.io_disk_n_steps, dump_every=config.io_screen_n_steps;

  LOGI << "Solver settings are ABS_ERR = " << abs_err << ", REL_ERR = " << rel_err;
  LOGI << "Beginning integration TIME_0 = " << time_start << " to TIME_N = " << time_end;

  typedef runge_kutta_cash_karp54< abundance_v > stepper_type;
  auto stepper = make_controlled(abs_err, rel_err, stepper_type());

  auto steps = integrate_adaptive(stepper, std::ref(*(this)), initial_abundances, time_start, time_end, dt0, cell_observer(solution_abundances, solution_times, store_every, dump_every)); 

  LOGI << "Integration finished with " << steps << " steps";

}

void
cell::print_abundances(const spec_v& following)
{
  std::string filename = "output" + std::to_string(cid) + ".dat";
  std::ofstream outf(filename);

  std::vector<int> following_idx;

  LOGI << "Outputting integration [" << cid << "] to file " << filename;

  outf << "integration of cell [" << cid << "]\n";
  outf << "time ";

  for (const auto &f : following)
  {
    auto idx = net->find_spec_idx(f);
    if (idx != -1)
    {
      following_idx.push_back(idx);
      outf << f << " ";
    }
  }
  outf << "\n";

  boost::format fmt("%1$5e ");
  for(auto i = 0; i < solution_times.size(); i++)
  {
    outf << fmt % solution_times[i];

    for(auto s_idx : following_idx)
      outf << fmt % solution_abundances[i][s_idx];

    outf << "\n";
  }

  LOGI << "File output complete";
  
} 

void
cell::operator() (const abundance_v &x, abundance_v &dxdt, const double t)
{
  using constants::kBeta;
  using constants::equPres;
  using constants::shapeFactor;

//  LOGD << "Beginning integration step (t = " << t << ")";
//  LOGI << "using " << net->reactions.size() << " reactions";

  double fi;
  double temperature;
  double volume, dvolume;

  temperature = env_temperature_spline(t);
  env_volume_spline.val_and_deriv(t, volume, dvolume); 

  std::fill(dxdt.begin(), dxdt.end(), 0.0);

  auto d_rho = -dvolume / volume;
  for(auto i = 0; i < dxdt.size(); i++)
  {
    dxdt[i] = d_rho * x[i];
  } 

  for(auto &reaction : net->reactions)
  {
    fi = 1.0;

    if (reaction.type == REACTION_TYPE_NUCLEATE)
    {

      bool nucleation_ok = true;
      for(const auto& r_idx : reaction.reacts_idx)
      {
        if (x[r_idx] < 1.0E-10)
        {
          nucleation_ok = false;
          break;
        }
      }

      if (!nucleation_ok) continue;
      auto pressure_equilibrium = equPres(reaction.alpha, reaction.beta, temperature);
      auto pressure = x[reaction.reacts_idx[0]] * kBeta(temperature);
      auto saturation = pressure / pressure_equilibrium;

      auto critical_size = net->nucl_rate_data[reaction.num].interpolate(1, temperature, saturation);
      if ( critical_size > 0.0 )
      {
        auto log_nucleation_rate = net->nucl_rate_data[reaction.num].interpolate(0, temperature, saturation);
        auto grains_nucleated = pow(10.0,log_nucleation_rate) * critical_size;

        for(const auto& r_idx : reaction.reacts_idx)
        {
          dxdt[r_idx] -= grains_nucleated;
        }
        for(const auto& p_idx : reaction.prods_idx)
        {
          dxdt[p_idx] += grains_nucleated;
        }
        
      }
    }
    else
    {
      for(const auto &r_idx : reaction.reacts_idx)
        fi *= x[r_idx];
      
      fi *= reaction.rate(temperature);
    
      for(const auto &r_idx : reaction.reacts_idx)
        dxdt[r_idx] -= fi;

      for(const auto &p_idx : reaction.prods_idx)
        dxdt[p_idx] += fi;
    } 
  } 
 
//  LOGD << "Integration step complete";
}

/*void
cell::jacobian(const abundance_v &x, jacobi_m &J, const double &t, abundance_v &dfdt)
{
  
} 
*/

