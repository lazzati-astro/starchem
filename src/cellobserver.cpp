#include <vector>
#include <string>

#include "network.h"

#include "cell.h"
#include "cellobserver.h"

CellObserver::CellObserver (  std::vector <abundance_v> &states,
                              std::vector<double> &times,
                              std::vector<cell_state> &statevars,
                              uint32_t store_every_n,
                              uint32_t dump_every_n,
                              size_t cid ,
                              const spec_v& foll,
                              const network* net )
    : m_states ( states ),
      m_times ( times ),
      m_vars ( statevars ),
      m_nstore ( store_every_n ),
      m_ndump ( dump_every_n ),
      n_called ( 0 ),
      following ( foll )
{
    ofname = "output/" + std::to_string ( cid ) + ".dat";
    ofs.open(ofname);

//    ofs << "integration of cell [" << cid << "]\n";
    ofs << "time temperature ";

    for ( const auto &f : following )
    {
        auto idx = net->get_species_index ( f );
        if ( idx != -1 )
        {
            following_idx.push_back ( idx );
            ofs << f << " ";
        }
    }
    ofs << "\n";
    ofs.close();
}

void CellObserver::dump_abundances ( )
{

    ofs.open( ofname, std::ofstream::app );
    boost::format fmt ( "%1$5e " );

    auto idx = 0;
    for ( auto idx = 0; idx < m_times.size(); ++idx )
    {
        ofs << fmt % m_times[idx];
        ofs << fmt % m_vars[idx].temperature;

        for ( auto s_idx : following_idx )
            ofs << fmt % (m_states[idx][s_idx] * m_vars[idx].volume);

        ofs << "\n";
    }

    m_states.clear();
    m_times.clear();
    m_vars.clear();
    ofs.close();

//	LOGI << "File output complete";
}


void CellObserver::operator() ( const abundance_v &x,
                                const cell_state &s,
                                const double t,
                                const double dt  )
{
    ++n_called;
    if ( n_called % m_nstore == 0 )
    {
        m_states.push_back ( x );
        m_times.push_back ( t );
        m_vars.push_back ( s );
    }

    if ( n_called % m_ndump == 0 )
    {
        std::cout << "t = " << t << " dt = " << dt << " T = " << s.temperature <<  " S = " << s.parts[1].saturation << " J = " << s.parts[1].nucleation_rate << std::endl;
        dump_abundances();
    }
}


