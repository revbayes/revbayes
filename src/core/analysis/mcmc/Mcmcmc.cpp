#include <iomanip>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstddef>
#include <functional>
#include <string>

#include "DagNode.h"
#include "MetropolisHastingsMove.h"
#include "Mcmcmc.h"
#include "Proposal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbException.h"
#include "RbMathLogic.h"
#include "Mcmc.h"
#include "Model.h"
#include "Monitor.h"
#include "MonteCarloAnalysisOptions.h"
#include "MonteCarloSampler.h"
#include "Move.h"
#include "RbFileManager.h"
#include "RbIterator.h"
#include "RbIteratorImpl.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "StringUtilities.h"

#ifdef RB_MPI
#include <mpi.h>
#endif

/*
 * Round trip statistics are inspired by the paper https://arxiv.org/pdf/cond-mat/0602085.pdf
 *
 * 1. The first idea is that simply having good acceptance probabilities for swapping adjacent temperatures is not good enough.
 *    We need chains to make a round-trip from cold -> hot -> cold.
 *    But it is possible to have a low rate of making round trips even if we have success swapping between adjacent temperatures.
 *
 * 2. The second idea is based on computing the fraction of the time that a chain at temperature j has most recently been at
         the coldest chain versus the hottest chain.
 *    The claim is that the fractions should be uniformly spaced.
 *    For example, if we have 5 chains, then the fraction should be 100% for the cold chain, and then 75%, 50%, 25% and 0%.
 *    In some cases, geometric heating produces very low round-trip rates, and this is because the fractions are NOT uniformly spaced.
 *
 * 3. The third idea is that one scenario where geometric heating works badly is when there is a "phase transition" at temperature T.
 *    This means that the posterior distribution on one side of T looks very different than the posterior distribution on
 *        the other side of T.
 *    This can happen when there are a very large number of low-probability states.
 *    As the temperature increases, the low-probability states increase in probability, and eventually overwhelm the high-probability
 *        states because there are so many of the lower-probability states.
 *    In such a case, we need closely-spaced temperatures near T, because the posterior changes rapidly as T changes.
 *    However, we can tolerate wider spacing between temperatures further away from T.
 *
 * BDR: I worry that diffusion behavior limits the number of round-trips as the number of chains increases.
 *      By diffusion behavior, I mean that a chain does not always move towards the coldest or hottest chain, but drifts randomly,
          sometimes stepping away.
 *      Does a chain take O(N^2) steps to move from cold -> hot, even if the temperatures are well placed?
 * 
 */

using namespace RevBayesCore;

using std::string;

Mcmcmc::Mcmcmc(const Model& m, const RbVector<Move> &mv, const RbVector<Monitor> &mn, std::string sT, size_t nc, size_t si, double dt, size_t ntries, bool th, double tht, std::string sm, std::string smo) : MonteCarloSampler( ),
    num_chains(nc),
    schedule_type(sT),
    current_generation(0),
    burnin_generation(0),
    generation(0),
    swap_interval(si),
    swap_interval2(0),
    active_chain_index( 0 ),
    delta( dt ),
    tune_heat(th),
    tune_heat_target(tht),
    useNeighborSwapping(true),
    useRandomSwapping(false),
    swap_mode(smo)
{
    
    // initialize container sizes
    chains = std::vector<Mcmc*>(num_chains, NULL);
    chain_values.resize(num_chains, 0.0);
    chain_heats.resize(num_chains, 0.0);
    chain_prev_boundary.resize(num_chains, boundary::intermediate);
    chain_half_trips.resize(num_chains, 0);
    heat_visitors.resize(num_chains, {0,0});
    pid_per_chain.resize(num_chains, 0);
    heat_ranks.resize(num_chains, 0);
    heat_temps.resize(num_chains, 0.0);
    
    num_attempted_swaps = std::vector< std::vector<size_t> > (num_chains, std::vector<size_t> (num_chains, 0));
    num_accepted_swaps = std::vector< std::vector<size_t> > (num_chains, std::vector<size_t> (num_chains, 0));
    
    chain_moves_tuningInfo = std::vector< std::vector<Mcmc::tuningInfo> > (num_chains);
    
    if (sm == "neighbor")
    {
        useNeighborSwapping = true;
        useRandomSwapping = false;
    }
    else if (sm == "random")
    {
        useNeighborSwapping = false;
        useRandomSwapping = true;
    }
    else if (sm == "both")
    {
        useNeighborSwapping = true;
        useRandomSwapping = true;
    }

    
    // assign chains to processors, instantiate Mcmc objects
    base_chain = new Mcmc(m, mv, mn, ntries);
    base_chain->setActivePID(this->pid, 1);
    
    // initialize the individual chains
    initializeChains();
}


Mcmcmc::Mcmcmc(const Mcmcmc &m) : MonteCarloSampler(m)
{
    
    delta               = m.delta;
    num_chains          = m.num_chains;
    heat_ranks          = m.heat_ranks;
    heat_temps          = m.heat_temps;
    swap_interval       = m.swap_interval;
    swap_interval2      = m.swap_interval2;
    swap_mode           = m.swap_mode;
    tune_heat           = m.tune_heat;
    tune_heat_target    = m.tune_heat_target;
    useNeighborSwapping = m.useNeighborSwapping;
    useRandomSwapping   = m.useRandomSwapping;
    
    active_chain_index  = m.active_chain_index;
    schedule_type       = m.schedule_type;
    pid_per_chain       = m.pid_per_chain;
    
    num_attempted_swaps = m.num_attempted_swaps;
    num_accepted_swaps  = m.num_accepted_swaps;
    generation          = m.generation;
    
    
    chains.clear();
    chains.resize(num_chains, NULL);
    for (size_t i = 0; i < num_chains; ++i)
    {
        if ( m.chains[i] != NULL)
        {
            chains[i]   = m.chains[i]->clone();
        }
        
    }
    
    chain_values            = m.chain_values;
    chain_heats             = m.chain_heats;
    chain_prev_boundary     = m.chain_prev_boundary;
    chain_half_trips        = m.chain_half_trips;
    heat_visitors           = m.heat_visitors;
    chain_moves_tuningInfo  = m.chain_moves_tuningInfo;
    
    burnin_generation       = m.burnin_generation;
    current_generation      = m.current_generation;
    base_chain              = m.base_chain->clone();
    
}

Mcmcmc::~Mcmcmc(void)
{
    for (size_t i = 0; i < chains.size(); ++i)
    {
        if (chains[i] != NULL)
        {
            delete chains[i];
        }
    }
    chains.clear();
    delete base_chain;
}


void Mcmcmc::addFileMonitorExtension(const std::string &s, bool dir)
{
    
    for (size_t i = 0; i < num_chains; ++i)
    {
        if ( chains[i] != NULL )
        {
            chains[i]->addFileMonitorExtension(s, dir);
        }
    }
    
}


void Mcmcmc::addMonitor(const Monitor &m)
{
    
    for (size_t i = 0; i < num_chains; ++i)
    {
        if ( chains[i] != NULL )
        {
            chains[i]->addMonitor( m );
        }
    }
    
}

size_t Mcmcmc::chainForHeatIndex(size_t i) const
{
    return std::find(chain_heats.begin(), chain_heats.end(), heatForIndex(i)) - chain_heats.begin();
/*
    // 1. Start with [0 .. num_chains-1]
    std::vector<size_t> tmp_heat_ranks;
    for(size_t i=0;i<num_chains;i++)
        tmp_heat_ranks.push_back(i);

    // 2. Chains with greater heat come first.
    sort(tmp_heat_ranks.begin(),
         tmp_heat_ranks.end(),
         [&](size_t j, size_t k) { return chain_heats[j] > chain_heats[k];} );

    // 3. Get the chain index for heat index i;
    return tmp_heat_ranks[i];
*/
}

void Mcmcmc::checkpoint( void ) const
{
    
    for (size_t i = 0; i < num_chains; ++i)
    {
        
        if ( chains[ chainForHeatIndex(i) ] != NULL )
        {
            // get the preliminary checkpoint file
            path f = chains[ chainForHeatIndex(i) ]->getCheckpointFile();
            
            std::string path_string = f.string();
            std::string test_string = "_chain_";
            path chain_file_name;
            
            // if the preliminary name does not contain "_chain_", append:
            size_t pos = path_string.find(test_string);
            if (pos == std::string::npos)
            {
                chain_file_name = appendToStem(f, test_string + std::to_string(i) );
            }
            // if it does contain it, replace the index:
            else
            {
                // if we have more than 9 chains, the number of characters in the index will be variable,
                // so we have to determine it first
                std::string expr = "_chain_*.";
                size_t idx_pos = expr.find("*");
                size_t tmp0 = path_string.find( expr.substr(0, idx_pos) );
                size_t tmp1 = path_string.find( expr.substr(idx_pos + 1) );
                size_t idx_length = tmp1 - tmp0 - idx_pos;
                
                path_string.replace(pos + test_string.length(), idx_length, std::to_string(i));
                chain_file_name = path_string;
            }
            
            chains[ chainForHeatIndex(i) ]->setCheckpointFile( chain_file_name );
            chains[ chainForHeatIndex(i) ]->checkpoint();
        }
        
    }
    
}


Mcmcmc* Mcmcmc::clone(void) const
{
    return new Mcmcmc(*this);
}


double Mcmcmc::computeBeta(double d, size_t idx)
{
    
    return 1.0 / (1.0+delta*idx);
}



void Mcmcmc::disableScreenMonitor( bool all, size_t rep )
{
    
    for (size_t i = 0; i < num_chains; ++i)
    {
        
        if ( chains[i] != NULL )
        {
            chains[i]->disableScreenMonitor(all, rep);
        }
        
    }
    
}


/**
 * Start the monitors at the beginning of a run which will simply delegate this call to each chain.
 */
void Mcmcmc::finishMonitors( size_t n_reps, MonteCarloAnalysisOptions::TraceCombinationTypes tc )
{
    
    // Monitor
    for (size_t i = 0; i < num_chains; ++i)
    {
        
        if ( chains[i] != NULL )
        {
            chains[i]->finishMonitors( n_reps, tc );
        }
    }
    
}



/**
 * Get the model instance.
 */
const Model& Mcmcmc::getModel( void ) const
{
    
    return base_chain->getModel();
}


double Mcmcmc::getModelLnProbability( bool likelihood_only )
{
    // we need to make sure that the vector chain_values is propagated properly
    // usualy we store the posteriors in there
    // now we might want to get the likelihoods only
    synchronizeValues(likelihood_only);
    
    // create the return value
    double rv = RbConstants::Double::neginf;
    
    for (size_t i=0; i<num_chains; ++i)
    {
        if ( isColdChain(i) )
        {
            rv = chain_values[i];
            break;
        }
    }
    
    // we need to make sure that the vector chain_values is propagated properly
    // usualy we store the posteriors in there
    // now we might want to get the likelihoods only
    synchronizeValues(false);
    
    return rv;
}


RbVector<Monitor>& Mcmcmc::getMonitors( void )
{
    RbVector<Monitor> *monitors = new RbVector<Monitor>();
    for (size_t i = 0; i < num_chains; ++i)
    {
        if ( chains[i] != NULL )
        {
            RbVector<Monitor>& m = chains[i]->getMonitors();
            for (size_t j = 0; j < m.size(); ++j)
            {
                monitors->push_back( m[j] );
            }
        }
    }
    return *monitors;
}


std::string Mcmcmc::getStrategyDescription( void ) const
{
    std::string description = "";
    std::stringstream stream;
    stream << "The MCMCMC simulator runs 1 cold chain and " << (num_chains-1) << " heated chains.\n";
    //    stream << chains[ chainsPerProcess[pid][0] ]->getStrategyDescription();
    size_t chain_index = 0;
    while ( chain_index < num_chains && chains[chain_index] == NULL ) ++chain_index;
    
    stream << chains[chain_index]->getStrategyDescription();
    description = stream.str();
    
    return description;
}

double Mcmcmc::heatForIndex(size_t i) const
{
    // Can we just maintain a sorted list of chain heats that is always valid?
    std::vector<double> tmp_chain_heats = chain_heats;
    sort(tmp_chain_heats.begin(), tmp_chain_heats.end(), std::greater<double>());
    return tmp_chain_heats[i];
}

double Mcmcmc::heatForChain(size_t i) const
{
    return chain_heats[i];
}

size_t Mcmcmc::heatIndexForChain(size_t j) const
{
    std::vector<double> tmp_chain_heats = chain_heats;
    std::sort(tmp_chain_heats.begin(), tmp_chain_heats.end(), std::greater<double>());

    return std::find(tmp_chain_heats.begin(), tmp_chain_heats.end(), chain_heats[j]) - tmp_chain_heats.begin();
}

bool Mcmcmc::isColdChain(size_t i) const
{
    return heatForChain(i) == 1.0;
}


void Mcmcmc::initializeChains(void)
{
    
    double processors_per_chain = double(num_processes) / double(num_chains);

    for (size_t i = 0; i < num_chains; ++i)
    {
        // all chains know heat-order and chain-processor schedules
        heat_ranks[i] = i;
        
        // get chain heat
        if (heat_temps[0] != 0.0)
        {
            chain_heats[i] = heat_temps[i];
        }
        else
        {
            double b = computeBeta(delta,i);
            chain_heats[i] = b;
        }
        
        chain_moves_tuningInfo[i]       = base_chain->getMovesTuningInfo();
        
        size_t active_pid_for_chain     = size_t( floor( i     * processors_per_chain ) + active_PID);
        size_t num_processer_for_chain  = size_t( floor( (i+1) * processors_per_chain ) + active_PID) - active_pid_for_chain;
        if ( num_processer_for_chain < 1 )
        {
            num_processer_for_chain = 1;
        }
        pid_per_chain[i] = active_pid_for_chain;
        

        // add chain to pid's chain vector (smaller memory footprint)
        if ( pid >= active_pid_for_chain && pid < (active_pid_for_chain + num_processer_for_chain) )
        {

            // create chains
            Mcmc* oneChain = new Mcmc( *base_chain );
            oneChain->setScheduleType( schedule_type );
            oneChain->setChainActive( i == 0 );
            oneChain->setChainPosteriorHeat( heatForChain(i) );
            oneChain->setChainIndex( i );
            oneChain->setActivePID( active_pid_for_chain, num_processer_for_chain );
            chains[i] = oneChain;
        }
        else
        {
            chains[i] = NULL;
        }
    }
    
}


void Mcmcmc::initializeSampler()
{
    
    // initialize each chain
    for (size_t i = 0; i < num_chains; ++i)
    {
        
        if ( chains[i] != NULL )
        {
            chains[i]->initializeSampler();
        }
    }
    
}



void Mcmcmc::initializeSamplerFromCheckpoint( void )
{
    
    for (size_t i = 0; i < num_chains; ++i)
    {
            
        if ( chains[i] != NULL )
        {
            chains[i]->initializeSamplerFromCheckpoint();
            setCurrentGeneration(chains[i]->getCurrentGeneration());
        }
        
    }
    
}


void Mcmcmc::monitor(unsigned long g)
{
    
    for (size_t i = 0; i < num_chains; ++i)
    {
        
        if ( chains[i] != NULL && chains[i]->isChainActive() )
        {
            chains[i]->monitor(g);
        }
        
    }
    
}

void Mcmcmc::nextCycle(bool advanceCycle)
{
    
    // run each chain for this process
    for (size_t i = 0; i < num_chains; ++i)
    {
        
        if ( chains[i] != NULL )
        {
            // @todo: #thread
            // This part should be done on several threads if possible
            // Sebastian: this call is very slow; a lot of work happens in nextCycle()

            // advance chain i by a single cycle
            chains[i]->nextCycle( advanceCycle );
        }
        
    } // loop over chains for this process
    
    if ( advanceCycle == true )
    {
        // advance gen counter
        ++current_generation;
    }
    else
    {
        ++burnin_generation;
    }
    
    if (useNeighborSwapping == true && useRandomSwapping == true)
    {
        if ((current_generation == 0 && burnin_generation % swap_interval == 0) || (current_generation > 0 && current_generation % swap_interval == 0))
        {
            
            // perform chain swap
            if (swap_mode == "single")
            {
                swapChains("neighbor");
            }
            else if (swap_mode == "multiple")
            {
                for (size_t i = 0; i < num_chains - 1; ++i)
                {
                    swapChains("neighbor");
                }
            }
            
        }
        if ((current_generation == 0 && burnin_generation % swap_interval2 == 0) || (current_generation > 0 && current_generation % swap_interval2 == 0))
        {
            
            // perform chain swap
            if (swap_mode == "single")
            {
                swapChains("random");
            }
            else if (swap_mode == "multiple")
            {
                for (size_t i = 0; i < num_chains - 1; ++i)
                {
                    for (size_t j = i + 1; j < num_chains; ++j)
                    {
                        swapChains("random");
                    }
                }
            }
            
        }
    }
    else if (useNeighborSwapping == true)
    {
        if ((current_generation == 0 && burnin_generation % swap_interval == 0) || (current_generation > 0 && current_generation % swap_interval == 0))
        {
            
            // perform chain swap
            if (swap_mode == "single")
            {
                swapChains("neighbor");
            }
            else if (swap_mode == "multiple")
            {
                for (size_t i = 0; i < num_chains - 1; ++i)
                {
                    swapChains("neighbor");
                }
            }
            
        }
    }
    else if (useRandomSwapping == true)
    {
        if ((current_generation == 0 && burnin_generation % swap_interval == 0) || (current_generation > 0 && current_generation % swap_interval == 0))
        {
            
            // perform chain swap
            if (swap_mode == "single")
            {
                swapChains("random");
            }
            else if (swap_mode == "multiple")
            {
                for (size_t i = 0; i < num_chains - 1; ++i)
                {
                    for (size_t j = i + 1; j < num_chains; ++j)
                    {
                        swapChains("random");
                    }
                }
            }
            
        }
    }
    
    for(size_t i=0;i<num_chains;i++)
    {
        size_t heat_rank = heatIndexForChain(i);
        if (chain_prev_boundary[i] == boundary::coldest)
            heat_visitors[heat_rank].first++;
        if (chain_prev_boundary[i] == boundary::hottest)
            heat_visitors[heat_rank].second++;
    }
    
    setCurrentGeneration(current_generation);
}


/**
 * Print the summary of the move.
 *
 * The summary just contains the current value of the tuning parameter.
 * It is printed to the stream that it passed in.
 *
 * \param[in]     o     The stream to which we print the summary.
 */
void Mcmcmc::printMoveSummary(std::ostream &o, size_t chainId, size_t moveId, Move &mv, bool current_period) const
{
    
    std::streamsize previousPrecision = o.precision();
    std::ios_base::fmtflags previousFlags = o.flags();
    
    o << std::fixed;
    o << std::setprecision(4);
    
    // print the name
    const std::string &mv_name = mv.getMoveName();
    size_t spaces = 40 - (mv_name.length() > 40 ? 40 : mv_name.length());
    o << mv_name;
    for (size_t i = 0; i < spaces; ++i)
    {
        o << " ";
    }
    o << " ";
    
    // print the DagNode name
    const std::string &dn_name = (*mv.getDagNodes().begin())->getName();
    spaces = 20 - (dn_name.length() > 20 ? 20 : dn_name.length());
    o << dn_name;
    for (size_t i = 0; i < spaces; ++i)
    {
        o << " ";
    }
    o << " ";
    
    // print the weight
    int w_length = 4;
    const double weight = mv.getUpdateWeight();
    if (weight > 0) w_length -= (int)log10(weight);
    for (int i = 0; i < w_length; ++i)
    {
        o << " ";
    }
    o << weight;
    o << " ";
    
    // print the number of tries
    int t_length = 9;
    size_t num_tried = chain_moves_tuningInfo[chainId][moveId].num_tried_total;
    if (current_period == true)
    {
        num_tried = chain_moves_tuningInfo[chainId][moveId].num_tried_current_period;
    }
    
    if (num_tried > 0) t_length -= (int)log10(num_tried);
    for (int i = 0; i < t_length; ++i)
    {
        o << " ";
    }
    o << num_tried;
    o << " ";
    
    // print the number of accepted
    int a_length = 9;
    size_t num_accepted = chain_moves_tuningInfo[chainId][moveId].num_accepted_total;
    if (current_period == true)
    {
        num_accepted = chain_moves_tuningInfo[chainId][moveId].num_accepted_current_period;
    }
    
    if (num_accepted > 0) a_length -= (int)log10(num_accepted);
    
    for (int i = 0; i < a_length; ++i)
    {
        o << " ";
    }
    o << num_accepted;
    o << " ";
    
    // print the acceptance ratio
    double ratio;
    if (num_tried == 0)
    {
        ratio = 0;
    }
    else
    {
        ratio = num_accepted / (double)num_tried;
    }
    int r_length = 5;
    
    for (int i = 0; i < r_length; ++i)
    {
        o << " ";
    }
    o << ratio;
    o << " ";
    
    
    if (RbMath::isNan(chain_moves_tuningInfo[chainId][moveId].tuning_parameter) == false)
    {
        MetropolisHastingsMove *move = dynamic_cast<MetropolisHastingsMove*>(&mv);
        
        if (move != NULL)
        {
            move->getProposal().printParameterSummary( o, true );
        }
        
        o << chain_moves_tuningInfo[chainId][moveId].tuning_parameter;
    }
    
    o << std::endl;
    o.flush();
    
    o.setf(previousFlags);
    o.precision(previousPrecision);
    
}


void Mcmcmc::printOperatorSummary(bool current_period)
{
    
    // send all chain heats to pid 0
    synchronizeHeats();
    
    // send all chain moves tuning information to pid 0
    synchronizeTuningInfo();
    
    if (process_active == true)
    {
        RbVector<Move> base_moves;
        size_t active_chainIdx = std::find(pid_per_chain.begin(), pid_per_chain.end(), active_PID) - pid_per_chain.begin();
        
        if (chains[active_chainIdx] != NULL)
        {
            base_moves = chains[active_chainIdx]->getMoves();
        }
        
        for (size_t i = 0; i < num_chains; ++i)
        {
            size_t chainIdx = chainForHeatIndex(i);
            
            std::cout << std::endl;
            std::cout << "chain_heats[" << i + 1 << "] = " << heatForIndex(i) << std::endl;
            std::cout.flush();
            
            // printing the moves summary
            std::cout << std::endl;
            std::cout << "                  Name                  | Param              |  Weight  |  Tried   | Accepted | Acc. Ratio| Parameters" << std::endl;
            std::cout << "===============================================================================================================================" << std::endl;
            
            for (size_t j = 0; j < chain_moves_tuningInfo[0].size(); ++j)
            {
                printMoveSummary(std::cout, chainIdx, j, base_moves[j], current_period);
            }
            
            std::cout << std::endl;
            std::cout.flush();
        }
        
        if (num_chains > 1)
        {
            printSwapSummary(std::cout);
            std::cout << std::endl;
            std::cout.flush();

            printTripSummary(std::cout);
            std::cout << std::endl;
            std::cout.flush();

            printHeatSummary(std::cout);
            std::cout << std::endl;
            std::cout.flush();
        }
    }

}


/**
 * Print the summary of the move.
 *
 * The summary just contains the current value of the tuning parameter.
 * It is printed to the stream that it passed in.
 *
 * \param[in]     o     The stream to which we print the summary.
 */
void Mcmcmc::printTripSummary(std::ostream &o) const
{
    std::streamsize previousPrecision = o.precision();
    std::ios_base::fmtflags previousFlags = o.flags();

    o << std::fixed;
    o << std::setprecision(4);

    o << std::endl;
    o << "MCMCMC chain round trips |  Number  |  Number/iteration   " << std::endl;
    o << "=======================================================" << std::endl;

    for (size_t i = 0; i < num_chains; ++i)
    {
        // i -> c -> h -> c  -> 2/2 = 1
        // i -> h -> c -> h  -> 1/2 = 0
        size_t round_trips = chain_half_trips[i];
        if (chain_prev_boundary[i] == boundary::hottest)
            round_trips--;
        round_trips /= 2;

        o<<std::setw(25)<<i<<std::setw(11)<<round_trips<<std::setw(19)<<double(round_trips)/(current_generation+burnin_generation)<< std::endl;
    }

    o.setf(previousFlags);
    o.precision(previousPrecision);
}


void Mcmcmc::printHeatSummary(std::ostream &o) const
{
    std::streamsize previousPrecision = o.precision();
    std::ios_base::fmtflags previousFlags = o.flags();

    o << std::fixed;
    o << std::setprecision(4);

    o << std::endl;
    o << "Heat rank |  Heat  | Temperature | Fraction" << std::endl;
    o << "===========================================" << std::endl;

    for (size_t i = 0; i < num_chains; ++i)
    {
        size_t c = heat_visitors[i].first;
        size_t h = heat_visitors[i].second;
        double B = heatForIndex(i);
        double T = 1.0/B;
        if (c + h > 0)
        {
            double fraction = double(c)/(c+h);
            o<<std::setw(10)<<i<<std::setw(9)<<B<<std::setw(14)<<T<<std::setw(10)<<fraction<<std::endl;
        }
        else
            o<<std::setw(10)<<i<<std::setw(9)<<B<<std::setw(14)<<T<<std::setw(10)<<"NA"<<std::endl;
    }

    o.setf(previousFlags);
    o.precision(previousPrecision);
}

/**
 * Print the summary of the move.
 *
 * The summary just contains the current value of the tuning parameter.
 * It is printed to the stream that it passed in.
 *
 * \param[in]     o     The stream to which we print the summary.
 */
void Mcmcmc::printSwapSummary(std::ostream &o) const
{
    
    std::streamsize previousPrecision = o.precision();
    std::ios_base::fmtflags previousFlags = o.flags();
    
    o << std::fixed;
    o << std::setprecision(4);
    
    std::cout << std::endl;
    
    if (useNeighborSwapping == true && useRandomSwapping == true)
    {
        std::cout << "MCMCMC chains swapping between|swapIntervalNeighbor|swapIntervalRandom| Tried | Accepted | Acc. Ratio |  HeatFrom  |  HeatTo   " << std::endl;
        std::cout << "===============================================================================================================================" << std::endl;
        
        for (size_t i = 0; i < num_chains - 1; ++i)
        {
            for (size_t j = i + 1; j < num_chains; ++j)
            {
                printSwapSummaryPair(o, i, j);
                printSwapSummaryPair(o, j, i);
            }
        }

    }
    else
    {
        std::cout << "MCMCMC chains swapping between|              swapInterval             | Tried | Accepted | Acc. Ratio |  HeatFrom  |  HeatTo   " << std::endl;
        std::cout << "===============================================================================================================================" << std::endl;
        
        if (useRandomSwapping == true)
        {
            for (size_t i = 0; i < num_chains - 1; ++i)
            {
                for (size_t j = i + 1; j < num_chains; ++j)
                {
                    printSwapSummaryPair(o, i, j);
                    printSwapSummaryPair(o, j, i);
                }
            }
        }
        else if (useNeighborSwapping == true)
        {
            for (size_t i = 0; i < num_chains - 1; ++i)
            {
                printSwapSummaryPair(o, i, i + 1);
            }
            
            for (size_t i = 1; i < num_chains; ++i)
            {
                printSwapSummaryPair(o, i, i - 1);
            }
        }

    }
    
    o.setf(previousFlags);
    o.precision(previousPrecision);

}


void Mcmcmc::printSwapSummaryPair(std::ostream &o, const size_t &row, const size_t &col) const
{
    // print the name
    o << row + 1;
    const std::string n1 = " to ";
    o << n1;
    o << col + 1;
    
    size_t n_length = n1.length() + (size_t)log10(row + 1) + (size_t)log10(col + 1);
    size_t spaces = 30 - (n_length > 30 ? 30 : n_length);
    for (size_t i = 0; i < spaces; ++i)
    {
        o << " ";
    }
    o << " ";
    
    // print the swap interval
    if (swap_interval2 > 0)
    {
        int si_length = 18;
        if (swap_interval > 0) si_length -= (int)log10(swap_interval);
        for (int i = 0; i < si_length; ++i)
        {
            o << " ";
        }
        o << swap_interval;
        o << " ";
        
        int si2_length = 16 - (int)log10(swap_interval2);
        for (int i = 0; i < si2_length; ++i)
        {
            o << " ";
        }
        o << swap_interval2;
        o << " ";
    }
    else if (swap_interval2 == 0)
    {
        int si_length = 36;
        if (swap_interval > 0) si_length -= (int)log10(swap_interval);
        for (int i = 0; i < si_length; ++i)
        {
            o << " ";
        }
        o << swap_interval;
        o << " ";
    }
    
    // print the number of tries
    int t_length = 6;
    if (num_attempted_swaps[row][col] > 0) t_length -= (int)log10(num_attempted_swaps[row][col]);
    for (int i = 0; i < t_length; ++i)
    {
        o << " ";
    }
    o << num_attempted_swaps[row][col];
    o << " ";
    
    // print the number of accepted
    int a_length = 9;
    if (num_accepted_swaps[row][col] > 0) a_length -= (int)log10(num_accepted_swaps[row][col]);
    for (int i = 0; i < a_length; ++i)
    {
        o << " ";
    }
    o << num_accepted_swaps[row][col];
    o << " ";
    
    // print the acceptance ratio
    double ratio = num_accepted_swaps[row][col] / (double)num_attempted_swaps[row][col];
    if (num_attempted_swaps[row][col] == 0) ratio = 0;
    
    int r_length = 6;
    for (int i = 0; i < r_length; ++i)
    {
        o << " ";
    }
    o << ratio;
    o << " ";
    
    // print the heat of the chain that swaps are proposed from
    int h_length = 6;
    for (int i = 0; i < h_length; ++i)
    {
        o << " ";
    }
    o << heatForIndex(row);
    o << " ";
    
    // print the heat of the chain that swaps are proposed to
    h_length = 5;
    for (int i = 0; i < h_length; ++i)
    {
        o << " ";
    }
    o << heatForIndex(col);
    o << " ";
    
    o << std::endl;
    
}


void Mcmcmc::redrawStartingValues( void )
{
    
    // initialize each chain
    for (size_t i = 0; i < num_chains; ++i)
    {
        RandomNumberGenerator *rng = GLOBAL_RNG;
        for (size_t j=0; j<10; ++j) rng->uniform01();
        
        if ( chains[i] != NULL )
        {
            chains[i]->redrawStartingValues();
        }
    }
    
}


void Mcmcmc::removeMonitors( void )
{
    
    for (size_t i = 0; i < num_chains; ++i)
    {
        
        if ( chains[i] != NULL )
        {
            chains[i]->removeMonitors();
        }
        
    }
    
}


void Mcmcmc::reset( void )
{
    
    // reset counters
    resetCounters();
    
    //    /* Reset the monitors */
    //    for (size_t i = 0; i < chainsPerProcess[pid].size(); ++i)
    //    {
    //        RbVector<Monitor>& monitors = chains[ chainsPerProcess[pid][i] ]->getMonitors();
    //        for (size_t i=0; i<monitors.size(); ++i)
    //        {
    //            monitors[i].reset();
    //        }
    //    }
    
    for (size_t i = 0; i < num_chains; ++i)
    {
        
        if ( chains[i] != NULL )
        {
            chains[i]->reset();
        }
        
    }
    
    
}


void Mcmcmc::resetCounters( void )
{
    for (size_t i = 0; i < num_chains; ++i)
    {
        for (size_t j = 0; j < num_chains; ++j)
        {
            num_attempted_swaps[i][j] = 0;
            num_accepted_swaps[i][j] = 0;
//            std::cout << "num_attempted_swaps[" << i << "][" << j << "]=" << num_attempted_swaps[i][j] << ", num_accepted_swaps[" << i << "][" << j << "]=" << num_accepted_swaps[i][j] << "; ";
        }
    }
//    std::cout << std::endl;
}


void Mcmcmc::setHeatsInitial(const std::vector<double>& ht)
{
    // 1. Check that there are num_chains heats.
    if (ht.size() != num_chains)
        throw RbException()<<"Heat vector contains "<<ht.size()<<" heats, but there are "<<num_chains<<" chains.";

    // 2. Check that the heats are non-zero
    for(size_t i=0; i<ht.size(); i++)
    {
        if (not (ht[i] > 0.0))
            throw RbException()<<"Heat "<<i+1<<" is "<<ht[i]<<", but should be > 0";
    }

    // 3. Check that the heats are sorted with the largest heat first (== smallest temperature first).
    for(size_t i=0; i+1 < ht.size(); i++)
        if (ht[i] < ht[i+1])
            throw RbException()<<"Heat "<<i+1<<" ("<<ht[i]<<") should be greater than heat "<<i+2<<" ("<<ht[i+1]<<")";

    // 4. Check that the first heat is 1
    if (ht[0] != 1.0)
        throw RbException()<<"Heat 1 should be 1.0 (the cold chain), but is actually "<<ht[0];

    heat_temps = ht;
}


void Mcmcmc::setSwapInterval2(const size_t &si2)
{
    swap_interval2 = si2;
}


/**
 * Set the heat of the likelihood of the current chain.
 * This heat is used in posterior posterior MCMC algorithms to
 * heat the likelihood
 * The heat is passed to the moves for the accept-reject mechanism.
 */
void Mcmcmc::setLikelihoodHeat(double h)
{
    
    for (size_t i = 0; i < num_chains; ++i)
    {
        if (chains[i] != NULL)
        {
            chains[i]->setLikelihoodHeat( h );
        }
        
    }
    
}


/**
  * Set the model by delegating the model to the chains.
  */
void Mcmcmc::setModel( Model *m, bool redraw )
{
    
    // set the models of the chains
    for (size_t i = 0; i < num_chains; ++i)
    {
        if ( chains[i] != NULL )
        {
            Model *m_clone = m->clone();
            chains[i]->setModel( m_clone, redraw );
        }
        
    }
    
    if ( base_chain != NULL )
    {
        Model *m_clone = m->clone();
        base_chain->setModel( m_clone, redraw );
    }
    
    delete m;
    
}


void Mcmcmc::setActivePIDSpecialized(size_t i, size_t n)
{
    
    // initialize container sizes
    for (size_t i = 0; i < chains.size(); ++i)
    {
        
        if (chains[i] != NULL)
        {
            delete chains[i];
        }
        
    }
    
    chains.clear();
    chain_values.clear();
    chain_heats.clear();
    chain_prev_boundary.clear();
    chain_half_trips.clear();
    heat_visitors.clear();
    heat_ranks.clear();
    chain_moves_tuningInfo.clear();
    
    chains.resize(num_chains);
    chain_values.resize(num_chains, 0.0);
    chain_heats.resize(num_chains, 0.0);
    chain_prev_boundary.resize(num_chains, boundary::intermediate);
    chain_half_trips.resize(num_chains, 0);
    heat_visitors.resize(num_chains, {0,0});
    pid_per_chain.resize(num_chains, 0);
    heat_ranks.resize(num_chains, 0);
    
    chain_moves_tuningInfo = std::vector< std::vector<Mcmc::tuningInfo> > (num_chains);
    
    initializeChains();
}


/**
 * Start the monitors at the beginning of a run which will simply delegate this call to each chain.
 */
void Mcmcmc::startMonitors(size_t num_cycles, bool reopen)
{
    
    // Monitor
    for (size_t i = 0; i < num_chains; ++i)
    {
        
        if ( chains[i] != NULL )
        {
            chains[i]->startMonitors( num_cycles, reopen );
        }
        
    }
    
}


void Mcmcmc::setCheckpointFile(const path &f)
{
    
    for (size_t j = 0; j < num_chains; ++j)
    {

        if ( chains[j] != NULL )
        {
            chains[j]->setCheckpointFile(f);
        }
        
    }
    
}


void Mcmcmc::synchronizeValues( bool likelihood_only )
{
    
    // synchronize chain values
    double results[num_chains];
    for (size_t j = 0; j < num_chains; ++j)
    {
        results[j] = 0.0;
    }
    for (size_t j = 0; j < num_chains; ++j)
    {
        
        if ( chains[j] != NULL )
        {
            results[j] = chains[j]->getModelLnProbability(likelihood_only);
//            std::cout << "results[" << j << "]=" << results[j] << ", ";
        }
        
    }
//    std::cout << std::endl;
    
#ifdef RB_MPI
    if ( active_PID != pid )
    {
        for (size_t i=0; i<num_chains; ++i)
        {
            if ( pid == pid_per_chain[i] )
            {
                MPI_Send(&results[i], 1, MPI_DOUBLE, active_PID, 0, MPI_COMM_WORLD);
            }
            
        }
        
    }
#endif
    
    if ( active_PID == pid )
    {
#ifdef RB_MPI
        
        for (size_t j = 0; j < num_chains; ++j)
        {
            
            // ignore self
            if (pid != pid_per_chain[j])
            {
                MPI_Status status;
                MPI_Recv(&results[j], 1, MPI_DOUBLE, pid_per_chain[j], 0, MPI_COMM_WORLD, &status);
            }
            
        }
#endif
        for (size_t i = 0; i < num_chains; ++i)
        {
            chain_values[i] = results[i];
        }
    }
    
#ifdef RB_MPI
    if ( active_PID == pid )
    {
        for (size_t i=1; i<num_processes; ++i)
        {
            for (size_t j=0; j<num_chains; ++j)
            {
                MPI_Send(&chain_values[j], 1, MPI_DOUBLE, active_PID+i, 0, MPI_COMM_WORLD);
            }
        }
    }
    else
    {
        for (size_t i=0; i<num_chains; ++i)
        {
            MPI_Status status;
            MPI_Recv(&chain_values[i], 1, MPI_DOUBLE, active_PID, 0, MPI_COMM_WORLD, &status);
        }
        
    }
#endif
    
}

void Mcmcmc::synchronizeHeats(void)
{
    
    // synchronize heat values
    double heats[num_chains];
    for (size_t j = 0; j < num_chains; ++j)
    {
        heats[j] = 0.0;
    }
    for (size_t j = 0; j < num_chains; ++j)
    {
        if (chains[j] != NULL)
        {
            heats[j] = chains[j]->getChainPosteriorHeat();
        }
    }
    
#ifdef RB_MPI
    // share the heats accross processes
    if ( active_PID != pid )
    {
        for (size_t i=0; i<num_chains; ++i)
        {
            if ( pid == pid_per_chain[i] )
            {
                MPI_Send(&heats[i], 1, MPI_DOUBLE, active_PID, 0, MPI_COMM_WORLD);
            }
            
        }
        
    }
#endif
    
    if ( active_PID == pid )
    {
#ifdef RB_MPI
        for (size_t j = 0; j < num_chains; ++j)
        {
            
            // ignore self
            if (pid != pid_per_chain[j])
            {
                MPI_Status status;
                MPI_Recv(&heats[j], 1, MPI_DOUBLE, pid_per_chain[j], 0, MPI_COMM_WORLD, &status);
            }
            
        }
        
#endif
        for (size_t i = 0; i < num_chains; ++i)
        {
            chain_heats[i] = heats[i];
        }
    }
#ifdef RB_MPI
    if ( active_PID == pid )
    {
        for (size_t i=1; i<num_processes; ++i)
        {
            
            for (size_t j = 0; j < num_chains; ++j)
            {
                MPI_Send(&chain_heats[j], 1, MPI_DOUBLE, active_PID+i, 0, MPI_COMM_WORLD);
            }
        }
    }
    else
    {
        for (size_t i=0; i<num_chains; ++i)
        {
            MPI_Status status;
            MPI_Recv(&chain_heats[i], 1, MPI_DOUBLE, active_PID, 0, MPI_COMM_WORLD, &status);
        }
        
    }
#endif
    
}


void Mcmcmc::synchronizeTuningInfo(void)
{
    
    // synchronize tuning information
    std::vector< std::vector<Mcmc::tuningInfo> > chain_mvs_ti = std::vector< std::vector<Mcmc::tuningInfo> > (num_chains);
    for (size_t j = 0; j < num_chains; ++j)
    {
        if (chains[j] != NULL)
        {
            chain_mvs_ti[j] = chains[j]->getMovesTuningInfo();
        }
        else
        {
            chain_mvs_ti[j] = base_chain->getMovesTuningInfo();
        }
    }
    
#ifdef RB_MPI
    // create MPI type
    const int num_items = 5;
    int block_lengths[num_items] = {1, 1, 1, 1, 1};
    MPI_Datatype types[num_items] = {MPI_UNSIGNED_LONG, MPI_UNSIGNED_LONG, MPI_UNSIGNED_LONG, MPI_UNSIGNED_LONG, MPI_DOUBLE};
    MPI_Datatype tmp_mvs_ti_types, mvs_ti_types;
    MPI_Aint offsets[num_items];
    
    offsets[0] = offsetof(Mcmc::tuningInfo, num_tried_current_period);
    offsets[1] = offsetof(Mcmc::tuningInfo, num_tried_total);
    offsets[2] = offsetof(Mcmc::tuningInfo, num_accepted_current_period);
    offsets[3] = offsetof(Mcmc::tuningInfo, num_accepted_total);
    offsets[4] = offsetof(Mcmc::tuningInfo, tuning_parameter);
    
    MPI_Type_create_struct(num_items, block_lengths, offsets, types, &tmp_mvs_ti_types);
    
    MPI_Aint lb, extent;
    MPI_Type_get_extent(tmp_mvs_ti_types, &lb, &extent);
    
    MPI_Type_create_resized(tmp_mvs_ti_types, lb, extent, &mvs_ti_types);
    MPI_Type_commit(&mvs_ti_types);
    
    // share the heats accross processes
    if ( active_PID != pid )
    {
        for (size_t i=0; i<num_chains; ++i)
        {
            if ( pid == pid_per_chain[i] )
            {
                MPI_Send(&chain_mvs_ti[i].front(), chain_mvs_ti[i].size(), mvs_ti_types, active_PID, 0, MPI_COMM_WORLD);
            }
            
        }
        
    }
#endif
    
    if ( active_PID == pid )
    {
#ifdef RB_MPI
        for (size_t j = 0; j < num_chains; ++j)
        {
            
            // ignore self
            if (pid != pid_per_chain[j])
            {
                MPI_Status status;
                MPI_Recv(&chain_mvs_ti[j].front(), chain_mvs_ti[j].size(), mvs_ti_types, pid_per_chain[j], 0, MPI_COMM_WORLD, &status);
            }
            
        }
        
#endif
        for (size_t i = 0; i < num_chains; ++i)
        {
            chain_moves_tuningInfo[i] = chain_mvs_ti[i];
        }
    }
    
#ifdef RB_MPI
    if ( active_PID == pid )
    {
        for (size_t i=1; i<num_processes; ++i)
        {
            
            for (size_t j = 0; j < num_chains; ++j)
            {
                MPI_Send(&chain_moves_tuningInfo[j].front(), chain_moves_tuningInfo[j].size(), mvs_ti_types, active_PID+i, 0, MPI_COMM_WORLD);
            }
        }
    }
    else
    {
        for (size_t i=0; i<num_chains; ++i)
        {
            MPI_Status status;
            MPI_Recv(&chain_moves_tuningInfo[i].front(), chain_moves_tuningInfo[i].size(), mvs_ti_types, active_PID, 0, MPI_COMM_WORLD, &status);
        }
        
    }
    
    // free the datatype
    MPI_Type_free(&mvs_ti_types);
    
#endif
    
}


// MJL: allow swapChains to take a swap function -- e.g. pairwise swap for 1..n-1
void Mcmcmc::swapChains(const std::string swap_method)
{
    
    // exit if there is only one chain
    if (num_chains < 2)
    {
        return;
    }
    
    // send all chain values to pid 0
    synchronizeValues( false );
    
    // send all chain heats to pid 0
    synchronizeHeats();
    
    // send all chain moves tuning information to pid 0
    synchronizeTuningInfo();
    
    // swap chains
    if (swap_method == "neighbor")
    {
        swapNeighborChains();
    }
    else if (swap_method == "random")
    {
        swapRandomChains();
    }

}


void Mcmcmc::swapMovesTuningInfo(RbVector<Move> &mvsj, RbVector<Move> &mvsk)
{
    
    RbIterator<Move> itk = mvsk.begin();
    for (RbIterator<Move> itj = mvsj.begin(); itj != mvsj.end(); ++itj, ++itk)
    {
        if (itj->getMoveName() != itk->getMoveName())
        {
            throw RbException( "The two moves whose tuning information is attempted to be swapped are not the same move as their names do not match." );
        }
        
        size_t tmp_numTriedCurrentPeriodj = itj->getNumberTriedCurrentPeriod();
        size_t tmp_numTriedCurrentPeriodk = itk->getNumberTriedCurrentPeriod();
        itj->setNumberTriedCurrentPeriod(tmp_numTriedCurrentPeriodk);
        itk->setNumberTriedCurrentPeriod(tmp_numTriedCurrentPeriodj);
        
        size_t tmp_numTriedTotalj = itj->getNumberTriedTotal();
        size_t tmp_numTriedTotalk = itk->getNumberTriedTotal();
        itj->setNumberTriedTotal(tmp_numTriedTotalk);
        itk->setNumberTriedTotal(tmp_numTriedTotalj);
        
        size_t tmp_numAcceptedCurrentPeriodj = itj->getNumberAcceptedCurrentPeriod();
        size_t tmp_numAcceptedCurrentPeriodk = itk->getNumberAcceptedCurrentPeriod();
        itj->setNumberAcceptedCurrentPeriod(tmp_numAcceptedCurrentPeriodk);
        itk->setNumberAcceptedCurrentPeriod(tmp_numAcceptedCurrentPeriodj);
        
        size_t tmp_numAcceptedTotalj = itj->getNumberAcceptedTotal();
        size_t tmp_numAcceptedTotalk = itk->getNumberAcceptedTotal();
        itj->setNumberAcceptedTotal(tmp_numAcceptedTotalk);
        itk->setNumberAcceptedTotal(tmp_numAcceptedTotalj);
        
        double tmp_tuningParameterj = itj->getMoveTuningParameter();
        double tmp_tuningParameterk = itk->getMoveTuningParameter();
        
        if ((std::isnan(tmp_tuningParameterj) == true && std::isnan(tmp_tuningParameterk) == false) || (std::isnan(tmp_tuningParameterj) == false && std::isnan(tmp_tuningParameterk) == true))
        {
            throw RbException( "The two moves whose tuning information is attempted to be swapped are not the same move as only one of them has tuning parameter." );
        }
        else if (std::isnan(tmp_tuningParameterj) == false)
        {
            itj->setMoveTuningParameter(tmp_tuningParameterk);
            itk->setMoveTuningParameter(tmp_tuningParameterj);
        }
    }
    
    if (itk != mvsk.end())
    {
        throw RbException( "The two moves objects whose tuning information is attempted to be swapped have different number of moves." );
    }
    
}


void Mcmcmc::swapNeighborChains(void)
{
    // randomly pick the indices of two chains
    
    size_t heat_rankj = size_t(GLOBAL_RNG->uniform01() * (num_chains - 1));
    size_t heat_rankk = heat_rankj;
    if (num_chains > 1)
    {
        heat_rankk = heat_rankj + 1;
    }
    
    if ( GLOBAL_RNG->uniform01() >= 0.5)
        std::swap( heat_rankj, heat_rankk );
    
    size_t j = chainForHeatIndex( heat_rankj );
    size_t k = chainForHeatIndex( heat_rankk );

    swapGivenChains(j, k);
}



void Mcmcmc::swapRandomChains(void)
{
    // randomly pick the indices of two chains
    size_t j = 0;
    size_t k = 0;
    
    // swap?
    bool accept = false;
    
    j = size_t(GLOBAL_RNG->uniform01() * num_chains);
    if (num_chains > 1)
    {
        do {
            k = size_t(GLOBAL_RNG->uniform01() * num_chains);
        }
        while(j == k);
    }
    
    swapGivenChains(j, k);
}

void Mcmcmc::updateTrips(size_t j)
{
    // I think this code assumes that there is only one coldest chain and one hottest chain.
    // For example, there are not two chains with heat == 1.0.
    // This is true, and would probably stay true.

    // The assumption is that we run this right after we have swapped j with another chain.
    // Therefore heatIndexForChain(j) is the NEW heat rank for chain j.

    if ( heatIndexForChain(j) == 0)
    {
        if (chain_prev_boundary[j] == boundary::hottest)
            chain_half_trips[j]++;

        chain_prev_boundary[j] = boundary::coldest;
    }
    else if (heatIndexForChain(j) == num_chains - 1)
    {
        if (chain_prev_boundary[j] == boundary::coldest)
            chain_half_trips[j]++;

        chain_prev_boundary[j] = boundary::hottest;
    }
}

void Mcmcmc::swapGivenChains(size_t j, size_t k, double lnProposalRatio)
{
    size_t heat_rankj = heatIndexForChain(j);
    size_t heat_rankk = heatIndexForChain(k);

    ++num_attempted_swaps[heat_rankj][heat_rankk];

    // compute exchange ratio
    double bj = chain_heats[j];
    double bk = chain_heats[k];
    double lnPj = chain_values[j];
    double lnPk = chain_values[k];
    double lnR = bj * (lnPk - lnPj) + bk * (lnPj - lnPk) + lnProposalRatio;
    //        std::cout << "bj=" << bj << ", bk=" << bk << ", lnPj=" << lnPj << ", lnPk=" << lnPk << ", lnR=" << lnR << std::endl;

    // determine whether we accept or reject the chain swap
    double u = GLOBAL_RNG->uniform01();
    bool accept = false;
    if (lnR >= 0)
    {
        accept = true;
    }
    else if (lnR < -100)
    {
        accept = false;
    }
    else if (u < exp(lnR))
    {
        accept = true;
    }
    else
    {
        accept = false;
    }


#ifdef RB_MPI
    if ( active_PID == pid )
    {
        for (size_t i = 1; i < num_processes; ++i)
        {
            MPI_Send(&j, 1, MPI_INT, active_PID+i, 0, MPI_COMM_WORLD);
            MPI_Send(&k, 1, MPI_INT, active_PID+i, 0, MPI_COMM_WORLD);
            MPI_Send(&accept, 1, MPI_C_BOOL, active_PID+i, 0, MPI_COMM_WORLD);
        }
    }
    else
    {
        MPI_Status status;
        MPI_Recv(&j, 1, MPI_INT, active_PID, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&k, 1, MPI_INT, active_PID, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&accept, 1, MPI_C_BOOL, active_PID, 0, MPI_COMM_WORLD, &status);
    }
#endif


    // on accept, swap beta values and active chains
    if (accept == true )
    {
        ++num_accepted_swaps[heat_rankj][heat_rankk];

        // swap active chain
        if (active_chain_index == j)
        {
            active_chain_index = k;
        }
        else if (active_chain_index == k)
        {
            active_chain_index = j;
        }

        std::swap( chain_heats[j], chain_heats[k] );
        std::swap( heat_ranks[j], heat_ranks[k] );
        std::swap( chain_moves_tuningInfo[j], chain_moves_tuningInfo[k] );

        for (size_t i=0; i<num_chains; ++i)
        {
            if ( chains[i] != NULL )
            {
                chains[i]->setChainPosteriorHeat( heatForChain(i) );
                chains[i]->setMovesTuningInfo(chain_moves_tuningInfo[i]);
                chains[i]->setChainActive( isColdChain(i) );
            }
        }

        // Update statistics on round trips (well, half trips) for the swapped chains.
        updateTrips(j);
        updateTrips(k);
    }

    // update the chains accross processes
    // this is necessary because only process 0 does the swap
    // all the other processes need to be told that there was a swap
    //    updateChainState(j);
    //    updateChainState(k);
}


void Mcmcmc::tune( void )
{
    // send all chain heats to pid 0
    synchronizeHeats();
    
    // send all chain moves tuning information to pid 0
    synchronizeTuningInfo();
    
    if (tune_heat == true && num_chains > 1)
    {
        std::vector<double> heats_diff(num_chains - 1, 0.0);
        
        for (size_t i = 1; i < num_chains; ++i)
        {
            heats_diff[i - 1] = heatForIndex(i-1) - heatForIndex(i);
        }
        
        for (size_t i = 1; i < num_chains; ++i)
        {
//            std::cout << ", num_accepted_swaps[" << i - 1 << "][" << i << "]=" << num_accepted_swaps[i - 1][i] << ", num_accepted_swaps[" << i << "][" << i - 1 << "]=" << num_accepted_swaps[i][i - 1] << std::endl;
//            std::cout << ", num_attempted_swaps[" << i - 1 << "][" << i << "]=" << num_attempted_swaps[i - 1][i] << ", num_attempted_swaps[" << i << "][" << i - 1 << "]=" << num_attempted_swaps[i][i - 1] << std::endl;
            
            size_t num_attempted_swaps_neighbor = num_attempted_swaps[i - 1][i] + num_attempted_swaps[i][i - 1];
            if ( num_attempted_swaps_neighbor > 2 ) {
                
                size_t num_accepted_swaps_neighbor = num_accepted_swaps[i - 1][i] + num_accepted_swaps[i][i - 1];
                double rate = num_accepted_swaps_neighbor / (double)num_attempted_swaps_neighbor;
                
                if ( rate > tune_heat_target )
                {
                    heats_diff[i - 1] *= (1.0 + (rate - tune_heat_target) / (1.0 - tune_heat_target));
                }
                else
                {
                    heats_diff[i - 1] /= (2.0 - rate / tune_heat_target);
                }
                
//                std::cout << "rate=" << rate << ", heats_diff[" << i - 1 << "]=" << heats_diff[i - 1] << std::endl;
            }
        }
        
        double heatMinBound = 0.01;
        size_t j = 1;
        size_t colderChainIdx = chainForHeatIndex(j-1);
        size_t hotterChainIdx = colderChainIdx;
        
        for (; j < num_chains; ++j)
        {
            hotterChainIdx = chainForHeatIndex(j);
            chain_heats[hotterChainIdx] = chain_heats[colderChainIdx] - heats_diff[j - 1];
            
            if (chain_heats[hotterChainIdx] < heatMinBound)
            {
                break;
            }
            else
            {
                colderChainIdx = hotterChainIdx;
            }
//            std::cout << "chain_heats[" << hotterChainIdx << "]=" << chain_heats[hotterChainIdx] << std::endl;
        }
        
        /* If the heat of a given hot chain is smaller than the minimum bound, we want to linearly interpolate
         * this heat and the heats of all the hotter chains to fall between the lowest heat that is greater
         * than the minimum bound and the minimum bound.
         * The formula for linearly interpolating between (x1,y1) and (x2,y2) is y(x) = y2 - (y2-y1)*[ (x-x1)/(x2-x1) ].
         *
         * Let "heat(i)" be a notational shortcut for chain_heats[ chainForHeatIndex(i) ], and let k be the heat
         * index of the coldest of such "bad" chains. Then, its heat can be calculated as
         *
         * log(heat(k)) = log(heat(k-1)) - [ log(heat(k-1)) - log(heatMinBound) ]*[ (k - (k-1)) / (num_chains - (k-1)) ]
         *              = log(heat(k-1)) - [ log(heat(k-1)) - log(heatMinBound) ]*[ 1 / (num_chains - (k-1)) ]
         *
         * Let us further denote heat(k-1) as "lowestGreaterThanMin", and k-1 as "m". Then,
         *
         * log(heat(k)) = log(lowestGreaterThanMin) - [ log(lowestGreaterThanMin) - log(heatMinBound) ]*[ 1 / (num_chains - m) ]
         *
         * Finally, let us denote the subtrahend on the right-hand side of the equation above by log(rho). Then,
         *
         *       rho = (lowestGreaterThanMin / heatMinBound)^[ 1 / (num_chains - m) ],
         *   heat(k) = lowestGreaterThanMin / rho
         * heat(k+1) = lowestGreaterThanMin / rho^2, etc.
         */
        
        std::vector<double> greaterThanMin;
        std::vector<size_t> badChainIdx;
        
        // std::cout << "The following should be sorted from largest to smallest:" << std::endl;
        for (size_t l = 0; l < num_chains; ++l)
        {
            // std::cout << chain_heats[ chainForHeatIndex(l) ] << std::endl;
            if ( chain_heats[ chainForHeatIndex(l) ] > heatMinBound )
            {
                greaterThanMin.push_back( chain_heats[ chainForHeatIndex(l) ] );
            }
            else
            {
                badChainIdx.push_back(l);
            }
        }
        
        // std::cout << "We have " << badChainIdx.size() << " bad chain(s)" << std::endl;
        if ( badChainIdx.size() != 0 )
        {
            double lowestGreaterThanMin = *std::min_element( greaterThanMin.begin(), greaterThanMin.end() );
            size_t k = *std::min_element( badChainIdx.begin(), badChainIdx.end() );
            size_t m = k - 1;
            
            double rho = pow(lowestGreaterThanMin / heatMinBound, 1.0 / (num_chains - m));
            // std::cout << "Current rho is: " << rho << std::endl;
            
            for (; k < num_chains; ++k)
            {
                // std::cout << "Attempting to set the heat of chain " << k << " to " << (lowestGreaterThanMin / pow(rho, k - m)) << std::endl;
                chain_heats[ chainForHeatIndex(k) ] = lowestGreaterThanMin / pow(rho, k - m);
            }
        }
        
        resetCounters();
    }

    
    for (size_t i=0; i<num_chains; ++i)
    {
        if ( chains[i] != NULL )
        {
            chains[i]->setChainPosteriorHeat( heatForChain(i) );
            chains[i]->setChainActive( isColdChain( i ) );
            
            chains[i]->tune();
        }
    }
    
}


void Mcmcmc::updateChainState(size_t j)
{
    
#ifdef RB_MPI
    // update heat
    if ( active_PID == pid )
    {
        for (size_t i = 1; i < num_processes; ++i)
        {
            MPI_Send(&chain_heats[j], 1, MPI_DOUBLE, active_PID+i, 0, MPI_COMM_WORLD);
        }
    }
    else
    {
        MPI_Status status;
        MPI_Recv(&chain_heats[j], 1, MPI_DOUBLE, active_PID, 0, MPI_COMM_WORLD, &status);
    }
#endif
    
    if ( chains[j] != NULL )
    {
        chains[j]->setChainPosteriorHeat( chain_heats[j] );
    }
    
    for (size_t i=0; i<num_chains; ++i)
    {
        if ( chains[i] != NULL )
        {
            chains[i]->setChainActive( isColdChain(i) );
        }
    }
    
}


/**
 * Start the monitors at the beginning of a run which will simply delegate this call to each chain.
 */
void Mcmcmc::writeMonitorHeaders( bool screen_monitor_only )
{
    
    // Monitor
    for (size_t i = 0; i < num_chains; ++i)
    {
        
        if ( chains[i] != NULL )
        {
            chains[i]->writeMonitorHeaders( screen_monitor_only );
        }
        
    }
    
}

