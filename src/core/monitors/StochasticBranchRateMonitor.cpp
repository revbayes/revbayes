#include "StochasticBranchRateMonitor.h"

#include <cstddef>
#include <ostream>
#include <vector>

#include "FastBirthDeathShiftProcess.h"
#include "GeneralizedLineageHeterogeneousBirthDeathSamplingProcess.h"
#include "StochasticNode.h"
#include "StateDependentSpeciationExtinctionProcess.h"
#include "Cloneable.h"
#include "Tree.h"
#include "TypedDistribution.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;


/* Constructor for state dependent birth death process */
StochasticBranchRateMonitor::StochasticBranchRateMonitor(StochasticNode<Tree>* ch, std::uint64_t g, const std::string &fname, const std::string &del) : VariableMonitor(ch, g, fname, del, false, false, false),
    cdbdp( ch )
{
    // the cdbdp is both the tree and character evolution model
    addVariable( cdbdp );
}



/**
 * Copy constructor.
 */
StochasticBranchRateMonitor::StochasticBranchRateMonitor( const StochasticBranchRateMonitor &m) : VariableMonitor( m ),
    cdbdp( m.cdbdp )
{
    
}


/**
 * Destructor.
 */
StochasticBranchRateMonitor::~StochasticBranchRateMonitor()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
StochasticBranchRateMonitor* StochasticBranchRateMonitor::clone(void) const
{
    
    return new StochasticBranchRateMonitor(*this);
}


/**
 * Monitor value at given generation.
 *
 * \param[in]   gen    The current generation.
 */
void StochasticBranchRateMonitor::monitorVariables(std::uint64_t gen)
{
    auto& separator = to<SeparatorFormat>(format)->separator;
    
    std::vector<double> speciation;
    std::vector<double> extinction;
    std::vector<double> sampling;
    std::vector<double> destructive_sampling;
    std::vector<std::int64_t>   n_shifts;

    // other variables for the birth-death-shift monitor 
    std::vector<double>        delta_speciation;
    std::vector<double>        delta_extinction;
    std::vector<std::int64_t>  n_speciation_shifts;
    std::vector<std::int64_t>  n_extinction_shifts;


    size_t num_nodes = cdbdp->getValue().getNumberOfNodes();
    std::vector<std::string> character_histories( num_nodes );

    StateDependentSpeciationExtinctionProcess *sse = dynamic_cast<StateDependentSpeciationExtinctionProcess*>( &cdbdp->getDistribution() );
    GeneralizedLineageHeterogeneousBirthDeathSamplingProcess *glhbdsp = dynamic_cast<GeneralizedLineageHeterogeneousBirthDeathSamplingProcess*>( &cdbdp->getDistribution() );
    FastBirthDeathShiftProcess *bds = dynamic_cast<FastBirthDeathShiftProcess*>( &cdbdp->getDistribution() );

    bool is_sse = sse != NULL;
    bool is_glhbdsp = glhbdsp != NULL;
    bool is_bds = bds != NULL;


    if (is_sse)
    {
        // draw stochastic character map
        sse->drawStochasticCharacterMap( character_histories );
        speciation = sse->getAverageSpeciationRatePerBranch();
        extinction = sse->getAverageExtinctionRatePerBranch();
        n_shifts   = sse->getNumberOfShiftEventsPerBranch();
    }

    if (is_glhbdsp)
    {
        glhbdsp->drawStochasticCharacterMap(character_histories, speciation, extinction, sampling, destructive_sampling, n_shifts);
    }

    if (is_bds)
    {
        bds->drawStochasticCharacterMap( character_histories );

        // get the variables we are interested in
        speciation = bds->getAverageSpeciationRatePerBranch();
        extinction = bds->getAverageExtinctionRatePerBranch();
        delta_speciation = bds->getDeltaSpeciationPerBranch();
        delta_extinction = bds->getDeltaExtinctionPerBranch();
        n_speciation_shifts = bds->getNumberOfSpeciationShiftEventsPerBranch();
        n_extinction_shifts = bds->getNumberOfExtinctionShiftEventsPerBranch();
    }

    bool found_distribution = is_sse | is_glhbdsp | is_bds;

    if (!found_distribution)
    {
        throw RbException("Could not cast the correct birth-death distribution (in StochasticBranchRateMonitor)");
    }

    
    // print to monitor file
    for (int i = 0; i < speciation.size(); i++)
    {
        // add a separator before every new element
        out_stream << separator;
        out_stream << speciation[i];
    }
    
    for (int i = 0; i < extinction.size(); i++)
    {
        out_stream << separator;
        out_stream << extinction[i];
    }

    if (is_sse | is_glhbdsp) {
     for (int i = 0; i < n_shifts.size(); i++)
        {
            out_stream << separator;
            out_stream << n_shifts[i];
        }   
    }else {
        for (int i = 0; i < delta_speciation.size(); i++)
        {
            out_stream << separator;
            out_stream << delta_speciation[i];
        }

        for (int i = 0; i < delta_extinction.size(); i++)
        {
            out_stream << separator;
            out_stream << delta_extinction[i];
        }

        
        for (int i = 0; i < n_speciation_shifts.size(); i++)
        {
            out_stream << separator;
            out_stream << n_speciation_shifts[i];
        }

        for (int i = 0; i < n_extinction_shifts.size(); i++)
        {
            out_stream << separator;
            out_stream << n_extinction_shifts[i];
        }
    }
}


/**
 * Print header for monitored values
 */
void StochasticBranchRateMonitor::printFileHeader()
{
    auto& separator = to<SeparatorFormat>(format)->separator;
    
    std::vector<double> speciation;
    std::vector<double> extinction;
    std::vector<double> sampling;
    std::vector<double> destructive_sampling;
    std::vector<std::int64_t>   n_shifts;

    // other variables for the birth-death-shift monitor 
    std::vector<double>        delta_speciation;
    std::vector<double>        delta_extinction;
    std::vector<std::int64_t>  n_speciation_shifts;
    std::vector<std::int64_t>  n_extinction_shifts;

    size_t num_nodes = cdbdp->getValue().getNumberOfNodes();
    std::vector<std::string> character_histories( num_nodes );

    StateDependentSpeciationExtinctionProcess *sse = dynamic_cast<StateDependentSpeciationExtinctionProcess*>( &cdbdp->getDistribution() );
    GeneralizedLineageHeterogeneousBirthDeathSamplingProcess *glhbdsp = dynamic_cast<GeneralizedLineageHeterogeneousBirthDeathSamplingProcess*>( &cdbdp->getDistribution() );
    FastBirthDeathShiftProcess *bds = dynamic_cast<FastBirthDeathShiftProcess*>( &cdbdp->getDistribution() );

    bool is_sse = sse != NULL;
    bool is_glhbdsp = glhbdsp != NULL;
    bool is_bds = bds != NULL;


    if (is_sse)
    {
        // draw stochastic character map
        sse->drawStochasticCharacterMap( character_histories );
        speciation = sse->getAverageSpeciationRatePerBranch();
        extinction = sse->getAverageExtinctionRatePerBranch();
        n_shifts   = sse->getNumberOfShiftEventsPerBranch();
    }

    if (is_glhbdsp)
    {
        glhbdsp->drawStochasticCharacterMap(character_histories, speciation, extinction, sampling, destructive_sampling, n_shifts);
    }

    if (is_bds)
    {
        bds->drawStochasticCharacterMap( character_histories );

        // get the variables we are interested in
        speciation = bds->getAverageSpeciationRatePerBranch();
        extinction = bds->getAverageExtinctionRatePerBranch();
        delta_speciation = bds->getDeltaSpeciationPerBranch();
        delta_extinction = bds->getDeltaExtinctionPerBranch();
        n_speciation_shifts = bds->getNumberOfSpeciationShiftEventsPerBranch();
        n_extinction_shifts = bds->getNumberOfExtinctionShiftEventsPerBranch();
    }
    
//    for (int i = 0; i < time_in_states.size(); i++)
//    {
//        out_stream << separator;
//        out_stream << i + 1;
//    }
    
    for (int i = 0; i < speciation.size(); i++)
    {
        out_stream << separator;
        out_stream << "avg_lambda[";
        out_stream << i + 1;
        out_stream << "]";
    }
    
    for (int i = 0; i < extinction.size(); i++)
    {
        out_stream << separator;
        out_stream << "avg_mu[";
        out_stream << i + 1;
        out_stream << "]";
    }
   
    if (is_sse | is_glhbdsp)
    {
        for (int i = 0; i < n_shifts.size(); i++)
        {
            out_stream << separator;
            out_stream << "num_shifts[";
            out_stream << i + 1;
            out_stream << "]";
        }
    }
    else
    {
        for (int i = 0; i < delta_speciation.size(); i++)
        {
            out_stream << separator;
            out_stream << "delta_lambda[";
            out_stream << i + 1;
            out_stream << "]";
        }

        for (int i = 0; i < delta_extinction.size(); i++)
        {
            out_stream << separator;
            out_stream << "delta_mu[";
            out_stream << i + 1;
            out_stream << "]";
        }

        
        for (int i = 0; i < n_speciation_shifts.size(); i++)
        {
            out_stream << separator;
            out_stream << "num_speciation_shifts[";
            out_stream << i + 1;
            out_stream << "]";
        }
    
        for (int i = 0; i < n_extinction_shifts.size(); i++)
        {
            out_stream << separator;
            out_stream << "num_extinction_shifts[";
            out_stream << i + 1;
            out_stream << "]";
        }
    }
    
}


void StochasticBranchRateMonitor::swapNode(DagNode *oldN, DagNode* newN)
{
    
    if ( oldN == cdbdp )
    {
        cdbdp = static_cast< StochasticNode<Tree> *>( newN );
    }
    
    VariableMonitor::swapNode( oldN, newN );
    
}


