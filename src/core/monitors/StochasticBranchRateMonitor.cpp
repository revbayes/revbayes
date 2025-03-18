#include "StochasticBranchRateMonitor.h"

#include <cstddef>
#include <ostream>
#include <vector>

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
    std::vector<long>   n_shifts;

    size_t num_nodes = cdbdp->getValue().getNumberOfNodes();
    std::vector<std::string> character_histories( num_nodes );

    StateDependentSpeciationExtinctionProcess *sse = dynamic_cast<StateDependentSpeciationExtinctionProcess*>( &cdbdp->getDistribution() );
    if ( sse != NULL )
    {
        // draw stochastic character map
        sse->drawStochasticCharacterMap( character_histories );
        speciation = sse->getAverageSpeciationRatePerBranch();
        extinction = sse->getAverageExtinctionRatePerBranch();
        n_shifts   = sse->getNumberOfShiftEventsPerBranch();
    }
    else
    {
    	GeneralizedLineageHeterogeneousBirthDeathSamplingProcess *glhbdsp = dynamic_cast<GeneralizedLineageHeterogeneousBirthDeathSamplingProcess*>( &cdbdp->getDistribution() );
    	if ( glhbdsp != NULL )
    	{
    		glhbdsp->drawStochasticCharacterMap(character_histories, speciation, extinction, sampling, destructive_sampling, n_shifts);
    	}

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
        // add a separator before every new element
        out_stream << separator;
        
        out_stream << extinction[i];
        
    }
    
    for (int i = 0; i < n_shifts.size(); i++)
    {
        // add a separator before every new element
        out_stream << separator;
        
        out_stream << n_shifts[i];
        
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
    std::vector<long>   n_shifts;

    size_t num_nodes = cdbdp->getValue().getNumberOfNodes();
    std::vector<std::string> character_histories( num_nodes );

    StateDependentSpeciationExtinctionProcess *sse = dynamic_cast<StateDependentSpeciationExtinctionProcess*>( &cdbdp->getDistribution() );
    if ( sse != NULL )
    {
        // draw stochastic character map
        sse->drawStochasticCharacterMap( character_histories );
        speciation = sse->getAverageSpeciationRatePerBranch();
        extinction = sse->getAverageExtinctionRatePerBranch();
        n_shifts   = sse->getNumberOfShiftEventsPerBranch();
    }
    else
    {
    	GeneralizedLineageHeterogeneousBirthDeathSamplingProcess *glhbdsp = dynamic_cast<GeneralizedLineageHeterogeneousBirthDeathSamplingProcess*>( &cdbdp->getDistribution() );
    	if ( glhbdsp != NULL )
    	{
    		glhbdsp->drawStochasticCharacterMap(character_histories, speciation, extinction, sampling, destructive_sampling, n_shifts);
    	}

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
    
    for (int i = 0; i < n_shifts.size(); i++)
    {
        out_stream << separator;
        out_stream << "num_shifts[";
        out_stream << i + 1;
        out_stream << "]";
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


