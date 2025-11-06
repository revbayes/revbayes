#include "StochasticBranchStateTimesMonitor.h"

#include <cstddef>
#include <vector>

#include "StochasticNode.h"
#include "StateDependentSpeciationExtinctionProcess.h"
#include "Cloneable.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;


/* Constructor for state dependent birth death process */
StochasticBranchStateTimesMonitor::StochasticBranchStateTimesMonitor(StochasticNode<Tree>* ch, std::uint64_t g, const std::string &fname, const std::string &del) : VariableMonitor(ch, g, fname, del, false, false, false),
    cdbdp( ch )
{
    // the cdbdp is both the tree and character evolution model
    addVariable( cdbdp );
    tree = static_cast< StochasticNode<Tree> *>( cdbdp );
}



/**
 * Copy constructor.
 */
StochasticBranchStateTimesMonitor::StochasticBranchStateTimesMonitor( const StochasticBranchStateTimesMonitor &m) : VariableMonitor( m ),
    tree( m.tree ),
    cdbdp( m.cdbdp )
{
    
}


/**
 * Destructor.
 */
StochasticBranchStateTimesMonitor::~StochasticBranchStateTimesMonitor()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
StochasticBranchStateTimesMonitor* StochasticBranchStateTimesMonitor::clone(void) const
{
    
    return new StochasticBranchStateTimesMonitor(*this);
}


/**
 * Monitor value at given generation.
 *
 * \param[in]   gen    The current generation.
 */
void StochasticBranchStateTimesMonitor::monitorVariables(std::uint64_t gen)
{
    auto& separator = to<SeparatorFormat>(format)->separator;

    StateDependentSpeciationExtinctionProcess *sse = dynamic_cast<StateDependentSpeciationExtinctionProcess*>( &cdbdp->getDistribution() );
    size_t num_nodes = tree->getValue().getNumberOfNodes();
    std::vector<std::string> character_histories( num_nodes );
    
    // draw stochastic character map
    sse->drawStochasticCharacterMap( character_histories );
    std::vector<double> time_in_states = sse->getTimeInStates(); 

    // print to monitor file
    for (int i = 0; i < time_in_states.size(); i++)
    {
        // add a separator before every new element
        out_stream << separator;
            
        out_stream << time_in_states[i];
        
    }
    
}


/**
 * Print header for monitored values
 */
void StochasticBranchStateTimesMonitor::printFileHeader()
{
    auto& separator = to<SeparatorFormat>(format)->separator;

    StateDependentSpeciationExtinctionProcess *sse = dynamic_cast<StateDependentSpeciationExtinctionProcess*>( &cdbdp->getDistribution() );
    std::vector<double> time_in_states = sse->getTimeInStates(); 
    
    for (int i = 0; i < time_in_states.size(); i++)
    {
        out_stream << separator;
        out_stream << i + 1;
    }
    
}


void StochasticBranchStateTimesMonitor::swapNode(DagNode *oldN, DagNode* newN)
{
    
    if ( oldN == tree )
    {
        tree = static_cast< TypedDagNode<Tree> *>( newN );
    }
    else if ( oldN == cdbdp )
    {
        cdbdp = static_cast< StochasticNode<Tree> *>( newN );
        tree = static_cast< StochasticNode<Tree> *>( newN );
    }
    
    VariableMonitor::swapNode( oldN, newN );
    
}


