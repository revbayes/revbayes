#include "BirthDeathShiftMonitor.h"

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
BirthDeathShiftMonitor::BirthDeathShiftMonitor(StochasticNode<Tree>* ch, unsigned long g, const std::string &fname, const std::string &del) : VariableMonitor(ch, g, fname, del, false, false, false),
    bdsp( ch )
{
    // the bdsp is both the tree and character evolution model
    addVariable( bdsp );
}



/**
 * Copy constructor.
 */
BirthDeathShiftMonitor::BirthDeathShiftMonitor( const BirthDeathShiftMonitor &m) : VariableMonitor( m ),
    bdsp( m.bdsp )
{
    
}


/**
 * Destructor.
 */
BirthDeathShiftMonitor::~BirthDeathShiftMonitor()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
BirthDeathShiftMonitor* BirthDeathShiftMonitor::clone(void) const
{
    
    return new BirthDeathShiftMonitor(*this);
}


/**
 * Monitor value at given generation.
 *
 * \param[in]   gen    The current generation.
 */
void BirthDeathShiftMonitor::monitorVariables(unsigned long gen)
{
    auto& separator = to<SeparatorFormat>(format)->separator;
    
	std::vector<double> speciation;
	std::vector<double> extinction;
	std::vector<double> delta_speciation;
	std::vector<double> delta_extinction;
	std::vector<double> sampling;
	std::vector<double> destructive_sampling;
    std::vector<long>   n_speciation_shifts;
    std::vector<long>   n_extinction_shifts;

    size_t num_nodes = bdsp->getValue().getNumberOfNodes();
    std::vector<std::string> character_histories( num_nodes );

    //StateDependentSpeciationExtinctionProcess *sse = dynamic_cast<StateDependentSpeciationExtinctionProcess*>( &bdsp->getDistribution() );
    FastBirthDeathShiftProcess *bds = dynamic_cast<FastBirthDeathShiftProcess*>( &bdsp->getDistribution() );

    // draw stochastic character map
    bds->drawStochasticCharacterMap( character_histories );

    // get the variables we are interested in
    speciation = bds->getAverageSpeciationRatePerBranch();
    extinction = bds->getAverageExtinctionRatePerBranch();
    delta_speciation = bds->getDeltaSpeciationPerBranch();
    delta_extinction = bds->getDeltaExtinctionPerBranch();
    n_speciation_shifts   = bds->getNumberOfSpeciationShiftEventsPerBranch();
    n_extinction_shifts   = bds->getNumberOfExtinctionShiftEventsPerBranch();
    
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


    for (int i = 0; i < delta_speciation.size(); i++)
    {
        // add a separator before every new element
        out_stream << separator;
        out_stream << delta_speciation[i];
    }

    for (int i = 0; i < delta_extinction.size(); i++)
    {
        // add a separator before every new element
        out_stream << separator;
        out_stream << delta_extinction[i];
    }

    
    for (int i = 0; i < n_speciation_shifts.size(); i++)
    {
        // add a separator before every new element
        out_stream << separator;
        out_stream << n_speciation_shifts[i];
    }

    for (int i = 0; i < n_extinction_shifts.size(); i++)
    {
        // add a separator before every new element
        out_stream << separator;
        out_stream << n_extinction_shifts[i];
    }
    
}


/**
 * Print header for monitored values
 */
void BirthDeathShiftMonitor::printFileHeader()
{
    auto& separator = to<SeparatorFormat>(format)->separator;
    
    std::vector<double> speciation;
    std::vector<double> extinction;
    std::vector<double> delta_speciation;
    std::vector<double> delta_extinction;
    std::vector<double> sampling;
    std::vector<double> destructive_sampling;
    std::vector<long>   n_speciation_shifts;
    std::vector<long>   n_extinction_shifts;

    size_t num_nodes = bdsp->getValue().getNumberOfNodes();
    std::vector<std::string> character_histories( num_nodes );

    FastBirthDeathShiftProcess *bds = dynamic_cast<FastBirthDeathShiftProcess*>( &bdsp->getDistribution() );
    // draw stochastic character map
    bds->drawStochasticCharacterMap( character_histories );
    speciation = bds->getAverageSpeciationRatePerBranch();
    extinction = bds->getAverageExtinctionRatePerBranch();
    delta_speciation = bds->getDeltaSpeciationPerBranch();
    delta_extinction = bds->getDeltaExtinctionPerBranch();
    n_speciation_shifts   = bds->getNumberOfSpeciationShiftEventsPerBranch();
    n_extinction_shifts   = bds->getNumberOfExtinctionShiftEventsPerBranch();
        
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


void BirthDeathShiftMonitor::swapNode(DagNode *oldN, DagNode* newN)
{
    
    if ( oldN == bdsp )
    {
        bdsp = static_cast< StochasticNode<Tree> *>( newN );
    }
    
    VariableMonitor::swapNode( oldN, newN );
    
}


