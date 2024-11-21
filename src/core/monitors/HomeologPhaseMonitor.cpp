#include <map>
#include <cstdint>
#include <ostream>
#include <string>
#include <utility>

#include "HomeologPhaseMonitor.h"
#include "Monitor.h"
#include "AbstractHomologousDiscreteCharacterData.h"
#include "StochasticNode.h"
#include "VariableMonitor.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;


/* Constructor for CTMC */
HomeologPhaseMonitor::HomeologPhaseMonitor(StochasticNode<AbstractHomologousDiscreteCharacterData>* ch, std::uint64_t g, const std::string &fname, const std::string &del) :
    VariableMonitor(ch, g, fname, del, false, false, false),
    ctmc( ch )
{
    
    addVariable( ctmc );
    
}


/**
 * Destructor.
 */
HomeologPhaseMonitor::~HomeologPhaseMonitor()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
HomeologPhaseMonitor* HomeologPhaseMonitor::clone(void) const
{
    
    return new HomeologPhaseMonitor(*this);
}



/**
 * Monitor value at given generation.
 *
 * \param[in]   gen    The current generation.
 */
void HomeologPhaseMonitor::monitorVariables(std::uint64_t gen)
{
    auto& separator = to<SeparatorFormat>(format)->separator;

    std::map<std::string, std::string> homeologMap = ctmc->getValue().getHomeologMap();
    for (std::map<std::string, std::string>::iterator it = homeologMap.begin(); it != homeologMap.end(); ++it)
    {
        out_stream << separator;
        out_stream << it->second;
    }
}


/**
 * Print header for monitored values
 */
void HomeologPhaseMonitor::printFileHeader()
{
    auto& separator = to<SeparatorFormat>(format)->separator;

    std::map<std::string, std::string> homeologMap = ctmc->getValue().getHomeologMap();
    for (std::map<std::string, std::string>::iterator it = homeologMap.begin(); it != homeologMap.end(); ++it)
    {
        out_stream << separator;
        out_stream << it->first;
    }
}


void HomeologPhaseMonitor::swapNode(DagNode *oldN, DagNode* newN)
{
	
	if ( oldN == ctmc )
	{
		ctmc = static_cast< StochasticNode<AbstractHomologousDiscreteCharacterData> *>( newN );
	}
    
    Monitor::swapNode( oldN, newN );
	
}


