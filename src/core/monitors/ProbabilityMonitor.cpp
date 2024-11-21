#include "ProbabilityMonitor.h"

#include <cstddef>
#include <cstdint>
#include <ostream>
#include <set>
#include <string>
#include <vector>

#include "DagNode.h"
#include "Model.h"
#include "Monitor.h"
#include "Cloneable.h"

using namespace RevBayesCore;

/* Constructor */
ProbabilityMonitor::ProbabilityMonitor(std::uint64_t g, const std::string &fname, const std::string &del) : VariableMonitor(std::vector<DagNode *>(),g,fname,del)
{
    
    flatten = false;
}


/**
 * Destructor.
 */
ProbabilityMonitor::~ProbabilityMonitor()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
ProbabilityMonitor* ProbabilityMonitor::clone(void) const
{
    
    return new ProbabilityMonitor(*this);
}


/**
 * Monitor value at generation gen
 */
void ProbabilityMonitor::monitorVariables(std::uint64_t gen)
{
    auto& separator = to<SeparatorFormat>(format)->separator;

    for (std::vector<DagNode*>::iterator i = nodes.begin(); i != nodes.end(); ++i)
    {
        // add a separator before every new element
        out_stream << separator;
        
        // get the node
        DagNode *the_node = *i;
        
        // print the value
        out_stream << the_node->getLnProbability();
    }
    
}



/**
 * Reset the currently monitored DAG nodes by extracting the DAG nodes from the model again
 * and store this in the set of DAG nodes.
 */
void ProbabilityMonitor::resetDagNodes( void )
{
    
    // for savety we empty our dag nodes
    while ( nodes.empty() == false )
    {
        removeVariable( *nodes.begin() );
    }
    
    if ( model != NULL )
    {
        // we only want to have each nodes once
        // this should by default happen by here we check again
        std::set<std::string> var_names;
        
        const std::vector<DagNode*> &n = model->getDagNodes();
        for (std::vector<DagNode*>::const_iterator it = n.begin(); it != n.end(); ++it)
        {
            
            DagNode *the_node = *it;
            
            // only simple numeric variable can be monitored (i.e. only integer and real numbers)
            if ( the_node->isStochastic() == true )
            {
                const std::string &name = the_node->getName();
                if ( var_names.find( name ) == var_names.end() )
                {
                    addVariable( the_node );
                    var_names.insert( name );
                }
                else
                {
#ifdef DEBUG_SEBASTIAN
                    std::cerr << "Trying to add variable with name '" << name << "' twice." << std::endl;
#endif
                }
                
            }
            
        }
        
    }
    
}


/**
 * Set the model from which this monitor will extract the variables.
 * This will automatically result into a reseting of the currently monitored variables.
 *
 * \param[in]   m    The new model.
 */
void ProbabilityMonitor::setModel(Model *m)
{
    
    // delegate call to super class
    Monitor::setModel(m);
    
    // reset the DAG nodes that should be monitored
    resetDagNodes();
    
    sortNodesByName();
}

