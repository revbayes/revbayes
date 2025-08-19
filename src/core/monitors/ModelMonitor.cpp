#include "ModelMonitor.h"

#include <cstddef>
#include <set>
#include <string>
#include <vector>

#include "DagNode.h"
#include "Model.h"
#include "Monitor.h"
#include "Cloneable.h"

using namespace RevBayesCore;

/* Constructor */
ModelMonitor::ModelMonitor(std::uint64_t g, const path &fname, const SampleFormat &f, std::set<std::string> exclude_list) :
    VariableMonitor(std::vector<DagNode *>(),g,fname,f),
    stochastic_nodes_only( false ),
    exclude(exclude_list)
{}


/**
 * Destructor.
 */
ModelMonitor::~ModelMonitor()
{}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself 
 */
ModelMonitor* ModelMonitor::clone(void) const 
{
    
    return new ModelMonitor(*this);
}



/**
 * Reset the currently monitored DAG nodes by extracting the DAG nodes from the model again 
 * and store this in the set of DAG nodes.
 */
void ModelMonitor::resetDagNodes( void )
{
    
    // for safety we empty our dag nodes
    while ( nodes.empty() == false )
    {
        removeVariable( *nodes.begin() );
    }
    
    if ( model != NULL )
    {
        // we only want to have each nodes once
        // this should by default happen but here we check again
        std::set<std::string> var_names;
        
        const std::vector<DagNode*> &n = model->getDagNodes();
        for (std::vector<DagNode*>::const_iterator it = n.begin(); it != n.end(); ++it) 
        {

            DagNode *the_node = *it;

            // only simple numeric variables can be monitored (i.e. only integer and real numbers)
            if ( the_node->isSimpleNumeric() && the_node->isClamped() == false )
            {
                if ( (!stochastic_nodes_only && !the_node->isConstant() && the_node->getName() != "" && !the_node->isHidden() && !the_node->isElementVariable() ) ||
                     ( the_node->isStochastic() && !the_node->isClamped() && the_node->isHidden() == false  && the_node->isElementVariable() == false ) )
                {
                    const std::string &name = the_node->getName();
                    if ( exclude.find(name) == exclude.end() && var_names.find( name ) == var_names.end() )
                    {
                        addVariable( the_node );
                        var_names.insert( name );
                    }
                    
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
void ModelMonitor::setModel(Model *m)
{
    
    // delegate call to super class
    Monitor::setModel(m);
    
    // reset the DAG nodes that should be monitored
    resetDagNodes();
    
    sortNodesByName();
}


/**
 * Set flag about whether to monitor only stochastic DAG nodes or deterministic ones as well. 
 *
 * \param[in]   tf   Flag if only stochastic nodes should be monitored.
 */
void ModelMonitor::setStochasticNodesOnly(bool tf) 
{
    
    stochastic_nodes_only = tf;
    
    // reset the DAG nodes that should be monitored
    resetDagNodes();
    
}

