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
        // We only want to have each node once. This should happen by default, but here we check again.
        std::set<std::string> added_var_names;
        std::set<std::string> vector_base_names;
        
        const std::vector<DagNode*> &n = model->getDagNodes();
        
        // Step 1: collect base names of vector variables
        for (DagNode* the_node : n)
        {
            // only simple numeric variables can be monitored (i.e. only integers and real numbers)
            if ( the_node->isSimpleNumeric() && the_node->isElementVariable() )
            {
                const std::string &name = the_node->getName();
                
                // extract vector name from the name of a vector element: everything before the last '['
                size_t bracket_pos = name.rfind('[');
                if (bracket_pos != std::string::npos)
                {
                    std::string base_name = name.substr(0, bracket_pos);
                    vector_base_names.insert(base_name);
                }
            }
        }
        
        // Step 2: add variables, including vector elements but not the vectors themselves
        for (DagNode* the_node : n)
        {
            if ( the_node->isSimpleNumeric() && !the_node->isClamped() )
            {
                const std::string &name = the_node->getName();
                
                // skip this node if it is a vector whose elements we are already collecting
                if ( vector_base_names.count(name) > 0 )
                {
                    continue;  // skip the non-element version
                }
                
                bool condition = !stochastic_nodes_only && !the_node->isConstant() && name != "";
                
                if ( !the_node->isHidden() && ( condition || the_node->isStochastic() ) )
                {
                    if ( exclude.find(name) == exclude.end() && added_var_names.find(name) == added_var_names.end() )
                    {
                        addVariable( the_node );
                        added_var_names.insert( name );
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
    
    sortNodesByName(true);
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

