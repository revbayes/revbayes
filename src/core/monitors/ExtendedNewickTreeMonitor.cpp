#include "ExtendedNewickTreeMonitor.h"

#include <cstddef>
#include <algorithm>
#include <string>

#include "DagNode.h"
#include "RbException.h"
#include "Cloneable.h"
#include "StringUtilities.h"
#include "TopologyNode.h"
#include "Tree.h"
#include "TypedDagNode.h"

namespace RevBayesCore { template <class valueType> class RbVector; }

using namespace RevBayesCore;


/* Constructor */
ExtendedNewickTreeMonitor::ExtendedNewickTreeMonitor(TypedDagNode<Tree> *t, const std::vector<DagNode*> &n, bool np, unsigned long g, const std::string &fname,
                                                     const std::string &del, bool pp, bool l, bool pr, bool ap) :
    VariableMonitor(t,g,fname,del,pp,l,pr,ap),
    isNodeParameter( np ),
    tree( t ),
    nodeVariables( n )
{
//    this->nodes.insert( tree );

    for (std::vector<DagNode*>::iterator it = nodeVariables.begin(); it != nodeVariables.end(); ++it)
    {
        this->nodes.push_back( *it );

        // tell the node that we have a reference to it (avoids deletion)
        (*it)->incrementReferenceCount();
    }
}


/* Clone the object */
ExtendedNewickTreeMonitor* ExtendedNewickTreeMonitor::clone(void) const
{

    return new ExtendedNewickTreeMonitor(*this);
}


/** Monitor value at generation gen */
void ExtendedNewickTreeMonitor::monitorVariables(unsigned long gen)
{
    auto& separator = to<SeparatorFormat>(format)->separator;

    out_stream << separator;

    tree->getValue().clearParameters();
    for (std::vector<DagNode*>::iterator it = nodeVariables.begin(); it != nodeVariables.end(); ++it)
    {
        DagNode *the_node = *it;
        const std::string &name = the_node->getName();
//        Container *c = dynamic_cast<Container *>( (*it)->getValue() );
        size_t numParams = the_node->getNumberOfElements();

        std::stringstream ss;
        the_node->printValue(ss,"\t", 0, true, false, true);
        std::string concatenatedValues = ss.str();
        std::vector<std::string> values;
        StringUtilities::stringSplit(concatenatedValues, "\t", values);
        for (size_t i = 0; i < numParams; ++i)
        {
            TopologyNode &node = tree->getValue().getNode( i );
            if ( isNodeParameter == true )
            {
                node.addNodeParameter( name, values[i]);
            }
            else
            {
                node.addBranchParameter( name, values[i]);
            }

        }
    }

    out_stream << tree->getValue();
    out_stream.flush();

}


/** Print header for monitored values */
void ExtendedNewickTreeMonitor::printFileHeader()
{
    auto& separator = to<SeparatorFormat>(format)->separator;

    // add a separator tree
    out_stream << separator << "Tree";

}


void ExtendedNewickTreeMonitor::swapNode(DagNode *oldN, DagNode *newN)
{

    TypedDagNode< RbVector<double> >* nodeVar = dynamic_cast< TypedDagNode< RbVector<double> > *>(oldN);
    if ( oldN == tree )
    {
        tree = static_cast< TypedDagNode< Tree > *>( newN );
    }
    else if ( nodeVar != NULL )
    {
        // error catching
//        if ( nodeVariables.find(nodeVar) == nodeVariables.end() )
//        it = find (myvector.begin(), myvector.end(), 30);
//        if (it != myvector.end())
        std::vector<DagNode*>::iterator it = find(nodeVariables.begin(), nodeVariables.end(), nodeVar);
        if (it == nodeVariables.end())
        {
            throw RbException() << "Cannot replace DAG node with name\"" << oldN->getName() << "\" in this extended newick monitor because the monitor doesn't hold this DAG node.";
        }
        *it = static_cast< TypedDagNode< RbVector<double> > *>(newN);
//        nodeVariables.erase( nodeVar );
//        nodeVariables.insert( static_cast< TypedDagNode< RbVector<double> > *>(newN) );
    }

    // delegate to base class
    VariableMonitor::swapNode(oldN, newN);
}
