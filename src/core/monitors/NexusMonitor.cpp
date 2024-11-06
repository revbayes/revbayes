#include "NexusMonitor.h"

#include <cstddef>
#include <ostream>

#include "DagNode.h"
#include "StringUtilities.h"
#include "Taxon.h"
#include "TopologyNode.h"
#include "Tree.h"
#include "TypedDagNode.h"

namespace RevBayesCore {

NexusMonitor::NexusMonitor(TypedDagNode<Tree> *t, const std::vector<DagNode *> &n, bool np, std::uint64_t g,
                           const std::string &fname, bool ap, bool taxa) :
    AbstractFileMonitor (t, g, fname, ap),
    isNodeParameter( np ),
    writeTaxa( taxa ),
    tree( t ),
    nodeVariables( n )
{
    for (size_t i=0; i < nodeVariables.size(); ++i) {
        DagNode* n = nodeVariables[i];
        this->nodes.push_back( n );

        // tell the node that we have a reference to it (avoids deletion)
        n->incrementReferenceCount();
    }
}

NexusMonitor* NexusMonitor::clone() const {
    return new NexusMonitor(*this);
}

void NexusMonitor::swapNode(DagNode *oldN, DagNode *newN) {

    TypedDagNode< RbVector<double> >* nodeVar = dynamic_cast< TypedDagNode< RbVector<double> > *>(oldN);
    if ( oldN == tree ) {
        tree = static_cast< TypedDagNode< Tree > *>( newN );
    }
    else if ( nodeVar != nullptr ) {
        std::vector<DagNode*>::iterator it = find(nodeVariables.begin(), nodeVariables.end(), nodeVar);
        if (it == nodeVariables.end()) {
            throw RbException("Cannot replace DAG node with name\"" + oldN->getName() + "\" in this nexus monitor because the monitor doesn't hold this DAG node.");
        }
        *it = static_cast< TypedDagNode< RbVector<double> > *>(newN);
    }

    // delegate to base class
    AbstractFileMonitor::swapNode(oldN, newN);
}


void NexusMonitor::printHeader() {
    if(!enabled) return;

    out_stream << "#NEXUS" << std::endl;
    out_stream << std::endl;

    if(writeTaxa) {
        std::vector<Taxon> taxa = tree->getValue().getTaxa();
        out_stream << "Begin taxa;" << std::endl;
        out_stream << "\tDimensions ntax=" << taxa.size() << ";" << std::endl;;
        out_stream << "\tTaxlabels" << std::endl;
        for (size_t i=0; i < taxa.size(); ++i) {
            Taxon label = taxa[i];
            out_stream << "\t\t" << label.getName() << std::endl;
        }
        out_stream << "\t\t;" << std::endl;
        out_stream << "End;" << std::endl;
    }

    out_stream << "Begin trees;" << std::endl;
    out_stream.flush();
}

void NexusMonitor::monitor(std::uint64_t gen) {
    if ( !enabled || gen % printgen != 0 ) return;

    out_stream.seekg(0, std::ios::end);
    out_stream << "tree TREE_" << gen << " = " << (tree->getValue().isRooted() ? "[&R]" : "[&U]");

    tree->getValue().clearParameters();
    for (size_t j=0; j < nodeVariables.size(); ++j) {
        DagNode* n = nodeVariables[j];

        const std::string &name = n->getName();
        size_t numParams = n->getNumberOfElements();

        std::stringstream ss;
        n->printValue(ss,"\t", 0, true, false, true);
        std::string concatenatedValues = ss.str();
        std::vector<std::string> values;
        StringUtilities::stringSplit(concatenatedValues, "\t", values);

        for (size_t i = 0; i < numParams; ++i) {
            TopologyNode &node = tree->getValue().getNode( i );
            if ( isNodeParameter == true ) node.addNodeParameter( name, values[i]);
            else node.addBranchParameter( name, values[i]);
        }
    }

    out_stream << tree->getValue() << std::endl;
    out_stream.flush();
}

void NexusMonitor::closeStream() {
    if(enabled) {
        out_stream << "End;" << std::endl;
        out_stream.flush();
    }
    AbstractFileMonitor::closeStream();
}

} /* namespace RevBayesCore */
