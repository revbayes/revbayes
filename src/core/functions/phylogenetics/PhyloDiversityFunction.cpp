#include "PhyloDiversityFunction.h"

#include <algorithm>
#include <cstddef>
#include <iosfwd>
#include <string>

#include "Clade.h"
#include "TypedFunction.h"
#include "RbBitSet.h"
#include "RbException.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "Taxon.h"
#include "TopologyNode.h"
#include "Tree.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }


using namespace RevBayesCore;


PhyloDiversityFunction::PhyloDiversityFunction(const TypedDagNode<Tree> *t, const Clade &c, const bool i, const TypedDagNode< RbVector< double > > *w) : TypedFunction<double>( new double(0.0) ),
    includeRoot( i ),
    sample( c ),
    tau( t ),
    tipWeights( w )
{
    addParameter( tau );
    this->num_taxa = sample.size();
}


PhyloDiversityFunction::~PhyloDiversityFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


PhyloDiversityFunction* PhyloDiversityFunction::clone( void ) const
{
    return new PhyloDiversityFunction( *this );
}


void PhyloDiversityFunction::update( void )
{
    if (tipWeights == NULL || tipWeights->getNumberOfElements() == 0)
    {
        weightedPD = false;
    }
    else
    {
        if (tipWeights->getNumberOfElements() == sample.size())
        {
            weightedPD = true;
            branchWeights = std::vector<double>(tau->getValue().getNodes().size(), 0.0);
            calculateBranchWeights( tau->getValue().getRoot().getIndex() );
        }
        else
        {
            throw RbException("The sample size must equal the number of weights.");
        }
    }
    
    size_t stopIndex = tau->getValue().getRoot().getIndex();
    
    if (!includeRoot)
    {
        bool found = false;
        stopIndex = 0;
        // get the tmrca of all the members of the sample
        const std::vector<TopologyNode*> &nodes = tau->getValue().getNodes();
        size_t minCladeSize = nodes.size() + 2;
        for (size_t i = tau->getValue().getNumberOfTips(); i < nodes.size(); ++i)
        {
            
            TopologyNode *node = nodes[i];
            std::vector<Taxon> taxa;
            node->getTaxa( taxa );
            size_t cladeSize = taxa.size();
            if ( cladeSize < minCladeSize && cladeSize >= num_taxa && node->containsClade( sample, false ) )
            {
                found = true;
                stopIndex = node->getIndex();
                minCladeSize = cladeSize;
                if ( num_taxa == cladeSize )
                {
                    break;
                }
            }
        }
        if ( found == false )
        {
            *value = 0.0;
            return;
        }
    }
    
    double pd = 0.0;
    std::vector<size_t> nodesVisited;
    std::vector<Taxon> taxa = sample.getTaxa();
    
    // loop through each of the sampled taxa
    for (size_t i = 0; i < num_taxa; ++i)
    {
        std::string taxon = taxa[i].getName();
        size_t j = tau->getValue().getTipIndex( taxon );
        if (j != ((size_t)-1))
        {
            pd += sumPDforNode(j, &nodesVisited, stopIndex);
        }
    }
    
    *value = pd;
}


/**
 * Function that recursively traverses the tree from the given node up to the stopIndex summing PD.
 **/
double PhyloDiversityFunction::sumPDforNode(size_t j, std::vector<size_t>* nodesVisited, size_t stopIndex)
{
    double pd = 0.0;
    // check to see if we have visited this node
    if ( std::find(nodesVisited->begin(), nodesVisited->end(), j) == nodesVisited->end() )
    {
        const TopologyNode& n = tau->getValue().getNode( j );
        if (weightedPD)
        {
            pd += tau->getValue().getBranchLengthForNode(n) * branchWeights[j];
        }
        else
        {
            pd += tau->getValue().getBranchLengthForNode(n);
        }
        nodesVisited->push_back( j );
        if (includeRoot)
        {
            if (j != tau->getValue().getRoot().getIndex())
            {
                size_t parentIndex = n.getParent().getIndex();
                pd += sumPDforNode(parentIndex, nodesVisited, stopIndex);
            }
        }
        else
        {
            if (n.getParent().getIndex() != stopIndex)
            {
                size_t parentIndex = n.getParent().getIndex();
                pd += sumPDforNode(parentIndex, nodesVisited, stopIndex);
            }
        }
    }
    return pd;
}


/**
 * Function that recursively traverses the tree in preorder calculating branch weights.
 **/
double PhyloDiversityFunction::calculateBranchWeights(size_t j)
{
    double weight = 0.0;
    const TopologyNode& n = tau->getValue().getNode( j );
    if (n.isTip())
    {
        std::vector<Taxon> taxa = sample.getTaxa();
        for (size_t i = 0; i < num_taxa; ++i)
        {
            if ( n.getName() == taxa[i].getName() )
            {
                branchWeights[j] = tipWeights->getValue()[i];
                weight = branchWeights[j];
                break;
            }
        }
    }
    else
    {
        std::vector<int> children = n.getChildrenIndices();
        double weightLeft = calculateBranchWeights( children[0] );
        double weightRight = calculateBranchWeights( children[1] );
        branchWeights[j] = weightLeft + weightRight / 2;
        weight = branchWeights[j];
    }
    return weight;
}


void PhyloDiversityFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == tau)
    {
        tau = static_cast<const TypedDagNode<Tree>* >( newP );
    }
}

