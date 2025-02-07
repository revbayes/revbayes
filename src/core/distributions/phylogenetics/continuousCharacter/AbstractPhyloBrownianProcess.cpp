#include "AbstractPhyloBrownianProcess.h"

#include <cmath>

#include "DistributionNormal.h"
#include "RandomNumberFactory.h"
#include "TopologyNode.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "Cloneable.h"
#include "ContinuousTaxonData.h"

namespace RevBayesCore { class RandomNumberGenerator; }
namespace RevBayesCore { class Tree; }
namespace RevBayesCore { template <class valueType> class TypedDagNode; }


using namespace RevBayesCore;


AbstractPhyloBrownianProcess::AbstractPhyloBrownianProcess(const TypedDagNode<Tree> *t, size_t n) : AbstractPhyloContinuousCharacterProcess( t, n )
{
    
    // now we need to reset the value
    // Sebastian: We cannot do this here because redraw has calls the pure virtual functions simulateRootCharacters and simulateRecursively
//    this->redrawValue();
}


/**
 * Destructor.
 */
AbstractPhyloBrownianProcess::~AbstractPhyloBrownianProcess( void )
{
    // We don't delete the params, because they might be used somewhere else too. The model needs to do that!
    
}


void AbstractPhyloBrownianProcess::simulateRecursively( const TopologyNode &node, std::vector< ContinuousTaxonData > &taxa)
{
    
    // get the children of the node
    const std::vector<TopologyNode*>& children = node.getChildren();
    
    // get the sequence of this node
    size_t node_index = node.getIndex();
    const ContinuousTaxonData &parent = taxa[ node_index ];
    
    // simulate the sequence for each child
    RandomNumberGenerator* rng = GLOBAL_RNG;
    for (std::vector< TopologyNode* >::const_iterator it = children.begin(); it != children.end(); ++it)
    {
        const TopologyNode &child = *(*it);
        
        // get the branch length for this child
        double branch_length = tau->getValue().getBranchLengthForNode(child);
        
        // get the branch specific rate
        double branch_sigma = sqrt( computeBranchTime( child.getIndex(), branch_length ) );
        
        ContinuousTaxonData &taxon = taxa[ child.getIndex() ];
        for ( size_t i = 0; i < num_sites; ++i )
        {
            // get the ancestral character for this site
            double parent_state = parent.getCharacter( i );
            
            // compute the standard deviation for this site
            double standDev = branch_sigma * computeSiteRate( i );
            
            // create the character
            double c = RbStatistics::Normal::rv( parent_state, standDev, *rng);
            
            // add the character to the sequence
            taxon.addCharacter( c );
        }
        
        if ( child.isTip() )
        {
            taxon.setTaxon( child.getTaxon() );
        }
        else
        {
            // recursively simulate the sequences
            simulateRecursively( child, taxa );
        }
        
    }
    
}
