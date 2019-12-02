#include "AbstractPhyloContinuousCharacterHistoryProcess.h"
#include "ConstantNode.h"
#include "DistributionNormal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "StochasticNode.h"
#include "TopologyNode.h"

#include <cmath>


using namespace RevBayesCore;


AbstractPhyloContinuousCharacterHistoryProcess::AbstractPhyloContinuousCharacterHistoryProcess(const TypedDagNode<Tree> *t, size_t n) : AbstractPhyloContinuousCharacterProcess( t, n )
{
 
    homogeneous_root_state      = new ConstantNode<double>("", new double(0.0) );
    heterogeneous_root_state    = NULL;
    
    addParameter( homogeneous_root_state );
    
}


/**
 * Destructor.
 */
AbstractPhyloContinuousCharacterHistoryProcess::~AbstractPhyloContinuousCharacterHistoryProcess( void )
{
    // We don't delete the params, because they might be used somewhere else too. The model needs to do that!
    
}

double AbstractPhyloContinuousCharacterHistoryProcess::computeRootState(size_t site_index)
{
    
    double root_state;
    if ( this->heterogeneous_root_state != NULL )
    {
        root_state = this->heterogeneous_root_state->getValue()[site_index];
    }
    else
    {
        root_state = this->homogeneous_root_state->getValue();
    }
    
    return root_state;
    
}

void AbstractPhyloContinuousCharacterHistoryProcess::redrawValue( void )
{
    
    // delete the old value first
    delete this->value;
    
    // create a new character data object
    this->value = new ContinuousCharacterData();
    
    // create a vector of taxon data
    std::vector< ContinuousTaxonData > taxa = std::vector< ContinuousTaxonData >( num_nodes, ContinuousTaxonData( Taxon("") ) );
    
    // simulate the root sequence
    ContinuousTaxonData &root = taxa[ tau->getValue().getRoot().getIndex() ];
    
    std::vector<double> root_states = simulateRootCharacters(num_sites);
    for ( size_t i = 0; i < num_sites; ++i )
    {
        // create the character
        double c = root_states[i];
        
        // add the character to the sequence
        root.addCharacter( c );
    }
    
    // recursively simulate the sequences
    simulateRecursively( tau->getValue().getRoot(), taxa );
    
    // add the taxon data to the character data
    for (size_t i = 0; i < tau->getValue().getNumberOfNodes(); ++i)
    {
        this->value->addTaxonData( taxa[i] );
    }
    
    // tell the derived classes
    this->resetValue();

}

void AbstractPhyloContinuousCharacterHistoryProcess::setRootState(const TypedDagNode<double> *s)
{
    
    // remove the old parameter first
    this->removeParameter( homogeneous_root_state );
    this->removeParameter( heterogeneous_root_state );
    homogeneous_root_state      = NULL;
    heterogeneous_root_state    = NULL;
    
    // set the value
    homogeneous_root_state = s;
    
    // add the new parameter
    this->addParameter( homogeneous_root_state );
    
    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }
    
}


void AbstractPhyloContinuousCharacterHistoryProcess::setRootState(const TypedDagNode<RbVector<double> > *s)
{
    
    // remove the old parameter first
    this->removeParameter( homogeneous_root_state );
    this->removeParameter( heterogeneous_root_state );
    homogeneous_root_state      = NULL;
    heterogeneous_root_state    = NULL;
    
    // set the value
    heterogeneous_root_state = s;
    
    // add the new parameter
    this->addParameter( heterogeneous_root_state );
    
    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }
    
}

std::vector<double> AbstractPhyloContinuousCharacterHistoryProcess::simulateRootCharacters(size_t n)
{
    
    std::vector<double> chars = std::vector<double>(num_sites, 0);
    for (size_t i=0; i<num_sites; ++i)
    {
        chars[i] = computeRootState(i);
    }
    
    return chars;
}

/** Swap a parameter of the distribution */
void AbstractPhyloContinuousCharacterHistoryProcess::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == homogeneous_root_state)
    {
        homogeneous_root_state = static_cast<const TypedDagNode< double >* >( newP );
    }
    else if (oldP == heterogeneous_root_state)
    {
        heterogeneous_root_state = static_cast<const TypedDagNode< RbVector< double > >* >( newP );
    }
    else
    {
        AbstractPhyloContinuousCharacterProcess::swapParameterInternal(oldP, newP);
    }
    
}

