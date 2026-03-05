#include "AbstractPhyloContinuousCharacterProcess.h"

#include <string>

#include "ConstantNode.h"
#include "StochasticNode.h"
#include "TopologyNode.h"
#include "Cloneable.h"
#include "ContinuousTaxonData.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "Taxon.h"
#include "Tree.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }


using namespace RevBayesCore;


AbstractPhyloContinuousCharacterProcess::AbstractPhyloContinuousCharacterProcess(const TypedDagNode<Tree> *t, size_t ns) : TypedDistribution< ContinuousCharacterData >(  new ContinuousCharacterData() ),
    num_nodes( t->getValue().getNumberOfNodes() ),
    num_sites( ns ),
    tau( t ),
    use_missing_data(false)
{
    
    // initialize with default parameters
    homogeneous_clock_rate      = new ConstantNode<double>("", new double(1.0) );
    heterogeneous_clock_rates   = NULL;
    homogeneous_site_rate       = new ConstantNode<double>("", new double(1.0) );
    heterogeneous_site_rates    = NULL;
    
    
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    this->addParameter( tau );
    this->addParameter( homogeneous_clock_rate );
    this->addParameter( homogeneous_site_rate );
    
}


/**
 * Destructor.
 */
AbstractPhyloContinuousCharacterProcess::~AbstractPhyloContinuousCharacterProcess( void )
{
    // We don't delete the params, because they might be used somewhere else too. The model needs to do that!
    
    
    // remove myself from the tree listeners
//    if ( this->tau != NULL )
//    {
//        
//        this->tau->getValue().getTreeChangeEventHandler().removeListener( this );
//        
//    }
}


double AbstractPhyloContinuousCharacterProcess::computeBranchTime( size_t nide_idx, double brlen )
{
    
    // get the clock rate for the branch
    double branch_time;
    if ( this->heterogeneous_clock_rates != NULL )
    {
        double sigma = this->heterogeneous_clock_rates->getValue()[nide_idx];
        branch_time = sigma * sigma * brlen;
    }
    else
    {
        double sigma = this->homogeneous_clock_rate->getValue();
        branch_time = sigma * sigma * brlen;
    }
    
    return branch_time;
}


double AbstractPhyloContinuousCharacterProcess::computeSiteRate(size_t siteIdx)
{
    
    // second, get the clock rate for the branch
    double sr;
    if ( this->heterogeneous_site_rates != NULL )
    {
        sr = this->heterogeneous_site_rates->getValue()[siteIdx];
    }
    else
    {
        sr = this->homogeneous_site_rate->getValue();
    }
    
    return sr;
}



void AbstractPhyloContinuousCharacterProcess::redrawValue( void )
{
    
    // delete the old value first
    delete this->value;
    
    // create a new character data object
    this->value = new ContinuousCharacterData();
    
    // create a vector of taxon data
    std::vector< ContinuousTaxonData > taxa = std::vector< ContinuousTaxonData >( num_nodes, ContinuousTaxonData( Taxon("") ) );
    
    // simulate the root sequence
    size_t root_index = tau->getValue().getRoot().getIndex();
    ContinuousTaxonData &root = taxa[ root_index ];
    
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
    
    // we call now our method to resample the tips
    // this is important if we have multiple samples (e.g. individuals) per species
    simulateTipSamples( taxa );
    
    // tell the derived classes
    this->resetValue();
}


void AbstractPhyloContinuousCharacterProcess::reInitialized( void )
{
    
    // tell the derived classes
    this->resetValue();
}


void AbstractPhyloContinuousCharacterProcess::resetValue( void )
{
    // nothing to do here.
    // derived class may overwrite this method
}


void AbstractPhyloContinuousCharacterProcess::setBranchRate(const TypedDagNode<double> *r)
{
    
    // remove the old parameter first
    this->removeParameter( homogeneous_clock_rate );
    this->removeParameter( heterogeneous_clock_rates );
    homogeneous_clock_rate      = NULL;
    heterogeneous_clock_rates   = NULL;
    
    
    // set the value
    homogeneous_clock_rate = r;
    
    // add the new parameter
    this->addParameter( homogeneous_clock_rate );
    
    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }
    
}


void AbstractPhyloContinuousCharacterProcess::setBranchRate(const TypedDagNode<RbVector<double> > *r)
{
    
    // remove the old parameter first
    this->removeParameter( homogeneous_clock_rate );
    this->removeParameter( heterogeneous_clock_rates );
    homogeneous_clock_rate      = NULL;
    heterogeneous_clock_rates   = NULL;
    
    
    // set the value
    heterogeneous_clock_rates = r;
    
    // add the new parameter
    this->addParameter( heterogeneous_clock_rates );
    
    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }
    
}


void AbstractPhyloContinuousCharacterProcess::setSiteRate(const TypedDagNode<double> *r)
{
    
    // remove the old parameter first
    this->removeParameter( homogeneous_site_rate );
    this->removeParameter( heterogeneous_site_rates );
    homogeneous_site_rate      = NULL;
    heterogeneous_site_rates   = NULL;
    
    
    // set the value
    homogeneous_site_rate = r;
    
    // add the new parameter
    this->addParameter( homogeneous_site_rate );
    
    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }
    
}


void AbstractPhyloContinuousCharacterProcess::setSiteRate(const TypedDagNode<RbVector<double> > *r)
{
    
    // remove the old parameter first
    this->removeParameter( homogeneous_site_rate );
    this->removeParameter( heterogeneous_site_rates );
    homogeneous_site_rate       = NULL;
    heterogeneous_site_rates    = NULL;
    
    this->num_sites = r->getValue().size();
    
    
    // set the value
    heterogeneous_site_rates = r;
    
    // add the new parameter
    this->addParameter( heterogeneous_site_rates );
    
    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }
    
}


void AbstractPhyloContinuousCharacterProcess::setValue(ContinuousCharacterData *v, bool force)
{
    
    // delegate to the parent class
    TypedDistribution< ContinuousCharacterData >::setValue(v, force);
    
    // reset the number of sites
    num_sites = v->getNumberOfIncludedCharacters();
    
    // tell the derived classes
    this->resetValue();
    
}


void AbstractPhyloContinuousCharacterProcess::simulateTipSamples( const std::vector< ContinuousTaxonData > &taxon_data )
{
    
    // add the taxon data to the character data
    for (size_t i = 0; i < tau->getValue().getNumberOfTips(); ++i)
    {
        this->value->addTaxonData( taxon_data[i] );
    }
    
}


/** Swap a parameter of the distribution */
void AbstractPhyloContinuousCharacterProcess::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == homogeneous_clock_rate)
    {
        homogeneous_clock_rate = static_cast<const TypedDagNode< double >* >( newP );
    }
    else if (oldP == heterogeneous_clock_rates)
    {
        heterogeneous_clock_rates = static_cast<const TypedDagNode< RbVector< double > >* >( newP );
    }
    else if (oldP == homogeneous_site_rate)
    {
        homogeneous_site_rate = static_cast<const TypedDagNode< double >* >( newP );
    }
    else if (oldP == heterogeneous_site_rates)
    {
        heterogeneous_site_rates = static_cast<const TypedDagNode< RbVector< double > >* >( newP );
    }
    else if (oldP == tau)
    {
        tau = static_cast<const TypedDagNode<Tree>* >( newP );
        num_nodes = tau->getValue().getNumberOfNodes();
    }
    
}


