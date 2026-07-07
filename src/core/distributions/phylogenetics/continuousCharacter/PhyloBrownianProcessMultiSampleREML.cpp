#include <cmath>
#include <cstddef>
#include <iosfwd>
#include <set>
#include <string>
#include <vector>

#include "DistributionNormal.h"
#include "PhyloBrownianProcessMultiSampleREML.h"
#include "RandomNumberFactory.h"
#include "RbException.h"
#include "StochasticNode.h"
#include "TopologyNode.h"
#include "AbstractPhyloBrownianProcess.h"
#include "ContinuousCharacterData.h"
#include "ContinuousTaxonData.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "Taxon.h"
#include "Tree.h"
#include "TreeChangeEventHandler.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { class RandomNumberGenerator; }


using namespace RevBayesCore;

PhyloBrownianProcessMultiSampleREML::PhyloBrownianProcessMultiSampleREML(const TypedDagNode<Tree> *tr, const TypedDagNode< RbVector< double > > *v, const std::vector<Taxon> &ta, size_t ns) : AbstractPhyloBrownianProcess( tr, ns ),
    within_species_variances( v ),
    node_likelihoods(this->num_nodes),
    taxa( ta )
{
    
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    this->addParameter( within_species_variances );
    
    // We don'e want tau to die before we die, or it can't remove us as listener
    tau->getValue().getTreeChangeEventHandler().addListener( this );
    
    // now we need to reset the value
    this->redrawValue();
    
    // we need to reset the means
    resetValue();
}


/**
 * Destructor. Because we added ourselves as a reference to tau when we added a listener to its
 * TreeChangeEventHandler, we need to remove ourselves as a reference and possibly delete tau
 * when we die. All other parameters are handled by others.
 */
PhyloBrownianProcessMultiSampleREML::~PhyloBrownianProcessMultiSampleREML( void )
{
    // We don't delete the params, because they might be used somewhere else too. The model needs to do that!
    
    // remove myself from the tree listeners
    if ( tau != NULL )
    {
        tau->getValue().getTreeChangeEventHandler().removeListener( this );
    }
    
}



PhyloBrownianProcessMultiSampleREML* PhyloBrownianProcessMultiSampleREML::clone( void ) const
{
    
    return new PhyloBrownianProcessMultiSampleREML( *this );
}


double PhyloBrownianProcessMultiSampleREML::computeLnProbability( void )
{
    // we need to check here if we still are listining to this tree for change events
    // the tree could have been replaced without telling us
    if ( tau->getValue().getTreeChangeEventHandler().isListening( this ) == false )
    {
        tau->getValue().getTreeChangeEventHandler().addListener( this );
        resetValue();
    }
    
    // compute the ln probability by recursively calling the probability calculation for each node
    const TopologyNode &root = this->tau->getValue().getRoot();
    
    // we start with the root and then traverse down the tree
    size_t root_index = root.getIndex();
    
    // only necessary if the root is actually dirty
    if ( not this->node_likelihoods.is_valid(root_index) )
    {
        recursiveComputeLnProbability( root, root_index );
    }
    
    // sum the partials up
    this->ln_prob = sumRootLikelihood();

    return this->ln_prob;
}


double PhyloBrownianProcessMultiSampleREML::computeMeanForSpecies(const std::string &name, size_t index)
{
    
    double mean = 0.0;
    double num_samples = 0.0;
    
    for (size_t i=0; i<taxa.size(); ++i)
    {
        
        const Taxon &t = taxa[i];
        if ( name == t.getSpeciesName() )
        {
            ContinuousTaxonData& taxon = this->value->getTaxonData( t.getName() );
            
            if ( taxon.isCharacterResolved(index) == true )
            {
                mean += taxon.getCharacter(index);

                ++num_samples;
            }
            
        }
        
    }
    
    // normalize
    if ( num_samples > 0 )
    {
        mean /= num_samples;
    }
    else
    {
        mean = RbConstants::Double::nan;
    }
    
    return mean;
}


void PhyloBrownianProcessMultiSampleREML::fireTreeChangeEvent( const TopologyNode &n, const unsigned& m )
{
    
    // call a recursive flagging of all node above (closer to the root) and including this node
    recursivelyFlagNodeDirty( n );
    
}


double PhyloBrownianProcessMultiSampleREML::getNumberOfSamplesForSpecies(const std::string &name)
{
    
    double num_samples = 0.0;
    
    for (size_t i=0; i<taxa.size(); ++i)
    {
        
        const Taxon &t = taxa[i];
        if ( name == t.getSpeciesName() )
        {
            ++num_samples;
        }
        
    }
    
    return num_samples;
}


double PhyloBrownianProcessMultiSampleREML::getWithinSpeciesVariance(const std::string &name)
{
    
    size_t index = this->tau->getValue().getTipIndex( name );
    
    return within_species_variances->getValue()[ index ];
}




void PhyloBrownianProcessMultiSampleREML::keepSpecialization( const DagNode* affecter )
{
    if (node_likelihoods.has_snapshot())
    {
        node_likelihoods.keep();
    }
}


void PhyloBrownianProcessMultiSampleREML::recursiveComputeLnProbability( const TopologyNode &node, size_t node_index )
{
    // check for recomputation
    if ( node.isTip() == true && not node_likelihoods.is_valid(node_index) )
    {
        NodeCache &node_cache = this->node_likelihoods.init_for_writing(node_index);
        node_cache.partial_likelihoods.assign(this->num_sites, 0.0);
        node_cache.means.assign(this->num_sites, 0.0);
        node_cache.variances.assign(this->num_sites, 0.0);
        node_cache.missing_data.assign(this->num_sites, false);

        std::vector<double> &mu_node  = node_cache.means;
        std::vector<double> &var_node = node_cache.variances;
        std::vector<double> &p_node   = node_cache.partial_likelihoods;

        const std::string &name = this->tau->getValue().getNode( node_index ).getName();

        std::vector<double> num_samples(this->num_sites, 0.0);
        
        double within_species_variance = getWithinSpeciesVariance(name);

        // initialize log probabilities to zero
        for (int i=0; i<this->num_sites; i++)
        {
             p_node[i] = 0;
        }

        // first we grab the specimen data
        std::vector<std::vector<double> > observations(this->num_sites);

        std::vector<Taxon> specimens; 

        for (Taxon &specimen : taxa)
        {
            if (specimen.getSpeciesName() == name)
            {
                string specimen_name = specimen.getName();
                ContinuousTaxonData& dt = this->value->getTaxonData( specimen_name );

                for (int char_index=0; char_index<this->num_sites; ++char_index)
                {
                    if ( dt.isCharacterResolved( site_indices[char_index] ) )
                    {
                        double x = dt.getCharacter( site_indices[char_index] );

                        if (std::isfinite(x) )
                        {
                            observations[char_index].push_back(x); 
                            ++num_samples[char_index];
                        }
                    }
                }
            }
        }

        // next we calculate the probability of the contrasts
        // and we calculate what mu and v should be for the species node
        //
        for (int char_index=0; char_index<this->num_sites; ++char_index)
        {
            node_cache.missing_data[char_index] = observations[char_index].empty();

            if ( !observations[char_index].empty() ){
                double mu = observations[char_index][0];
                double var = within_species_variance;

                for (size_t n = 1; n < num_samples[char_index]; n++){
                    double mu_left = mu;
                    double var_left = var;

                    double mu_right = observations[char_index][n];
                    double var_right = within_species_variance;

                    double contrast = mu_left - mu_right;
                    double var_contrast = var_left + var_right;

                    double lnl = RbStatistics::Normal::lnPdf(0, std::sqrt(var_contrast), contrast);

                    if ( RbMath::isFinite(lnl) == false )
                    {
                        std::cerr << "Issue in computing contrast probability in PhyloBrownianProcessMultiSampleREML." << std::endl;
                    }

                    p_node[char_index] += lnl;

                    mu = (mu_left * var_right + mu_right * var_left) / ( var_left + var_right );
                    var = (var_left * var_right) / (var_left + var_right);
                }

                // store the mean and variance variables for the tree recursion
                mu_node[char_index] = mu; 
                var_node[char_index] = var; 
            }
            else 
            {
                mu_node[char_index] = RbConstants::Double::nan;
                var_node[char_index] = 0.0;
            }
        }
    }
    else if ( node.isTip() == false && not node_likelihoods.is_valid(node_index) )
    {
        NodeCache &node_cache = this->node_likelihoods.init_for_writing(node_index);
        node_cache.partial_likelihoods.assign(this->num_sites, 0.0);
        node_cache.means.assign(this->num_sites, 0.0);
        node_cache.variances.assign(this->num_sites, 0.0);
        node_cache.missing_data.assign(this->num_sites, false);

        std::vector<double> &mu_node = node_cache.means;
        std::vector<double> &v_node  = node_cache.variances;
        std::vector<double> &p_node  = node_cache.partial_likelihoods;

        // get the number of children
        size_t num_children = node.getNumberOfChildren();
        if (num_children != 2 )
        {
            throw RbException("internal node in the phylogeny does not have two descendants (in PhyloBrownianProcessMultiSampleREML), not supported");
        }

        const TopologyNode &left = node.getChild(0);
        size_t left_index = left.getIndex();
        recursiveComputeLnProbability( left, left_index );

        const TopologyNode &right = node.getChild(1);
        size_t right_index = right.getIndex();
        recursiveComputeLnProbability( right, right_index );

        const NodeCache &left_cache  = this->node_likelihoods[left_index];
        const NodeCache &right_cache = this->node_likelihoods[right_index];
        const std::vector<double> &mu_left  = left_cache.means;
        const std::vector<double> &mu_right = right_cache.means;

        const std::vector<double> &v_left   = left_cache.variances;
        const std::vector<double> &v_right  = right_cache.variances;
        
        const std::vector<double> &p_left   = left_cache.partial_likelihoods;
        const std::vector<double> &p_right  = right_cache.partial_likelihoods;
       
        size_t num_sites = this->num_sites;

        for (size_t char_index = 0; char_index < num_sites; char_index++)
        {
            bool left_missing = left_cache.missing_data[char_index];
            bool right_missing = right_cache.missing_data[char_index];

            if ( left_missing && right_missing )
            {
                node_cache.missing_data[char_index] = true;
                mu_node[char_index] = RbConstants::Double::nan;
                v_node[char_index] = 0.0;
                p_node[char_index] = 0.0;
            }
            else if ( left_missing && !right_missing )
            {
                node_cache.missing_data[char_index] = false;
                mu_node[char_index] = mu_right[char_index];
                double var_right = this->computeBranchTime(right_index, right.getBranchLength());
                v_node[char_index] = v_right[char_index] + var_right;
                p_node[char_index] = p_right[char_index];
            }
            else if ( !left_missing && right_missing )
            {
                node_cache.missing_data[char_index] = false;
                mu_node[char_index] = mu_left[char_index];
                double var_left = this->computeBranchTime(left_index, left.getBranchLength());
                v_node[char_index] = v_left[char_index] + var_left;
                p_node[char_index] = p_left[char_index];
            }
            else 
            {
                node_cache.missing_data[char_index] = false;

                // merging rule
                // D_node(y) = D_left(y) * D_right(y)
                double var_left = this->computeBranchTime(left_index, left.getBranchLength()) + v_left[char_index];
                double var_right = this->computeBranchTime(right_index, right.getBranchLength()) + v_right[char_index];

                mu_node[char_index] = (mu_left[char_index]*var_right + mu_right[char_index]*var_left) / (var_left+var_right);
                v_node[char_index] = (var_left*var_right) / (var_left+var_right);

                double contrast = mu_left[char_index] - mu_right[char_index];

               
                double stdev = sqrt(var_left + var_right);
                double p = RbStatistics::Normal::lnPdf(0, stdev, contrast);

                p_node[char_index] = p_left[char_index] + p_right[char_index] + p;
            }
        }
       
    } // end if we need to compute something for this node.
}


void PhyloBrownianProcessMultiSampleREML::recursivelyFlagNodeDirty( const TopologyNode &n )
{
    
    // we need to flag this node and all ancestral nodes for recomputation
    size_t index = n.getIndex();
    
    // if this node is already invalid, then all ancestral nodes must have been invalidated too
    if ( node_likelihoods.is_valid(index) )
    {
        // the root doesn't have an ancestor
        if ( n.isRoot() == false )
        {
            recursivelyFlagNodeDirty( n.getParent() );
        }

        node_likelihoods.invalidate(index);
    }
    
}


/*
 * Invalidate the computation that uses a branch transform.
 * Branch parameters belong to the edge ending at n, and are applied by n's parent.
 */
void PhyloBrownianProcessMultiSampleREML::invalidateBranchAndAncestors( const TopologyNode &n )
{
    if (n.isRoot())
    {
        recursivelyFlagNodeDirty(n);
    }
    else
    {
        recursivelyFlagNodeDirty(n.getParent());
    }
}



void PhyloBrownianProcessMultiSampleREML::redrawValue( void )
{
    
    // delete the old value first
    delete this->value;
    
    // create a new character data object
    this->value = new ContinuousCharacterData();
    
    // create a vector of taxon data
    std::vector< ContinuousTaxonData > taxon_data = std::vector< ContinuousTaxonData >( num_nodes, ContinuousTaxonData( Taxon("") ) );
    
    // simulate the root sequence
    ContinuousTaxonData &root = taxon_data[ tau->getValue().getRoot().getIndex() ];
    
    std::vector<double> root_states = simulateRootCharacters(num_sites);
    for ( size_t i = 0; i < num_sites; ++i )
    {
        // create the character
        double c = root_states[i];
        
        // add the character to the sequence
        root.addCharacter( c );
    }
    
    // recursively simulate the sequences
    simulateRecursively( tau->getValue().getRoot(), taxon_data );
    
    // Get the random number generator
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    // add the taxon data to the character data
    for (size_t i = 0; i < tau->getValue().getNumberOfTips(); ++i)
    {
        const std::string &species_name = tau->getValue().getNode(i).getName();
        const ContinuousTaxonData &species_data = taxon_data[i];
        double species_sigma = sqrt( getWithinSpeciesVariance( species_name ) );
        
        for ( size_t j=0; j<taxa.size(); ++j )
        {
            
            const Taxon &t = taxa[j];
            if ( species_name == t.getSpeciesName() )
            {
                ContinuousTaxonData individual_data = ContinuousTaxonData( t );
                
                for ( size_t k = 0; k < num_sites; ++k )
                {
                    
                    // get the ancestral character for this site
                    double parent_state = species_data.getCharacter(k);
                    
                    // compute the standard deviation for this site
                    double stand_dev = species_sigma * computeSiteRate( k );
                    
                    // create the character
                    double c = RbStatistics::Normal::rv( parent_state, stand_dev, *rng);

                    // add the character to the sequence
                    individual_data.addCharacter( c );
                }
                
                this->value->addTaxonData( individual_data );
            }
            
        }
    }
    
    // tell the derived classes
    this->resetValue();
}


void PhyloBrownianProcessMultiSampleREML::resetValue( void )
{
    
    this->num_nodes = tau->getValue().getNumberOfNodes();
    node_likelihoods.resize(this->num_nodes);

    // create a vector with the correct site indices
    // some of the sites may have been excluded
    site_indices = std::vector<size_t>(this->num_sites,0);
    size_t site_index = 0;
    for (size_t i = 0; i < this->num_sites; ++i)
    {
        while ( this->value->isCharacterExcluded(site_index) )
        {
            ++site_index;
            if ( site_index >= this->value->getNumberOfCharacters()  )
            {
                throw RbException( "The character matrix cannot set to this variable because it does not have enough included characters." );
            }
        }
        site_indices[i] = site_index;
        ++site_index;
    }
    
    // first we check for missing data
    std::vector<std::vector<bool> > tip_missing(this->num_nodes, std::vector<bool>(this->num_sites, true));
    for (size_t site = 0; site < this->num_sites; ++site)
    {
        for (size_t i=0; i<taxa.size(); ++i)
        {
            const Taxon &t = taxa[i];
            const ContinuousTaxonData& taxon = this->value->getTaxonData( t.getName() );
            double c = taxon.getCharacter(site_indices[site]);
            
            if ( taxon.isCharacterResolved(site_indices[site]) == true && RbMath::isFinite(c) == true )
            {
                size_t species_index = tau->getValue().getTipIndex( t.getSpeciesName() );
                tip_missing[species_index][site] = false;
            }
        }
    }
    use_missing_data = false;
    std::vector<TopologyNode*> nodes = this->tau->getValue().getNodes();
    for (size_t site = 0; site < this->num_sites; ++site)
    {
        
        for (std::vector<TopologyNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it)
        {
            if ( (*it)->isTip() )
            {
                size_t species_index = (*it)->getIndex();
                
                if ( tip_missing[species_index][site] )
                {
                    use_missing_data = true;
                    break;
                }
                
            }
        }
    }
    
}


void PhyloBrownianProcessMultiSampleREML::restoreSpecialization( const DagNode* affecter )
{
    if (node_likelihoods.has_snapshot())
    {
        node_likelihoods.restore();
    }
}


std::vector<double> PhyloBrownianProcessMultiSampleREML::simulateRootCharacters(size_t n)
{
    
    std::vector<double> chars = std::vector<double>(num_sites, 0);
    for (size_t i=0; i<num_sites; ++i)
    {
        chars[i] = 0.0;
    }
    
    return chars;
}


double PhyloBrownianProcessMultiSampleREML::sumRootLikelihood( void )
{
    // get the root node
    const TopologyNode &root = this->tau->getValue().getRoot();
    
    // get the index of the root node
    size_t node_index = root.getIndex();
    
    // get the pointers to the partial likelihoods of the left and right subtree
    const std::vector<double> &p_node = this->node_likelihoods[node_index].partial_likelihoods;
    
    // sum the log-likelihoods for all sites together
    double sum_partial_probs = 0.0;
    for (size_t site = 0; site < this->num_sites; ++site)
    {
        sum_partial_probs += p_node[site];
    }
    
    return sum_partial_probs;
}


void PhyloBrownianProcessMultiSampleREML::snapshotSpecialization( void )
{
    node_likelihoods.snapshot();
}


/*
 * Mark multisample Brownian REML likelihood caches dirty after a dependency changes.
 * Snapshot state is stored by IndexedSnapshotCache before this invalidation hook runs.
 */
void PhyloBrownianProcessMultiSampleREML::invalidateSpecialization( const DagNode* affecter, bool touchAll )
{
    
    // if the topology wasn't the culprit for the touch, then we just flag everything as dirty
    if ( affecter == this->heterogeneous_clock_rates )
    {
        
        const std::set<size_t> &indices = this->heterogeneous_clock_rates->getTouchedElementIndices();
        
        // maybe all of them have been touched or the flags haven't been set properly
        if ( indices.size() == 0 )
        {
            // just flag everyting for recomputation
            touchAll = true;
        }
        else
        {
            const std::vector<TopologyNode *> &nodes = this->tau->getValue().getNodes();
            // flag recomputation only for the nodes
            for (std::set<size_t>::iterator it = indices.begin(); it != indices.end(); ++it)
            {
                if ( *it >= nodes.size() )
                {
                    touchAll = true;
                    break;
                }

                this->invalidateBranchAndAncestors( *nodes[*it] );
            }
        }
    }

    if ( affecter == this->within_species_variances )
    {
        
        const std::set<size_t> &indices = this->within_species_variances->getTouchedElementIndices();
        
        // maybe all of them have been touched or the flags haven't been set properly
        if ( indices.size() == 0 )
        {
            // just flag everyting for recomputation
            touchAll = true;
        }
        else
        {
            const std::vector<TopologyNode *> &nodes = this->tau->getValue().getNodes();
            // The variance vector is indexed by species/tip order; fall back to full invalidation
            // if that mapping is not available for a touched index.
            for (std::set<size_t>::iterator it = indices.begin(); it != indices.end(); ++it)
            {
                if ( *it >= nodes.size() || nodes[*it]->isTip() == false )
                {
                    touchAll = true;
                    break;
                }

                this->recursivelyFlagNodeDirty(*nodes[*it]);
            }
        }
        
    }

    if ( affecter == this->dag_node )
    {
        resetValue();
    }
    else if ( affecter != this->heterogeneous_clock_rates && affecter != this->within_species_variances && affecter != this->tau ) // if the topology wasn't the culprit for the touch, then we just flag everything as dirty
    {
        touchAll = true;
    }
    
    if ( touchAll )
    {
        node_likelihoods.invalidate_all();
    }
    
}


/** Swap a parameter of the distribution */
void PhyloBrownianProcessMultiSampleREML::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == this->tau)
    {
        this->tau->getValue().getTreeChangeEventHandler().removeListener( this );
        AbstractPhyloBrownianProcess::swapParameterInternal(oldP, newP);
        this->tau->getValue().getTreeChangeEventHandler().addListener( this );
    }
    else if ( oldP == this->within_species_variances )
    {
        within_species_variances = static_cast<const TypedDagNode< RbVector< double > >* >( newP );
    }
    else
    {
        AbstractPhyloBrownianProcess::swapParameterInternal(oldP, newP);
    }
    
}
