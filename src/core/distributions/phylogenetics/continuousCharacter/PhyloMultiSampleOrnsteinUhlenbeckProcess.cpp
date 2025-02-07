#include "PhyloMultiSampleOrnsteinUhlenbeckProcess.h"

#include <cmath>
#include <cstddef>
#include <string>

#include "ConstantNode.h"
#include "DistributionNormal.h"
#include "DistributionMultivariateNormal.h"
#include "RandomNumberFactory.h"
#include "RbException.h"
#include "StochasticNode.h"
#include "TreeUtilities.h"
#include "TopologyNode.h"
#include "Cloneable.h"
#include "ContinuousCharacterData.h"
#include "ContinuousTaxonData.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "Tree.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { class RandomNumberGenerator; }


using namespace RevBayesCore;

PhyloMultiSampleOrnsteinUhlenbeckProcess::PhyloMultiSampleOrnsteinUhlenbeckProcess(const TypedDagNode<Tree> *t, const TypedDagNode< RbVector< double > > *v, const std::vector<Taxon> &ta, size_t ns) : AbstractPhyloContinuousCharacterProcess( t, ns ),
    within_species_variances( v ),
    num_species( t->getValue().getNumberOfTips() ),
    num_individuals( ta.size() ),
    taxa( ta ),
    obs( std::vector<std::vector<double> >(this->num_sites, std::vector<double>(num_individuals, 0.0) ) ),
    means( new std::vector<double>(num_individuals, 0.0) ),
    phylogenetic_covariance_matrix( new MatrixReal(num_individuals, num_individuals) ),
    inverse_phylogenetic_covariance_matrix( num_individuals, num_individuals ),
    changed_covariance(false),
    needs_covariance_recomputation( true ),
    needs_scale_recomputation( true )
{
    // initialize default parameters
    root_state                  = new ConstantNode<double>("", new double(0.0) );
    homogeneous_alpha           = new ConstantNode<double>("", new double(0.0) );
    homogeneous_sigma           = new ConstantNode<double>("", new double(1.0) );
    homogeneous_theta           = new ConstantNode<double>("", new double(0.0) );
    heterogeneous_alpha         = NULL;
    heterogeneous_sigma         = NULL;
    heterogeneous_theta         = NULL;
    
    
    // add parameters
    addParameter( within_species_variances );
    addParameter( homogeneous_alpha );
    addParameter( homogeneous_sigma );
    addParameter( homogeneous_theta );
    
    num_individuals_per_species = std::vector<size_t>(num_species,0);
    const std::vector<TopologyNode *> &nodes = tau->getValue().getNodes();
    for (size_t i=0; i<num_species; ++i)
    {
        const std::string &species_name = nodes[i]->getSpeciesName();
        num_individuals_per_species[i] = getNumberOfSamplesForSpecies( species_name );
    }
    
    inverse_phylogenetic_covariance_matrix.setCholesky( true );
    
    
    // now we need to reset the value
    this->redrawValue();
}

PhyloMultiSampleOrnsteinUhlenbeckProcess::PhyloMultiSampleOrnsteinUhlenbeckProcess(const PhyloMultiSampleOrnsteinUhlenbeckProcess &p) : AbstractPhyloContinuousCharacterProcess( p ),
    root_state( p.root_state ),
    homogeneous_alpha( p.homogeneous_alpha ),
    homogeneous_sigma( p.homogeneous_sigma ),
    homogeneous_theta( p.homogeneous_theta ),
    heterogeneous_alpha( p.heterogeneous_alpha ),
    heterogeneous_sigma( p.heterogeneous_sigma ),
    heterogeneous_theta( p.heterogeneous_theta ),
    within_species_variances( p.within_species_variances ),
    num_species( p.num_species ),
    num_individuals( p.num_individuals ),
    num_individuals_per_species( p.num_individuals_per_species ),
    taxa( p.taxa ),
    obs( p.obs ),
    means( new std::vector<double>( *p.means ) ),
    phylogenetic_covariance_matrix( p.phylogenetic_covariance_matrix->clone() ),
    inverse_phylogenetic_covariance_matrix( p.inverse_phylogenetic_covariance_matrix ),
    changed_covariance( p.changed_covariance ),
    needs_covariance_recomputation( p.needs_covariance_recomputation ),
    needs_scale_recomputation( p.needs_scale_recomputation )
{
    
}



/**
 * Destructor. Because we added ourselves as a reference to tau when we added a listener to its
 * TreeChangeEventHandler, we need to remove ourselves as a reference and possibly delete tau
 * when we die. All other parameters are handled by others.
 */
PhyloMultiSampleOrnsteinUhlenbeckProcess::~PhyloMultiSampleOrnsteinUhlenbeckProcess( void )
{
    
    delete means;
    delete phylogenetic_covariance_matrix;
    
}



PhyloMultiSampleOrnsteinUhlenbeckProcess& PhyloMultiSampleOrnsteinUhlenbeckProcess::operator=(const PhyloMultiSampleOrnsteinUhlenbeckProcess &p)
{
    
    if ( this != &p )
    {
        AbstractPhyloContinuousCharacterProcess::operator=( p );
        
        delete means;
        delete phylogenetic_covariance_matrix;
        
        root_state                              = p.root_state;
        homogeneous_alpha                       = p.homogeneous_alpha;
        homogeneous_sigma                       = p.homogeneous_sigma;
        homogeneous_theta                       = p.homogeneous_theta;
        heterogeneous_alpha                     = p.heterogeneous_alpha;
        heterogeneous_sigma                     = p.heterogeneous_sigma;
        heterogeneous_theta                     = p.heterogeneous_theta;
        within_species_variances                = p.within_species_variances;
        num_species                             = p.num_species;
        obs                                     = p.obs;
        means                                   = new std::vector<double>( *p.means );
        phylogenetic_covariance_matrix          = p.phylogenetic_covariance_matrix->clone();
        inverse_phylogenetic_covariance_matrix  = p.inverse_phylogenetic_covariance_matrix;
        changed_covariance                      = p.changed_covariance;
        needs_covariance_recomputation          = p.needs_covariance_recomputation;
        needs_scale_recomputation               = p.needs_scale_recomputation;
    }
    
    return *this;
}



PhyloMultiSampleOrnsteinUhlenbeckProcess* PhyloMultiSampleOrnsteinUhlenbeckProcess::clone( void ) const
{
    
    return new PhyloMultiSampleOrnsteinUhlenbeckProcess( *this );
}


void PhyloMultiSampleOrnsteinUhlenbeckProcess::computeCovariance(MatrixReal &individual_covariance)
{
    const Tree& tree = tau->getValue();

    // create a temporary matrix for the species covariances
    MatrixReal species_covariance = MatrixReal(num_species,num_species);
    
    const TopologyNode &root = this->tau->getValue().getRoot();
    size_t root_index = root.getIndex();
    
    std::vector<double> variance = std::vector<double>(this->tau->getValue().getNumberOfNodes(), 0.0);
    
    // set the variance at the root
    // we always set this to 0.0 because we expect a root-state parameter
    variance[root_index] = 0.0;
    
    // calculate variance on daughter branches
    for (size_t i=0; i<root.getNumberOfChildren(); ++i)
    {
        const TopologyNode &child = root.getChild( i );
        computeVarianceRecursive(child, variance);
    }
    
    // copy variances into covariance diagonals
    for (size_t i=0; i<num_species; ++i)
    {
        species_covariance[i][i] = variance[i];
    }
    
    // calculate the covariance between all leaf pairs
    for (size_t left_index=0; left_index<num_species; ++left_index)
    {
        const TopologyNode *left_tip_node = &this->tau->getValue().getTipNode( left_index );
        
        for (size_t right_index=(left_index+1); right_index<num_species; ++right_index)
        {
            
            const TopologyNode *right_tip_node = &this->tau->getValue().getTipNode( right_index );
            
            // get mrca
            size_t mrca_index = TreeUtilities::getMrcaIndex(*left_tip_node, *right_tip_node);
            
            // get sum of alpha*branchlength for non-shared branches
            double sum_AT = 0.0;
            
            // first the computation for the left subtree
            size_t current_index = left_index;
            const TopologyNode *current_node = left_tip_node;
            while (current_index != mrca_index)
            {
                sum_AT += computeBranchAlpha(current_index) * tree.getBranchLengthForNode(*current_node);
                current_node = &current_node->getParent();
                current_index = current_node->getIndex();
            }
            
            // first the computation for the left subtree
            current_index = right_index;
            current_node = right_tip_node;
            while (current_index != mrca_index)
            {
                sum_AT += computeBranchAlpha(current_index) * tree.getBranchLengthForNode(*current_node);
                current_node = &current_node->getParent();
                current_index = current_node->getIndex();
            }
            
            species_covariance[left_index][right_index] = variance[mrca_index] * exp(-sum_AT);
        }
        
    }
    
    expandCovariance( individual_covariance, species_covariance );
    
    // make the matrix symmetric and fill remaining half.
    for (size_t i=0; i<num_individuals; ++i)
    {
        
        for (size_t j=0; j<num_individuals; ++j)
        {
            individual_covariance[j][i] = individual_covariance[i][j];
        }
        
    }
    
    
}



void PhyloMultiSampleOrnsteinUhlenbeckProcess::computeExpectation(std::vector<double> &individual_expectations)
{
    // create a temporary matrix for the species covariances
    std::vector<double> species_expectation = std::vector<double>(num_species, 0);
    
    const TopologyNode &root = this->tau->getValue().getRoot();
    
    // set the expectation at the root to the provided root state
    double my_expectation = computeRootState();
    
    // now compute the expectation for all other nodes recursively,
    // starting with the children of the root
    for (size_t i=0; i<root.getNumberOfChildren(); ++i)
    {
        const TopologyNode &child = root.getChild( i );
        computeExpectationRecursive(child, my_expectation, species_expectation);
    }
    
    // expand the expectation to multiple samples per species
    expandExpectation( individual_expectations, species_expectation );
    
}



void PhyloMultiSampleOrnsteinUhlenbeckProcess::computeExpectationRecursive(const TopologyNode &node, double parent_expectation, std::vector<double> &expectations)
{
    // get parameters for this branch
    size_t node_index       = node.getIndex();
    double alpha            = computeBranchAlpha(node_index);
    double theta            = computeBranchTheta(node_index);
    double bl               = tau->getValue().getBranchLengthForNode(node);
    
    double eAT = exp(-1.0 * alpha * bl);
    double my_expectation = 0.0;
    if ( eAT > 1E-10 )
    {
        my_expectation = parent_expectation*eAT + theta*(1.0-eAT);
    }
    else
    {
        my_expectation = parent_expectation;
    }
    
    if ( node.isTip() )
    {
        expectations[node_index] = my_expectation;
    }
    else
    {
        // keep on computing the expectation for all descendant nodes recursively
        for (size_t i=0; i<node.getNumberOfChildren(); ++i)
        {
            const TopologyNode &child = node.getChild( i );
            computeExpectationRecursive(child, my_expectation, expectations);
        }
    }
    
}


void PhyloMultiSampleOrnsteinUhlenbeckProcess::computeVarianceRecursive(const TopologyNode &node, std::vector<double> &variance)
{
    const Tree& tree = tau->getValue();

    // get parameters for this branch
    size_t node_index       = node.getIndex();
    size_t parent_index     = node.getParent().getIndex();
    double alpha            = computeBranchAlpha(node_index);
    double sigma            = computeBranchSigma(node_index);
    double sigma_square     = sigma * sigma;
    double bl               = tree.getBranchLengthForNode(node);
    
    if ( alpha > 1E-10 )
    {
        double e2AT = exp(-2.0 * alpha * bl);
        variance[node_index] = (sigma_square / (2.0*alpha))*(1.0 - e2AT) + variance[parent_index]*e2AT;
    }
    else
    {
        variance[node_index] = sigma_square*bl + variance[parent_index];
    }
    
    for (size_t i=0; i<node.getNumberOfChildren(); ++i)
    {
        const TopologyNode &child = node.getChild( i );
        computeVarianceRecursive(child, variance);
    }
    
}



double PhyloMultiSampleOrnsteinUhlenbeckProcess::computeLnProbability( void )
{
    
    // first, compute the expectations for all tips and the variance-covariance matrix
    computeExpectation( *means );
    computeCovariance( *phylogenetic_covariance_matrix );
    
    inverse_phylogenetic_covariance_matrix = phylogenetic_covariance_matrix->computeInverse();
    
    // we need to make sure that we can use the Cholesky decomposition
    inverse_phylogenetic_covariance_matrix.setCholesky( true );
    
    // sum the probability for each site (column) up
    this->ln_prob = sumRootLikelihood();
    
    return this->ln_prob;
}


double PhyloMultiSampleOrnsteinUhlenbeckProcess::computeBranchAlpha(size_t branch_idx) const
{
    
    // get the selection rate for the branch
    double a;
    if ( this->heterogeneous_alpha != NULL )
    {
        a = this->heterogeneous_alpha->getValue()[branch_idx];
    }
    else
    {
        a = this->homogeneous_alpha->getValue();
    }
    
    return a;
}


double PhyloMultiSampleOrnsteinUhlenbeckProcess::computeBranchSigma(size_t branch_idx) const
{
    
    // get the drift rate for the branch
    double s;
    if ( this->heterogeneous_sigma != NULL )
    {
        s = this->heterogeneous_sigma->getValue()[branch_idx];
    }
    else
    {
        s = this->homogeneous_sigma->getValue();
    }
    
    return s;
}


double PhyloMultiSampleOrnsteinUhlenbeckProcess::computeBranchTheta(size_t branch_idx) const
{
    
    // get the optimum (theta) rate for the branch
    double t;
    if ( this->heterogeneous_theta != NULL )
    {
        t = this->heterogeneous_theta->getValue()[branch_idx];
    }
    else
    {
        t = this->homogeneous_theta->getValue();
    }
    
    return t;
}


double PhyloMultiSampleOrnsteinUhlenbeckProcess::computeRootState( void ) const
{
    
    // get the root-state parameter
    double root_state = this->root_state->getValue();
    
    return root_state;
}


void PhyloMultiSampleOrnsteinUhlenbeckProcess::expandExpectation( std::vector<double> &individual_expectation, const std::vector<double> &species_expectation )
{

    // expand the species to individuals
    size_t index_expectation_individual = 0;
    for (size_t i=0; i<num_species; ++i)
    {

        for (size_t j=0; j<num_individuals_per_species[i]; ++j)
        {
            individual_expectation[index_expectation_individual] = species_expectation[i];
            ++index_expectation_individual;
        }
        
    }
    
}


void PhyloMultiSampleOrnsteinUhlenbeckProcess::expandCovariance( MatrixReal &individual_covariance, const MatrixReal &species_covariance )
{

    // expand the species to individuals
    size_t index_covariance_first_species = 0;
    for (size_t i=0; i<num_species; ++i)
    {
        size_t index_covariance_second_species = index_covariance_first_species;
        
        size_t num_indiv_in_first_pop = num_individuals_per_species[i];
        for (size_t j=i; j<num_species; ++j)
        {
            size_t num_indiv_in_second_pop = num_individuals_per_species[j];
            for (size_t ii=0; ii<num_indiv_in_first_pop; ++ii)
            {
                for (size_t jj=0; jj<num_indiv_in_second_pop; ++jj)
                {
                    individual_covariance[index_covariance_first_species+ii][index_covariance_second_species+jj] = species_covariance[i][j];
                }
            }
            index_covariance_second_species += num_indiv_in_second_pop;
        }
        index_covariance_first_species += num_indiv_in_first_pop;
    }

    // add variance within each species
    index_covariance_first_species = 0;
    const RbVector<double> &variances_per_species = within_species_variances->getValue();
    for (size_t i=0; i<num_species; ++i)
    {
        
        double var_in_species_i = variances_per_species[i];
        for (size_t j=0; j<num_individuals_per_species[i]; ++j)
        {
            individual_covariance[index_covariance_first_species][index_covariance_first_species] += var_in_species_i;
            ++index_covariance_first_species;
        }
        
    }
    
    
}


double PhyloMultiSampleOrnsteinUhlenbeckProcess::getNumberOfSamplesForSpecies(const std::string &name)
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


double PhyloMultiSampleOrnsteinUhlenbeckProcess::getWithinSpeciesVariance(const std::string &name)
{
    
    size_t index = this->tau->getValue().getTipIndex( name );
    
    return within_species_variances->getValue()[ index ];
}




void PhyloMultiSampleOrnsteinUhlenbeckProcess::keepSpecialization( const DagNode* affecter )
{
    
    // reset the flags
    changed_covariance = false;
    needs_covariance_recomputation = false;
    
}


void PhyloMultiSampleOrnsteinUhlenbeckProcess::simulateTipSamples( const std::vector< ContinuousTaxonData > &taxon_data )
{
    
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
                    double stand_dev = species_sigma;
                    
                    // create the character
                    double c = RbStatistics::Normal::rv( parent_state, stand_dev, *rng);
                    
                    // add the character to the sequence
                    individual_data.addCharacter( c );
                }
                
                this->value->addTaxonData( individual_data );
            }
            
        }
        
    }
    
}


void PhyloMultiSampleOrnsteinUhlenbeckProcess::resetValue( void )
{
    
    // create a vector with the correct site indices
    // some of the sites may have been excluded
    std::vector<size_t> site_indices = std::vector<size_t>(this->num_sites,0);
    size_t site_index = 0;
    for (size_t i = 0; i < this->num_sites; ++i)
    {
        while ( this->value->isCharacterExcluded(site_index) )
        {
            site_index++;
            if ( site_index >= this->value->getNumberOfCharacters()  )
            {
                throw RbException( "The character matrix cannot set to this variable because it does not have enough included characters." );
            }
        }
        site_indices[i] = site_index;
        site_index++;
    }
    
    obs = std::vector<std::vector<double> >(this->num_sites, std::vector<double>(num_individuals, 0.0) );
    
    std::vector<TopologyNode*> nodes = this->tau->getValue().getNodes();
    std::vector<size_t> used_indiv_per_species = std::vector<size_t>(num_species,0);
    std::vector<size_t> species_offsets = std::vector<size_t>(num_species,0);
    for (size_t i=1; i<num_species; ++i)
    {
        species_offsets[i] = species_offsets[i-1] + num_individuals_per_species[i-1];
    }

    for (size_t site = 0; site < this->num_sites; ++site)
    {
        
        for (std::vector<Taxon>::iterator it = taxa.begin(); it != taxa.end(); ++it)
        {
            const std::string &species_name = it->getSpeciesName();
            size_t species_index = this->tau->getValue().getTipIndex( species_name );
            size_t index = species_offsets[species_index] + used_indiv_per_species[species_index];
            
            ContinuousTaxonData& taxon = this->value->getTaxonData( it->getName() );
            double &c = taxon.getCharacter(site_indices[site]);
            obs[site][index] = c;
            
            ++used_indiv_per_species[species_index];
        }
    }
    
    // reset the means vectors
    delete means;
    means = new std::vector<double>(num_individuals, 0.0);
    
    
    
    // finally we set all the flags for recomputation
    needs_covariance_recomputation = true;
    needs_scale_recomputation = true;
    
}


void PhyloMultiSampleOrnsteinUhlenbeckProcess::recursiveComputeRootToTipDistance(std::vector<double> &distances, double distance_from_root, const RevBayesCore::TopologyNode &node, size_t node_index)
{
    const Tree& tree = tau->getValue();

    if ( node.isRoot() == false )
    {
        // get my scaled branch length
        double v = this->computeBranchTime(node_index, tree.getBranchLengthForNode(node) );
        
        if ( node.isTip() )
        {
            distances[node_index] = distance_from_root + v;
        }
        else
        {
            const TopologyNode &left = node.getChild(0);
            size_t left_index = left.getIndex();
            recursiveComputeRootToTipDistance(distances, distance_from_root+v, left, left_index );
            
            const TopologyNode &right = node.getChild(1);
            size_t right_index = right.getIndex();
            recursiveComputeRootToTipDistance(distances, distance_from_root+v, right, right_index );
            
        }
        
    }
    else
    {
        
        for (size_t i = 0; i < node.getNumberOfChildren(); ++i)
        {
            const TopologyNode &child = node.getChild(i);
            size_t childIndex = child.getIndex();
            recursiveComputeRootToTipDistance(distances, 0, child, childIndex );
        }
        
    }
    
}


std::set<size_t> PhyloMultiSampleOrnsteinUhlenbeckProcess::recursiveComputeDistanceMatrix(MatrixReal &m, const TopologyNode &node, size_t node_index)
{
    const Tree& tree = tau->getValue();

    // I need to know all my children
    std::set<size_t> children;
    
    if ( node.isRoot() == false )
    {
        // get my scaled branch length
        double v = this->computeBranchTime(node_index, tree.getBranchLengthForNode(node) );
        
        if ( node.isTip() )
        {
            children.insert( node_index );
            m[node_index][node_index] += v;
        }
        else
        {
            const TopologyNode &left = node.getChild(0);
            size_t left_index = left.getIndex();
            children = recursiveComputeDistanceMatrix(m, left, left_index );
            
            const TopologyNode &right = node.getChild(1);
            size_t right_index = right.getIndex();
            std::set<size_t> childrenRight = recursiveComputeDistanceMatrix(m, right, right_index );
            
            children.insert(childrenRight.begin(), childrenRight.end());
            
            // now we loop over all combination of the children pairs to add their variance terms
            for (std::set<size_t>::iterator i_itr = children.begin(); i_itr != children.end(); ++i_itr)
            {
                for (std::set<size_t>::iterator j_itr = children.begin(); j_itr != children.end(); ++j_itr)
                {
                    m[*i_itr][*j_itr] += v;
                }
            }
            
        }
        
    }
    else
    {
        
        for (size_t i = 0; i < node.getNumberOfChildren(); ++i)
        {
            const TopologyNode &child = node.getChild(i);
            size_t childIndex = child.getIndex();
            std::set<size_t> childrenRight = recursiveComputeDistanceMatrix(m, child, childIndex );
        }
        
    }
    
    return children;
    
}



void PhyloMultiSampleOrnsteinUhlenbeckProcess::restoreSpecialization( const DagNode* affecter )
{
    
    // reset the flags
    if ( changed_covariance == true )
    {
        changed_covariance = false;
        needs_covariance_recomputation = false;
        
    }
    
}


void PhyloMultiSampleOrnsteinUhlenbeckProcess::setAlpha(const TypedDagNode<double> *a)
{
    
    // remove the old parameter first
    this->removeParameter( homogeneous_alpha );
    this->removeParameter( heterogeneous_alpha );
    homogeneous_alpha      = NULL;
    heterogeneous_alpha    = NULL;
    
    
    // set the value
    homogeneous_alpha = a;
    
    // add the new parameter
    this->addParameter( homogeneous_alpha );
    
    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }
    
}


void PhyloMultiSampleOrnsteinUhlenbeckProcess::setAlpha(const TypedDagNode<RbVector<double> > *a)
{
    
    // remove the old parameter first
    this->removeParameter( homogeneous_alpha );
    this->removeParameter( heterogeneous_alpha );
    homogeneous_alpha      = NULL;
    heterogeneous_alpha    = NULL;
    
    
    // set the value
    heterogeneous_alpha = a;
    
    // add the new parameter
    this->addParameter( heterogeneous_alpha );
    
    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }
    
}


void PhyloMultiSampleOrnsteinUhlenbeckProcess::setRootState(const TypedDagNode<double> *s)
{
    
    // remove the old parameter first
    this->removeParameter( root_state );
    root_state = s;
    
    // add the new parameter
    this->addParameter( root_state );
    
    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }
    
}


void PhyloMultiSampleOrnsteinUhlenbeckProcess::setSigma(const TypedDagNode<double> *s)
{
    
    // remove the old parameter first
    this->removeParameter( homogeneous_sigma );
    this->removeParameter( heterogeneous_sigma );
    homogeneous_sigma      = NULL;
    heterogeneous_sigma    = NULL;
    
    
    // set the value
    homogeneous_sigma = s;
    
    // add the new parameter
    this->addParameter( homogeneous_sigma );
    
    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }
    
}


void PhyloMultiSampleOrnsteinUhlenbeckProcess::setSigma(const TypedDagNode<RbVector<double> > *s)
{
    
    // remove the old parameter first
    this->removeParameter( homogeneous_sigma );
    this->removeParameter( heterogeneous_sigma );
    homogeneous_sigma      = NULL;
    heterogeneous_sigma    = NULL;
    
    
    // set the value
    heterogeneous_sigma = s;
    
    // add the new parameter
    this->addParameter( heterogeneous_sigma );
    
    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }
    
}


void PhyloMultiSampleOrnsteinUhlenbeckProcess::setTheta(const TypedDagNode<double> *t)
{
    
    // remove the old parameter first
    this->removeParameter( homogeneous_theta );
    this->removeParameter( heterogeneous_theta );
    homogeneous_theta      = NULL;
    heterogeneous_theta    = NULL;
    
    
    // set the value
    homogeneous_theta = t;
    
    // add the new parameter
    this->addParameter( homogeneous_theta );
    
    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }
    
}


void PhyloMultiSampleOrnsteinUhlenbeckProcess::setTheta(const TypedDagNode<RbVector<double> > *t)
{
    
    // remove the old parameter first
    this->removeParameter( homogeneous_theta );
    this->removeParameter( heterogeneous_theta );
    homogeneous_theta      = NULL;
    heterogeneous_theta    = NULL;
    
    
    // set the value
    heterogeneous_theta = t;
    
    // add the new parameter
    this->addParameter( heterogeneous_theta );
    
    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }
    
}


void PhyloMultiSampleOrnsteinUhlenbeckProcess::simulateRecursively( const TopologyNode &node, std::vector< ContinuousTaxonData > &taxa)
{
    const Tree& tree = tau->getValue();

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
        double branch_length = tree.getBranchLengthForNode(child);
        
        // get the branch specific rate
        double branch_time = computeBranchTime( child.getIndex(), branch_length );
        
        // get the branch specific rate
        double branch_sigma = computeBranchSigma( child.getIndex() );
        
        // get the branch specific optimum (theta)
        double branch_theta = computeBranchTheta( child.getIndex() );
        
        // get the branch specific optimum (theta)
        double branch_alpha = computeBranchAlpha( child.getIndex() );
        
        ContinuousTaxonData &taxon = taxa[ child.getIndex() ];
        for ( size_t i = 0; i < num_sites; ++i )
        {
            // get the ancestral character for this site
            double parent_state = parent.getCharacter( i );
            
            // compute the standard deviation for this site
            double branch_rate = branch_time * computeSiteRate( i );
            
            double e = exp(-branch_alpha * branch_rate);
            double e2 = exp(-2 * branch_alpha * branch_rate);
            double m = e * parent_state + (1 - e) * branch_theta;
            
            double stand_dev = 0.0;
            if ( branch_alpha > 1E-10 )
            {
                double sigma_square = branch_sigma * branch_sigma;
                stand_dev = sqrt( (sigma_square / (2.0*branch_alpha)*(1.0 - e2)) );
            }
            else
            {
                // compute the standard deviation for this site
                stand_dev = branch_sigma * sqrt(branch_rate);
            }
            
            // create the character
            double c = RbStatistics::Normal::rv( m, stand_dev, *rng);
            
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


std::vector<double> PhyloMultiSampleOrnsteinUhlenbeckProcess::simulateRootCharacters(size_t n)
{
    
    std::vector<double> chars = std::vector<double>(num_sites, 0);
    for (size_t i=0; i<num_sites; ++i)
    {
        chars[i] = computeRootState();
    }
    
    return chars;
}


double PhyloMultiSampleOrnsteinUhlenbeckProcess::sumRootLikelihood( void )
{
    
    // sum the log-likelihoods for all sites together
    double sum_site_probs = 0.0;
    for (size_t site = 0; site < this->num_sites; ++site)
    {
        sum_site_probs += RbStatistics::MultivariateNormal::lnPdfPrecision( (*means), inverse_phylogenetic_covariance_matrix, obs[site], 1.0);
    }
    
    return sum_site_probs;
}


/** Swap a parameter of the distribution */
void PhyloMultiSampleOrnsteinUhlenbeckProcess::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == root_state)
    {
        root_state = static_cast<const TypedDagNode< double >* >( newP );
    }
    
    if (oldP == homogeneous_alpha)
    {
        homogeneous_alpha = static_cast<const TypedDagNode< double >* >( newP );
    }
    else if (oldP == heterogeneous_alpha)
    {
        heterogeneous_alpha = static_cast<const TypedDagNode< RbVector< double > >* >( newP );
    }
    
    if (oldP == homogeneous_sigma)
    {
        homogeneous_sigma = static_cast<const TypedDagNode< double >* >( newP );
    }
    else if (oldP == heterogeneous_sigma)
    {
        heterogeneous_sigma = static_cast<const TypedDagNode< RbVector< double > >* >( newP );
    }
    
    if (oldP == homogeneous_theta)
    {
        homogeneous_theta = static_cast<const TypedDagNode< double >* >( newP );
    }
    else if (oldP == heterogeneous_theta)
    {
        heterogeneous_theta = static_cast<const TypedDagNode< RbVector< double > >* >( newP );
    }
    
    if (oldP == within_species_variances)
    {
        within_species_variances = static_cast<const TypedDagNode< RbVector< double > >* >( newP );
    }
    
    this->AbstractPhyloContinuousCharacterProcess::swapParameterInternal(oldP, newP);
    
    
    num_individuals_per_species = std::vector<size_t>(num_species,0);
    const std::vector<TopologyNode *> &nodes = tau->getValue().getNodes();
    for (size_t i=0; i<num_species; ++i)
    {
        const std::string &species_name = nodes[i]->getSpeciesName();
        num_individuals_per_species[i] = getNumberOfSamplesForSpecies( species_name );
    }
    
}


void PhyloMultiSampleOrnsteinUhlenbeckProcess::touchSpecialization( const DagNode* affecter, bool touchAll )
{
    
    // if the topology wasn't the culprit for the touch, then we just flag everything as dirty
    if ( affecter == root_state )
    {
        
        
    }
    else
    {
        needs_covariance_recomputation = true;
        if ( changed_covariance == false )
        {
            //            MatrixReal *tmp = phylogenetic_covariance_matrix;
            //            phylogenetic_covariance_matrix = stored_phylogenetic_covariance_matrix;
            //            stored_phylogenetic_covariance_matrix = tmp;
        }
        changed_covariance = true;
        
    }
    //    else if ( affecter != this->tau ) // if the topology wasn't the culprit for the touch, then we just flag everything as dirty
    //    {
    //        touchAll = true;
    //    }
    
}



