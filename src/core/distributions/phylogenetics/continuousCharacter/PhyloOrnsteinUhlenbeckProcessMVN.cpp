#include "PhyloOrnsteinUhlenbeckProcessMVN.h"

#include <cmath>
#include <cstddef>
#include <string>
#include <iosfwd>

#include "ConstantNode.h"
#include "DistributionNormal.h"
#include "DistributionMultivariateNormal.h"
#include "RandomNumberFactory.h"
#include "RbException.h"
#include "StochasticNode.h"
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

PhyloOrnsteinUhlenbeckProcessMVN::PhyloOrnsteinUhlenbeckProcessMVN(const TypedDagNode<Tree> *t, size_t ns) :
    AbstractPhyloContinuousCharacterProcess( t, ns ),
    num_tips( t->getValue().getNumberOfTips() ),
    obs( std::vector<std::vector<double> >(this->num_sites, std::vector<double>(num_tips, 0.0) ) ),
    means( new std::vector<std::vector<double> >(this->num_sites, std::vector<double>(num_tips, 0.0) ) ),
    stored_means( new std::vector<std::vector<double> >(this->num_sites, std::vector<double>(num_tips, 0.0) ) ),
    phylogenetic_covariance_matrix( new MatrixReal(num_tips, num_tips) ),
    stored_phylogenetic_covariance_matrix( new MatrixReal(num_tips, num_tips) ),
    inverse_phylogenetic_covariance_matrix( num_tips, num_tips ),
    changed_covariance(false),
    needs_covariance_recomputation( true ),
    needs_scale_recomputation( true )
{
    // initialize default parameters
    homogeneous_root_state      = new ConstantNode<double>("", new double(0.0) );
    homogeneous_alpha           = new ConstantNode<double>("", new double(1.0) );
    homogeneous_sigma           = new ConstantNode<double>("", new double(1.0) );
    homogeneous_theta           = new ConstantNode<double>("", new double(0.0) );
    heterogeneous_root_state    = NULL;
    heterogeneous_alpha         = NULL;
    heterogeneous_sigma         = NULL;
    heterogeneous_theta         = NULL;
    

    // add parameters
    addParameter( homogeneous_root_state );
    addParameter( homogeneous_alpha );
    addParameter( homogeneous_sigma );
    addParameter( homogeneous_theta );
    
    
    // now we need to reset the value
    this->redrawValue();
}


/**
 * Destructor. Because we added ourselves as a reference to tau when we added a listener to its
 * TreeChangeEventHandler, we need to remove ourselves as a reference and possibly delete tau
 * when we die. All other parameters are handled by others.
 */
PhyloOrnsteinUhlenbeckProcessMVN::~PhyloOrnsteinUhlenbeckProcessMVN( void )
{

//    delete phylogenetic_covariance_matrix;
//    delete stored_phylogenetic_covariance_matrix;
}



PhyloOrnsteinUhlenbeckProcessMVN* PhyloOrnsteinUhlenbeckProcessMVN::clone( void ) const
{
    
    return new PhyloOrnsteinUhlenbeckProcessMVN( *this );
}


double PhyloOrnsteinUhlenbeckProcessMVN::computeLnProbability( void )
{
    
    // compute the ln probability by recursively calling the probability calculation for each node
    const TopologyNode &root = this->tau->getValue().getRoot();
    
    // we start with the root and then traverse down the tree
    size_t root_index = root.getIndex();
    
    if ( needs_covariance_recomputation == true )
    {
        std::vector<double> distances = std::vector<double>(num_tips,0.0);
        recursiveComputeRootToTipDistance(distances, 0.0, root, root_index);
        MatrixReal c_matrix = MatrixReal(num_tips,num_tips);
        for (size_t i=0; i<num_tips; ++i)
        {
            c_matrix[i][i] = distances[i] + distances[i];
            for (size_t j=i+1; j<num_tips; ++j)
            {
                c_matrix[i][j] = distances[i] + distances[j];
                c_matrix[j][i] = distances[i] + distances[j];
            }
        }
        
        MatrixReal shared_distances_matrix = MatrixReal(num_tips, num_tips);
        recursiveComputeDistanceMatrix(shared_distances_matrix, root, root_index);
        
        double sel   = computeBranchAlpha(0);
        double drift = computeBranchSigma(0);
        double opt   = computeBranchTheta(0);
        double var = (drift*drift)/sel * 0.5;
        
        MatrixReal &m = *phylogenetic_covariance_matrix;
        for (size_t i=0; i<num_tips; ++i)
        {
            for (size_t j=0; j<num_tips; ++j)
            {
                m[i][j] = var * exp(-sel*c_matrix[i][j]) * (exp(2*sel*shared_distances_matrix[i][j])-1.0);
            }
        }
        inverse_phylogenetic_covariance_matrix = phylogenetic_covariance_matrix->computeInverse();

        // now compute the means
        for (size_t i=0; i<num_sites; ++i)
        {
            for (size_t j=0; j<num_tips; ++j)
            {
                (*means)[i][j] = opt * (1.0-exp(-sel*distances[j])) + computeRootState(i)*exp(-sel*distances[j]);
            }
        }
        
        needs_covariance_recomputation = false;
    }
    
    // sum the partials up
    this->ln_prob = sumRootLikelihood();
    
    return this->ln_prob;
}


double PhyloOrnsteinUhlenbeckProcessMVN::computeBranchAlpha(size_t branch_idx) const
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


double PhyloOrnsteinUhlenbeckProcessMVN::computeBranchSigma(size_t branch_idx) const
{
    
    // get the selection rate for the branch
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


double PhyloOrnsteinUhlenbeckProcessMVN::computeBranchTheta(size_t branch_idx) const
{
    
    // get the selection rate for the branch
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


double PhyloOrnsteinUhlenbeckProcessMVN::computeRootState(size_t siteIdx) const
{
    
    // second, get the clock rate for the branch
    double root_state;
    if ( this->heterogeneous_root_state != NULL )
    {
        root_state = this->heterogeneous_root_state->getValue()[siteIdx];
    }
    else
    {
        root_state = this->homogeneous_root_state->getValue();
    }
    
    return root_state;
}




void PhyloOrnsteinUhlenbeckProcessMVN::keepSpecialization( const DagNode* affecter )
{
    
    // reset the flags
    changed_covariance = false;
    needs_covariance_recomputation = false;
    
}


void PhyloOrnsteinUhlenbeckProcessMVN::resetValue( void )
{
    
    // create a vector with the correct site indices
    // some of the sites may have been excluded
    std::vector<size_t> siteIndices = std::vector<size_t>(this->num_sites,0);
    size_t siteIndex = 0;
    for (size_t i = 0; i < this->num_sites; ++i)
    {
        while ( this->value->isCharacterExcluded(siteIndex) )
        {
            siteIndex++;
            if ( siteIndex >= this->value->getNumberOfCharacters()  )
            {
                throw RbException( "The character matrix cannot set to this variable because it does not have enough included characters." );
            }
        }
        siteIndices[i] = siteIndex;
        siteIndex++;
    }
    
    obs = std::vector<std::vector<double> >(this->num_sites, std::vector<double>(num_tips, 0.0) );
    
    std::vector<TopologyNode*> nodes = this->tau->getValue().getNodes();
    for (size_t site = 0; site < this->num_sites; ++site)
    {
        
        for (std::vector<TopologyNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it)
        {
            if ( (*it)->isTip() )
            {
                ContinuousTaxonData& taxon = this->value->getTaxonData( (*it)->getName() );
                double &c = taxon.getCharacter(siteIndices[site]);
                obs[site][(*it)->getIndex()] = c;
            }
        }
    }
    
    // reset the means vectors
    delete means;
    delete stored_means;
    means           = new std::vector<std::vector<double> >(this->num_sites, std::vector<double>(num_tips, 0.0) );
    stored_means    = new std::vector<std::vector<double> >(this->num_sites, std::vector<double>(num_tips, 0.0) );
    
    
    
    // finally we set all the flags for recomputation
    needs_covariance_recomputation = true;
    needs_scale_recomputation = true;
    
}


void PhyloOrnsteinUhlenbeckProcessMVN::recursiveComputeRootToTipDistance(std::vector<double> &distances, double distance_from_root, const RevBayesCore::TopologyNode &node, size_t node_index)
{
    
    if ( node.isRoot() == false )
    {
        // get my scaled branch length
        double v = this->computeBranchTime(node_index, tau->getValue().getBranchLengthForNode(node));
        
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


std::set<size_t> PhyloOrnsteinUhlenbeckProcessMVN::recursiveComputeDistanceMatrix(MatrixReal &m, const TopologyNode &node, size_t node_index)
{
    
    // I need to know all my children
    std::set<size_t> children;
    
    if ( node.isRoot() == false )
    {
        // get my scaled branch length
        double v = this->computeBranchTime(node_index, tau->getValue().getBranchLengthForNode(node));
        
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



void PhyloOrnsteinUhlenbeckProcessMVN::restoreSpecialization( const DagNode* affecter )
{
    
    // reset the flags
    if ( changed_covariance == true )
    {
        changed_covariance = false;
        needs_covariance_recomputation = false;
        
        MatrixReal *tmp = phylogenetic_covariance_matrix;
        phylogenetic_covariance_matrix = stored_phylogenetic_covariance_matrix;
        stored_phylogenetic_covariance_matrix = tmp;
        
    }
    
}


void PhyloOrnsteinUhlenbeckProcessMVN::setAlpha(const TypedDagNode<double> *a)
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


void PhyloOrnsteinUhlenbeckProcessMVN::setAlpha(const TypedDagNode<RbVector<double> > *a)
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


void PhyloOrnsteinUhlenbeckProcessMVN::setRootState(const TypedDagNode<double> *s)
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


void PhyloOrnsteinUhlenbeckProcessMVN::setRootState(const TypedDagNode<RbVector<double> > *s)
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


void PhyloOrnsteinUhlenbeckProcessMVN::setSigma(const TypedDagNode<double> *s)
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


void PhyloOrnsteinUhlenbeckProcessMVN::setSigma(const TypedDagNode<RbVector<double> > *s)
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


void PhyloOrnsteinUhlenbeckProcessMVN::setTheta(const TypedDagNode<double> *t)
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


void PhyloOrnsteinUhlenbeckProcessMVN::setTheta(const TypedDagNode<RbVector<double> > *t)
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


void PhyloOrnsteinUhlenbeckProcessMVN::simulateRecursively( const TopologyNode &node, std::vector< ContinuousTaxonData > &taxa)
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
        double branch_time = sqrt( computeBranchTime( child.getIndex(), branch_length ) );
        
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
            double standDev = branch_sigma * sqrt((1 - e2) / 2 / branch_alpha);
            
            // create the character
            double c = RbStatistics::Normal::rv( m, standDev, *rng);
            
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


std::vector<double> PhyloOrnsteinUhlenbeckProcessMVN::simulateRootCharacters(size_t n)
{
    
    std::vector<double> chars = std::vector<double>(num_sites, 0);
    for (size_t i=0; i<num_sites; ++i)
    {
        chars[i] = computeRootState(i);
    }
    
    return chars;
}


double PhyloOrnsteinUhlenbeckProcessMVN::sumRootLikelihood( void )
{
    
    // sum the log-likelihoods for all sites together
    double sumPartialProbs = 0.0;
    for (size_t site = 0; site < this->num_sites; ++site)
    {
        double sigma = this->computeSiteRate(site);
        //        sumPartialProbs += RbStatistics::MultivariateNormal::lnPdfCovariance(m, *phylogeneticCovarianceMatrix, obs[site], sigma*sigma);
        sumPartialProbs += RbStatistics::MultivariateNormal::lnPdfPrecision( (*means)[site], inverse_phylogenetic_covariance_matrix, obs[site], sigma*sigma);
    }
    
    return sumPartialProbs;
}


/** Swap a parameter of the distribution */
void PhyloOrnsteinUhlenbeckProcessMVN::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == homogeneous_root_state)
    {
        homogeneous_root_state = static_cast<const TypedDagNode< double >* >( newP );
    }
    else if (oldP == heterogeneous_root_state)
    {
        heterogeneous_root_state = static_cast<const TypedDagNode< RbVector< double > >* >( newP );
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

    this->AbstractPhyloContinuousCharacterProcess::swapParameterInternal(oldP, newP);
    
}


void PhyloOrnsteinUhlenbeckProcessMVN::touchSpecialization( const DagNode* affecter, bool touchAll )
{
    
    // if the topology wasn't the culprit for the touch, then we just flag everything as dirty
    if ( affecter == homogeneous_root_state )
    {
        
        
    }
    else if ( affecter == heterogeneous_root_state )
    {
        
        
    }
    else
    {
        needs_covariance_recomputation = true;
        if ( changed_covariance == false )
        {
            MatrixReal *tmp = phylogenetic_covariance_matrix;
            phylogenetic_covariance_matrix = stored_phylogenetic_covariance_matrix;
            stored_phylogenetic_covariance_matrix = tmp;
        }
        changed_covariance = true;
        
    }
//    else if ( affecter != this->tau ) // if the topology wasn't the culprit for the touch, then we just flag everything as dirty
//    {
//        touchAll = true;
//    }
    
}


