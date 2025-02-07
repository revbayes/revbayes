#include <cstddef>
#include <iosfwd>
#include <set>
#include <string>
#include <vector>

#include "ConstantNode.h"
#include "PhyloBrownianProcessMVN.h"
#include "DistributionMultivariateNormal.h"
#include "RbException.h"
#include "StochasticNode.h"
#include "TopologyNode.h"
#include "AbstractPhyloBrownianProcess.h"
#include "ContinuousCharacterData.h"
#include "ContinuousTaxonData.h"
#include "MatrixReal.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "StringUtilities.h"
#include "Tree.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }


using namespace RevBayesCore;

PhyloBrownianProcessMVN::PhyloBrownianProcessMVN(const TypedDagNode<Tree> *t, size_t ns) : AbstractPhyloBrownianProcess( t, ns ),
    num_tips( t->getValue().getNumberOfTips() ),
    obs( std::vector<std::vector<double> >(this->num_sites, std::vector<double>(num_tips, 0.0) ) ),
    phylogenetic_covariance_matrix( new MatrixReal(num_tips, num_tips) ),
    stored_phylogenetic_covariance_matrix( new MatrixReal(num_tips, num_tips) ),
    changed_covariance(false),
    needs_covariance_recomputation( true ),
    needs_scale_recomputation( true )
{
    homogeneous_root_state      = new ConstantNode<double>("", new double(0.0) );
    heterogeneous_root_state    = NULL;

    addParameter( homogeneous_root_state );
    phylogenetic_covariance_matrix->setCholesky( true );
    stored_phylogenetic_covariance_matrix->setCholesky( true );
    
    // now we need to reset the value
    this->redrawValue();
}


PhyloBrownianProcessMVN::PhyloBrownianProcessMVN(const PhyloBrownianProcessMVN &p) : AbstractPhyloBrownianProcess( p ),
    homogeneous_root_state( p.homogeneous_root_state ),
    heterogeneous_root_state( p.heterogeneous_root_state ),
    num_tips( p.num_tips ),
    obs( p.obs ),
    phylogenetic_covariance_matrix( p.phylogenetic_covariance_matrix->clone() ),
    stored_phylogenetic_covariance_matrix( p.stored_phylogenetic_covariance_matrix->clone() ),
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
PhyloBrownianProcessMVN::~PhyloBrownianProcessMVN( void )
{
    // We don't delete the params, because they might be used somewhere else too. The model needs to do that!
    delete phylogenetic_covariance_matrix;
    delete stored_phylogenetic_covariance_matrix;
    
}



PhyloBrownianProcessMVN& PhyloBrownianProcessMVN::operator=(const PhyloBrownianProcessMVN &p)

{
    
    if ( this != &p )
    {
        AbstractPhyloBrownianProcess::operator=( p );
        
        delete phylogenetic_covariance_matrix;
        delete stored_phylogenetic_covariance_matrix;
        
        homogeneous_root_state                      = p.homogeneous_root_state;
        heterogeneous_root_state                    = p.heterogeneous_root_state;
        num_tips                                    = p.num_tips;
        obs                                         = p.obs;
        phylogenetic_covariance_matrix              = p.phylogenetic_covariance_matrix->clone();
        stored_phylogenetic_covariance_matrix       = p.stored_phylogenetic_covariance_matrix->clone();
        inverse_phylogenetic_covariance_matrix      = p.inverse_phylogenetic_covariance_matrix;
        changed_covariance                          = p.changed_covariance;
        needs_covariance_recomputation              = p.needs_covariance_recomputation;
        needs_scale_recomputation                   = p.needs_scale_recomputation;
    }
    
    return *this;
}



PhyloBrownianProcessMVN* PhyloBrownianProcessMVN::clone( void ) const
{
    
    return new PhyloBrownianProcessMVN( *this );
}


double PhyloBrownianProcessMVN::computeLnProbability( void )
{
    
    // compute the ln probability by recursively calling the probability calculation for each node
    const TopologyNode &root = this->tau->getValue().getRoot();
    
    // we start with the root and then traverse down the tree
    size_t rootIndex = root.getIndex();
    
    if ( needs_covariance_recomputation == true )
    {
        // perhaps there is a more efficient way to reset the matrix to 0.
        delete phylogenetic_covariance_matrix;
        phylogenetic_covariance_matrix = new MatrixReal(num_tips, num_tips);
        phylogenetic_covariance_matrix->setCholesky( true );
        recursiveComputeCovarianceMatrix(*phylogenetic_covariance_matrix, root, rootIndex);
        needs_covariance_recomputation = false;
    }
    
    // sum the partials up
    this->ln_prob = sumRootLikelihood();
    
    return this->ln_prob;
}


double PhyloBrownianProcessMVN::computeRootState(size_t siteIdx)
{
    
    // second, get the clock rate for the branch
    double rootState;
    if ( this->heterogeneous_root_state != NULL )
    {
        rootState = this->heterogeneous_root_state->getValue()[siteIdx];
    }
    else
    {
        rootState = this->homogeneous_root_state->getValue();
    }
    
    return rootState;
}




void PhyloBrownianProcessMVN::keepSpecialization( const DagNode* affecter )
{
    
    // reset the flags
    changed_covariance = false;
    needs_covariance_recomputation = false;
    
}


void PhyloBrownianProcessMVN::resetValue( void )
{
    
    // create a vector with the correct site indices
    // some of the sites may have been excluded
    std::vector<size_t> site_indices = std::vector<size_t>(this->num_sites,0);
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
    
    obs = std::vector<std::vector<double> >(this->num_sites, std::vector<double>(num_tips, 0.0) );
    
    std::vector<TopologyNode*> nodes = this->tau->getValue().getNodes();
    for (size_t site = 0; site < this->num_sites; ++site)
    {
        
        for (std::vector<TopologyNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it)
        {
            if ( (*it)->isTip() )
            {
                ContinuousTaxonData& taxon = this->value->getTaxonData( (*it)->getName() );
                double &c = taxon.getCharacter(site_indices[site]);
                obs[site][(*it)->getIndex()] = c;
            }
        }
    }
    
    
    // finally we set all the flags for recomputation
    needs_covariance_recomputation = true;
    needs_scale_recomputation = true;
    
}


std::set<size_t> PhyloBrownianProcessMVN::recursiveComputeCovarianceMatrix(MatrixReal &m, const TopologyNode &node, size_t node_index)
{
    
    // I need to know all my children
    std::set<size_t> children;
    
    if ( node.isRoot() == false )
    {
        // get my scaled branch length
        double v = this->computeBranchTime(node_index, tau->getValue().getBranchLengthForNode(node) );
        
        if ( node.isTip() )
        {
            children.insert( node_index );
            m[node_index][node_index] += v;
        }
        else
        {
            const TopologyNode &left = node.getChild(0);
            size_t left_index = left.getIndex();
            children = recursiveComputeCovarianceMatrix(m, left, left_index );
            
            const TopologyNode &right = node.getChild(1);
            size_t right_index = right.getIndex();
            std::set<size_t> childrenRight = recursiveComputeCovarianceMatrix(m, right, right_index );
            
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
    else // this is the root node
    {
        
        for (size_t i = 0; i < node.getNumberOfChildren(); ++i)
        {
            const TopologyNode &child = node.getChild(i);
            size_t childIndex = child.getIndex();
            std::set<size_t> childrenRight = recursiveComputeCovarianceMatrix(m, child, childIndex );
        }
        
    }
    
    return children;
    
}



void PhyloBrownianProcessMVN::restoreSpecialization( const DagNode* affecter )
{
    
    // reset the flags
    if ( changed_covariance == true )
    {
        changed_covariance = false;
        needs_covariance_recomputation = false;
        
        MatrixReal *tmp = phylogenetic_covariance_matrix;
        phylogenetic_covariance_matrix = stored_phylogenetic_covariance_matrix;
        stored_phylogenetic_covariance_matrix = tmp;

        // The inverse needs to be recalculated.
        inverse_phylogenetic_covariance_matrix = {};
    }
    
}


void PhyloBrownianProcessMVN::setRootState(const TypedDagNode<double> *s)
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


void PhyloBrownianProcessMVN::setRootState(const TypedDagNode<RbVector<double> > *s)
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


std::vector<double> PhyloBrownianProcessMVN::simulateRootCharacters(size_t n)
{
    
    std::vector<double> chars = std::vector<double>(num_sites, 0);
    for (size_t i=0; i<num_sites; ++i)
    {
        chars[i] = computeRootState(i);
    }
    
    return chars;
}


double PhyloBrownianProcessMVN::sumRootLikelihood( void )
{
    if (not inverse_phylogenetic_covariance_matrix)
    {
        inverse_phylogenetic_covariance_matrix = phylogenetic_covariance_matrix->computeInverse();
        inverse_phylogenetic_covariance_matrix->setCholesky( true );
    }
    
    // sum the log-likelihoods for all sites together
    double sum_site_probs = 0.0;
    for (size_t site = 0; site < this->num_sites; ++site)
    {
        std::vector<double> m = std::vector<double>(num_tips, computeRootState(site) );
        
        double sr = this->computeSiteRate(site);
        sum_site_probs += RbStatistics::MultivariateNormal::lnPdfPrecision(m, *inverse_phylogenetic_covariance_matrix, obs[site], sr*sr);
    }
    
    return sum_site_probs;
}


void PhyloBrownianProcessMVN::touchSpecialization( const DagNode* affecter, bool touchAll )
{
    // changing the root state doesn't affect the covariance matrix.
    if ( affecter == homogeneous_root_state or affecter == heterogeneous_root_state )
        return;

    needs_covariance_recomputation = true;
    if ( changed_covariance == false )
    {
        MatrixReal *tmp = phylogenetic_covariance_matrix;
        phylogenetic_covariance_matrix = stored_phylogenetic_covariance_matrix;
        stored_phylogenetic_covariance_matrix = tmp;
        inverse_phylogenetic_covariance_matrix = {};
    }
    changed_covariance = true;
}


/** Swap a parameter of the distribution */
void PhyloBrownianProcessMVN::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
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
        AbstractPhyloBrownianProcess::swapParameterInternal(oldP, newP);
    }
    
}


