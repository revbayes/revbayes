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
    obs( std::vector<std::vector<double> >(this->num_sites, std::vector<double>(num_tips, 0.0) ) )
{
    homogeneous_root_state      = new ConstantNode<double>("", new double(0.0) );
    heterogeneous_root_state    = NULL;

    addParameter( homogeneous_root_state );
    
    // now we need to reset the value
    this->redrawValue();
}


PhyloBrownianProcessMVN::PhyloBrownianProcessMVN(const PhyloBrownianProcessMVN &p) : AbstractPhyloBrownianProcess( p ),
    homogeneous_root_state( p.homogeneous_root_state ),
    heterogeneous_root_state( p.heterogeneous_root_state ),
    num_tips( p.num_tips ),
    obs( p.obs ),
    covariance_cache( p.covariance_cache )
{
    
}


/**
 * Destructor. Because we added ourselves as a reference to tau when we added a listener to its
 * TreeChangeEventHandler, we need to remove ourselves as a reference and possibly delete tau
 * when we die. All other parameters are handled by others.
 */
PhyloBrownianProcessMVN::~PhyloBrownianProcessMVN( void )
{
    
}



PhyloBrownianProcessMVN& PhyloBrownianProcessMVN::operator=(const PhyloBrownianProcessMVN &p)

{
    
    if ( this != &p )
    {
        AbstractPhyloBrownianProcess::operator=( p );

        homogeneous_root_state                      = p.homogeneous_root_state;
        heterogeneous_root_state                    = p.heterogeneous_root_state;
        num_tips                                    = p.num_tips;
        obs                                         = p.obs;
        covariance_cache                            = p.covariance_cache;
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
    
    if ( covariance_cache.is_valid() == false )
    {
        // perhaps there is a more efficient way to reset the matrix to 0.
        CovarianceCache &cache = covariance_cache.init_for_writing();
        cache.covariance = MatrixReal(num_tips, num_tips);
        cache.covariance.setCholesky( true );
        recursiveComputeCovarianceMatrix(cache.covariance, root, rootIndex);
        cache.inverse_covariance = {};
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
    if ( covariance_cache.has_snapshot() )
        covariance_cache.keep();
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
    
    
    covariance_cache.clear();
    
}


std::set<size_t> PhyloBrownianProcessMVN::recursiveComputeCovarianceMatrix(MatrixReal &m, const TopologyNode &node, size_t node_index)
{
    
    // I need to know all my children
    std::set<size_t> children;
    
    if ( node.isRoot() == false )
    {
        // get my scaled branch length
        double v = this->computeBranchTime(node_index, node.getBranchLength() );
        
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
    if ( covariance_cache.has_snapshot() )
        covariance_cache.restore();
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
    CovarianceCache &cache = covariance_cache.get_mutable();
    if (not cache.inverse_covariance)
    {
        cache.inverse_covariance = cache.covariance.computeInverse();
        cache.inverse_covariance->setCholesky( true );
    }
    
    // sum the log-likelihoods for all sites together
    double sum_site_probs = 0.0;
    for (size_t site = 0; site < this->num_sites; ++site)
    {
        std::vector<double> m = std::vector<double>(num_tips, computeRootState(site) );
        
        double sr = this->computeSiteRate(site);
        sum_site_probs += RbStatistics::MultivariateNormal::lnPdfPrecision(m, *cache.inverse_covariance, obs[site], sr*sr);
    }
    
    return sum_site_probs;
}

void PhyloBrownianProcessMVN::snapshotSpecialization( void )
{
    covariance_cache.snapshot();
}


/*
 * Mark the Brownian MVN covariance cache invalid after non-root dependencies change.
 * Root-state proposals alter only the mean vector built while summing the likelihood.
 */
void PhyloBrownianProcessMVN::invalidateSpecialization( const DagNode* affecter, bool touchAll )
{
    // changing the root state doesn't affect the covariance matrix.
    if ( affecter == homogeneous_root_state or affecter == heterogeneous_root_state )
        return;

    covariance_cache.invalidate();
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

