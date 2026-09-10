#include <cmath>
#include <cstddef>
#include <set>
#include <vector>

#include "DistributionMultivariateNormal.h"
#include "PhyloMultivariateBrownianProcessREML.h"
#include "RandomNumberFactory.h"
#include "RbException.h"
#include "RbMathLogic.h"
#include "TopologyNode.h"
#include "AbstractPhyloBrownianProcess.h"
#include "Cloneable.h"
#include "ContinuousCharacterData.h"
#include "ContinuousTaxonData.h"
#include "MatrixReal.h"
#include "Tree.h"
#include "TreeChangeEventHandler.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { class RandomNumberGenerator; }


using namespace RevBayesCore;

namespace {

    MatrixReal extractSubmatrix( const MatrixReal &full, const std::vector<size_t> &idx )
    {
        size_t k = idx.size();
        MatrixReal sub(k, k, 0.0);
        for (size_t i = 0; i < k; ++i)
        {
            for (size_t j = 0; j < k; ++j)
            {
                sub[i][j] = full[idx[i]][idx[j]];
            }
        }
        return sub;
    }

    std::vector<double> extractSubvector( const std::vector<double> &full, const std::vector<size_t> &idx )
    {
        size_t k = idx.size();
        std::vector<double> sub(k, 0.0);
        for (size_t i = 0; i < k; ++i)
        {
            sub[i] = full[idx[i]];
        }
        return sub;
    }

}

PhyloMultivariateBrownianProcessREML::PhyloMultivariateBrownianProcessREML(const TypedDagNode<Tree> *t, const TypedDagNode<MatrixReal> *c, size_t ns) :
    AbstractPhyloBrownianProcess( t, ns ),
    partial_likelihoods( std::vector<std::vector<double> >(2, std::vector<double>(this->num_nodes, 0) ) ),
    contrasts( std::vector<std::vector<std::vector<double> > >(2, std::vector<std::vector<double> >(this->num_nodes, std::vector<double>(this->num_sites, 0.0) ) ) ),
    contrast_uncertainty( std::vector<std::vector<MatrixReal> >(2, std::vector<MatrixReal>(this->num_nodes, MatrixReal(this->num_sites, this->num_sites, 0.0) ) ) ),
    active_likelihood( std::vector<size_t>(this->num_nodes, 0) ),
    independent_contrasts( std::vector<std::vector<double> >(this->num_nodes, std::vector<double>(this->num_sites, 0.0) ) ),
    independent_contrasts_sds( std::vector<double>(this->num_nodes, 0.0) ),
    observed_dims( std::vector<std::vector<bool> >(this->num_nodes, std::vector<bool>(this->num_sites, false) ) ),
    is_pristine( std::vector<bool>(this->num_nodes, true) ),
    pristine_time( std::vector<std::vector<double> >(2, std::vector<double>(this->num_nodes, 0.0) ) ),
    changed_nodes( std::vector<bool>(this->num_nodes, false) ),
    dirty_nodes( std::vector<bool>(this->num_nodes, true) ),
    rate_matrix( c ),
    active_matrix(0),
    precision_matrices( std::vector<MatrixReal>( 2, MatrixReal(num_sites) ) )
{
    
    // add the parameters to our set
    this->addParameter( rate_matrix );
    
    // make sure the rate matrix is inverted using Cholesky decomposition
    rate_matrix->getValue().setCholesky(true);
    
    // compute the inverse variance-covariance matrix (the precision matrix)
    precision_matrices[0] = rate_matrix->getValue().computeInverse();
    precision_matrices[0].setCholesky(true);
    precision_matrices[1] = rate_matrix->getValue().computeInverse();
    precision_matrices[1].setCholesky(true);
    
    // We don't want tau to die before we die, or it can't remove us as listener
    tau->getValue().getTreeChangeEventHandler().addListener( this );

    // now we need to reset the values
    this->redrawValue();
    
    // we need to reset the contrasts
    resetValue();

}


/**
 * Destructor. Because we added ourselves as a reference to tau when we added a listener to its
 * TreeChangeEventHandler, we need to remove ourselves as a reference and possibly delete tau
 * when we die. All other parameters are handled by others.
 */
PhyloMultivariateBrownianProcessREML::~PhyloMultivariateBrownianProcessREML( void )
{
    // We don't delete the params, because they might be used somewhere else too. The model needs to do that!
    
    // remove myself from the tree listeners
    if ( tau != NULL )
    {
        tau->getValue().getTreeChangeEventHandler().removeListener( this );
    }
    
}



PhyloMultivariateBrownianProcessREML* PhyloMultivariateBrownianProcessREML::clone( void ) const
{
    
    return new PhyloMultivariateBrownianProcessREML( *this );
}


double PhyloMultivariateBrownianProcessREML::computeLnProbability( void )
{
    
    // we need to check here if we still are listining to this tree for change events
    // the tree could have been replaced without telling us
    if ( tau->getValue().getTreeChangeEventHandler().isListening( this ) == false )
    {
        tau->getValue().getTreeChangeEventHandler().addListener( this );
        dirty_nodes = std::vector<bool>(tau->getValue().getNumberOfNodes(), true);
    }
    
    // compute the ln probability by recursively calling the probability calculation for each node
    const TopologyNode &root = this->tau->getValue().getRoot();
    
    // we start with the root and then traverse down the tree
    size_t rootIndex = root.getIndex();
    
    // only necessary if the root is actually dirty
    if ( this->dirty_nodes[rootIndex] )
        recursiveComputeLnProbability( root, rootIndex );

    // return the likelihood at the root
    this->ln_prob = this->partial_likelihoods[this->active_likelihood[rootIndex]][rootIndex];
        
    
    return this->ln_prob;
    
}



void PhyloMultivariateBrownianProcessREML::fireTreeChangeEvent( const TopologyNode &n, const unsigned& m )
{
    
    // call a recursive flagging of all node above (closer to the root) and including this node
    recursivelyFlagNodeDirty( n );
    
}



std::vector<std::vector<double> > PhyloMultivariateBrownianProcessREML::getContrasts()
{

    // make sure the necessary quantities are clean
    computeLnProbability();
    
    // clean out the contrasts vector
    independent_contrasts.clear();
    independent_contrasts = std::vector<std::vector<double> >(this->num_nodes, std::vector<double>(this->num_sites, 0.0) );

    independent_contrasts_sds.clear();
    independent_contrasts_sds = std::vector<double>(this->num_nodes, 0.0);

    // compute the constrasts by recursively calling the contrast calculation for each node
    const TopologyNode &root = this->tau->getValue().getRoot();

    // we start with the root and then traverse down the tree
    size_t rootIndex = root.getIndex();

    // compute the contrast at the root, and recursively compute the remaining contrasts up the tree
    recursiveComputeContrasts( root, rootIndex );
    
    // remove the first n standardized contrasts, where n is the number of tips
    size_t num_tips = this->tau->getValue().getNumberOfTips();
    
    independent_contrasts.erase(independent_contrasts.begin(), independent_contrasts.begin() + num_tips);
    independent_contrasts_sds.erase(independent_contrasts_sds.begin(), independent_contrasts_sds.begin() + num_tips);

    return independent_contrasts;

}


void PhyloMultivariateBrownianProcessREML::executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<RbVector<double> > &rv) const
{

    if ( n == "contrasts" )
    {
        std::vector<std::vector<double> > c = const_cast<PhyloMultivariateBrownianProcessREML*>(this)->getContrasts();

        rv.clear();
        rv.resize( c.size() );
        for (size_t i = 0; i < c.size(); ++i)
        {
            rv[i] = RbVector<double>( c[i] );
        }
    }
    else
    {
        throw RbException() << "The phylogenetic multivariate Brownian motion process does not have a member method called '" << n << "'.";
    }

}


void PhyloMultivariateBrownianProcessREML::keepSpecialization( const DagNode* affecter )
{
    
    // reset all flags
    for (std::vector<bool>::iterator it = this->dirty_nodes.begin(); it != this->dirty_nodes.end(); ++it)
    {
        (*it) = false;
    }
    
    for (std::vector<bool>::iterator it = this->changed_nodes.begin(); it != this->changed_nodes.end(); ++it)
    {
        (*it) = false;
    }
    
}


void PhyloMultivariateBrownianProcessREML::recursiveComputeLnProbability( const TopologyNode &node, size_t node_index )
{

    // check for recomputation
    if ( node.isTip() == false && dirty_nodes[node_index] )
    {
        // mark as computed
        dirty_nodes[node_index] = false;

        double               &p_node   = this->partial_likelihoods[this->active_likelihood[node_index]][node_index];
        std::vector<double>  &mu_node  = this->contrasts[this->active_likelihood[node_index]][node_index];
        MatrixReal           &Omega_node = this->contrast_uncertainty[this->active_likelihood[node_index]][node_index];
        std::vector<bool>    &obs_node = this->observed_dims[node_index];
        double               &t_node   = this->pristine_time[this->active_likelihood[node_index]][node_index];

        const MatrixReal &Sigma = this->rate_matrix->getValue();

        // get the number of children
        size_t num_children = node.getNumberOfChildren();

        for (size_t j = 1; j < num_children; ++j)
        {

            size_t left_index = node_index;
            const TopologyNode *left = &node;
            if ( j == 1 )
            {
                left = &node.getChild(0);
                left_index = left->getIndex();
                recursiveComputeLnProbability( *left, left_index );
            }

            const TopologyNode &right = node.getChild(j);
            size_t right_index = right.getIndex();
            recursiveComputeLnProbability( right, right_index );

            const double &p_left  = this->partial_likelihoods[this->active_likelihood[left_index]][left_index];
            const double &p_right = this->partial_likelihoods[this->active_likelihood[right_index]][right_index];

            // get the per node and site contrasts
            const std::vector<double> &mu_left  = this->contrasts[this->active_likelihood[left_index]][left_index];
            const std::vector<double> &mu_right = this->contrasts[this->active_likelihood[right_index]][right_index];

            // get the propagated uncertainties
            const MatrixReal &Omega_left  = this->contrast_uncertainty[this->active_likelihood[left_index]][left_index];
            const MatrixReal &Omega_right = this->contrast_uncertainty[this->active_likelihood[right_index]][right_index];

            const std::vector<bool> &obs_left  = this->observed_dims[left_index];
            const std::vector<bool> &obs_right = this->observed_dims[right_index];

            bool pristine_left    = this->is_pristine[left_index];
            bool pristine_right   = this->is_pristine[right_index];
            double t_left_start   = this->pristine_time[this->active_likelihood[left_index]][left_index];
            double t_right_start  = this->pristine_time[this->active_likelihood[right_index]][right_index];

            // get the scaled branch lengths
            double v_left  = 0;
            if ( j == 1 )
            {
                v_left = this->computeBranchTime(left_index, left->getBranchLength());
            }
            double v_right = this->computeBranchTime(right_index, right.getBranchLength());

            size_t nL = 0;
            for (size_t i = 0; i < this->num_sites; ++i)
            {
                if ( obs_left[i] )
                {
                    ++nL;
                }
            }
            size_t nR = 0;
            for (size_t i = 0; i < this->num_sites; ++i)
            {
                if ( obs_right[i] )
                {
                    ++nR;
                }
            }

            if ( nR == 0 )
            {
                p_node = p_left + p_right;

                if ( nL == 0 )
                {
                    Omega_node = MatrixReal(this->num_sites, this->num_sites, 0.0);
                    for (size_t i = 0; i < this->num_sites; ++i)
                    {
                        obs_node[i] = false;
                    }
                    is_pristine[node_index] = true;
                    t_node = 0.0;
                }
                else if ( pristine_left )
                {
                    t_node = v_left + t_left_start;
                    for (size_t i = 0; i < this->num_sites; ++i)
                    {
                        obs_node[i] = obs_left[i];
                        if ( obs_left[i] )
                        {
                            mu_node[i] = mu_left[i];
                        }
                    }
                    Omega_node = Sigma * t_node;
                    is_pristine[node_index] = true;
                }
                else
                {
                    Omega_node = Omega_left;
                    for (size_t i = 0; i < this->num_sites; ++i)
                    {
                        obs_node[i] = obs_left[i];
                        if ( !obs_left[i] )
                        {
                            continue;
                        }
                        mu_node[i] = mu_left[i];
                        for (size_t jx = 0; jx < this->num_sites; ++jx)
                        {
                            if ( obs_left[jx] )
                            {
                                Omega_node[i][jx] = Omega_left[i][jx] + v_left * Sigma[i][jx];
                            }
                        }
                    }
                    is_pristine[node_index] = false;
                }
                continue;
            }
            else if ( nL == 0 )
            {
                p_node = p_left + p_right;

                if ( pristine_right )
                {
                    t_node = v_right + t_right_start;
                    for (size_t i = 0; i < this->num_sites; ++i)
                    {
                        obs_node[i] = obs_right[i];
                        if ( obs_right[i] )
                        {
                            mu_node[i] = mu_right[i];
                        }
                    }
                    Omega_node = Sigma * t_node;
                    is_pristine[node_index] = true;
                }
                else
                {
                    Omega_node = Omega_right;
                    for (size_t i = 0; i < this->num_sites; ++i)
                    {
                        obs_node[i] = obs_right[i];
                        if ( !obs_right[i] )
                        {
                            continue;
                        }
                        mu_node[i] = mu_right[i];
                        for (size_t jx = 0; jx < this->num_sites; ++jx)
                        {
                            if ( obs_right[jx] )
                            {
                                Omega_node[i][jx] = Omega_right[i][jx] + v_right * Sigma[i][jx];
                            }
                        }
                    }
                    is_pristine[node_index] = false;
                }
                continue;
            }
            else if ( nL == this->num_sites && nR == this->num_sites && pristine_left && pristine_right )
            {
                double t_left_prop  = v_left  + t_left_start;
                double t_right_prop = v_right + t_right_start;
                t_node = (t_left_prop * t_right_prop) / (t_left_prop + t_right_prop);
                double branch_length = t_left_prop + t_right_prop;

                std::vector<double> these_contrasts(this->num_sites);
                std::vector<double> means(this->num_sites, 0.0);
                for (size_t i = 0; i < this->num_sites; ++i)
                {
                    these_contrasts[i] = mu_left[i] - mu_right[i];
                    mu_node[i] = (mu_left[i] * t_right_prop + mu_right[i] * t_left_prop) / (t_left_prop + t_right_prop);
                    obs_node[i] = true;
                }

                double lnl_contrast = RbStatistics::MultivariateNormal::lnPdfPrecision(means, precision_matrices[active_matrix], these_contrasts, branch_length);
                p_node = lnl_contrast + p_left + p_right;

                Omega_node = Sigma * t_node;
                is_pristine[node_index] = true;
                continue;
            }

            is_pristine[node_index] = false;

            std::vector<size_t> idx_left, idx_right, idx_common, idx_node;
            std::vector<size_t> node_pos(this->num_sites, 0);
            for (size_t i = 0; i < this->num_sites; ++i)
            {
                bool L = obs_left[i];
                bool R = obs_right[i];
                if ( L )
                {
                    idx_left.push_back(i);
                }
                if ( R )
                {
                    idx_right.push_back(i);
                }
                if ( L || R )
                {
                    node_pos[i] = idx_node.size();
                    idx_node.push_back(i);
                }
                if ( L && R )
                {
                    idx_common.push_back(i);
                }
            }

            size_t nN = idx_node.size();
            size_t nC = idx_common.size();

            MatrixReal Omega_left_prop( nL, nL, 0.0 );
            if ( nL > 0 )
            {
                Omega_left_prop = extractSubmatrix(Omega_left, idx_left) + extractSubmatrix(Sigma, idx_left) * v_left;
            }
            MatrixReal Omega_right_prop( nR, nR, 0.0 );
            if ( nR > 0 )
            {
                Omega_right_prop = extractSubmatrix(Omega_right, idx_right) + extractSubmatrix(Sigma, idx_right) * v_right;
            }

            std::vector<size_t> left_pos(this->num_sites, 0);
            for (size_t i = 0; i < nL; ++i)
            {
                left_pos[idx_left[i]] = i;
            }
            std::vector<size_t> right_pos(this->num_sites, 0);
            for (size_t i = 0; i < nR; ++i)
            {
                right_pos[idx_right[i]] = i;
            }
            double lnl_contrast = 0.0;
            if ( nC > 0 )
            {
                std::vector<double> contrast_vec(nC, 0.0);
                std::vector<double> zeros(nC, 0.0);
                MatrixReal Omega_common(nC, nC, 0.0);
                for (size_t a = 0; a < nC; ++a)
                {
                    size_t gi = idx_common[a];
                    contrast_vec[a] = mu_left[gi] - mu_right[gi];
                    size_t li = left_pos[gi];
                    size_t ri = right_pos[gi];
                    for (size_t b = 0; b < nC; ++b)
                    {
                        size_t gj = idx_common[b];
                        Omega_common[a][b] = Omega_left_prop[li][left_pos[gj]] + Omega_right_prop[ri][right_pos[gj]];
                    }
                }

                lnl_contrast = RbStatistics::MultivariateNormal::lnPdfCovariance(zeros, Omega_common, contrast_vec, 1.0);
            }
            p_node = lnl_contrast + p_left + p_right;

            Omega_node = MatrixReal(this->num_sites, this->num_sites, 0.0);
            for (size_t i = 0; i < this->num_sites; ++i)
            {
                obs_node[i] = false;
            }

            if ( nN > 0 )
            {
                MatrixReal Precision_node(nN, nN, 0.0);
                std::vector<double> eta_node(nN, 0.0);

                if ( nL > 0 )
                {
                    Omega_left_prop.setCholesky(true);
                    MatrixReal Precision_left = Omega_left_prop.computeInverse();
                    std::vector<double> mu_left_sub = extractSubvector(mu_left, idx_left);

                    for (size_t i = 0; i < nL; ++i)
                    {
                        size_t pi = node_pos[idx_left[i]];
                        double eta_i = 0.0;
                        for (size_t j = 0; j < nL; ++j)
                        {
                            double pij = Precision_left[i][j];
                            eta_i += pij * mu_left_sub[j];
                            Precision_node[pi][node_pos[idx_left[j]]] += pij;
                        }
                        eta_node[pi] += eta_i;
                    }
                }

                if ( nR > 0 )
                {
                    Omega_right_prop.setCholesky(true);
                    MatrixReal Precision_right = Omega_right_prop.computeInverse();
                    std::vector<double> mu_right_sub = extractSubvector(mu_right, idx_right);

                    for (size_t i = 0; i < nR; ++i)
                    {
                        size_t pi = node_pos[idx_right[i]];
                        double eta_i = 0.0;
                        for (size_t j = 0; j < nR; ++j)
                        {
                            double pij = Precision_right[i][j];
                            eta_i += pij * mu_right_sub[j];
                            Precision_node[pi][node_pos[idx_right[j]]] += pij;
                        }
                        eta_node[pi] += eta_i;
                    }
                }

                Precision_node.setCholesky(true);
                MatrixReal Omega_node_sub = Precision_node.computeInverse();

                for (size_t i = 0; i < nN; ++i)
                {
                    double mu_i = 0.0;
                    for (size_t j = 0; j < nN; ++j)
                    {
                        mu_i += Omega_node_sub[i][j] * eta_node[j];
                    }
                    mu_node[idx_node[i]] = mu_i;

                    for (size_t j = 0; j < nN; ++j)
                    {
                        Omega_node[idx_node[i]][idx_node[j]] = Omega_node_sub[i][j];
                    }
                    obs_node[idx_node[i]] = true;
                }
            }

        } // end for-loop over all children

    } // end if we need to compute something for this node.

}


void PhyloMultivariateBrownianProcessREML::recursiveComputeContrasts( const TopologyNode &node, size_t node_index )
{
    
    // nothing to do for tips
    if ( node.isTip() == false )
    {

        std::vector<double> &mu_node = this->independent_contrasts[node_index];
        double              &sd_node = this->independent_contrasts_sds[node_index];

        const MatrixReal &Sigma = this->rate_matrix->getValue();

        // get the number of children
        size_t num_children = node.getNumberOfChildren();

        for (size_t j = 1; j < num_children; ++j)
        {

            size_t left_index = node_index;
            const TopologyNode *left = &node;
            if ( j == 1 )
            {
                left = &node.getChild(0);
                left_index = left->getIndex();
                recursiveComputeContrasts( *left, left_index );
            }

            const TopologyNode &right = node.getChild(j);
            size_t right_index = right.getIndex();
            recursiveComputeContrasts( right, right_index );

            // get the per node and site contrasts
            const std::vector<double> &mu_left  = this->contrasts[this->active_likelihood[left_index]][left_index];
            const std::vector<double> &mu_right = this->contrasts[this->active_likelihood[right_index]][right_index];

            const MatrixReal &Omega_left  = this->contrast_uncertainty[this->active_likelihood[left_index]][left_index];
            const MatrixReal &Omega_right = this->contrast_uncertainty[this->active_likelihood[right_index]][right_index];
            const std::vector<bool> &obs_left  = this->observed_dims[left_index];
            const std::vector<bool> &obs_right = this->observed_dims[right_index];

            double delta_left   = 0.0;
            double delta_right  = 0.0;
            size_t n_obs_left   = 0;
            size_t n_obs_right  = 0;
            for (size_t i = 0; i < this->num_sites; ++i)
            {
                if ( obs_left[i] )
                {
                    delta_left += Omega_left[i][i] / Sigma[i][i];
                    ++n_obs_left;
                }
                if ( obs_right[i] )
                {
                    delta_right += Omega_right[i][i] / Sigma[i][i];
                    ++n_obs_right;
                }
            }
            if ( n_obs_left > 0 )
            {
                delta_left /= n_obs_left;
            }
            if ( n_obs_right > 0 )
            {
                delta_right /= n_obs_right;
            }

            // get the scaled branch lengths
            double v_left  = 0;
            if ( j == 1 )
            {
                v_left = this->computeBranchTime(left_index, left->getBranchLength());
            }
            double v_right = this->computeBranchTime(right_index, right.getBranchLength());
            
            // add the propagated uncertainty to the branch lengths
            double t_left        = v_left  + delta_left;
            double t_right       = v_right + delta_right;
            double branch_length = t_left + t_right;
            
            sd_node = pow(branch_length, 0.5);
            for (size_t i = 0; i < this->num_sites; ++i)
            {
//                mu_node[i] = (mu_left[i] - mu_right[i]) / pow(branch_length, 0.5);
//                mu_node[i] = (mu_left[i] - mu_right[i]) / branch_length;
                mu_node[i] = mu_left[i] - mu_right[i];
            }
            
        } // end for-loop over all children
        
    } // end if we need to compute something for this node.
    
}


void PhyloMultivariateBrownianProcessREML::recursivelyFlagNodeDirty( const TopologyNode &n )
{
    
    // we need to flag this node and all ancestral nodes for recomputation
    size_t index = n.getIndex();
    
    // if this node is already dirty, the also all the ancestral nodes must have been flagged as dirty
    if ( !dirty_nodes[index] )
    {
        // the root doesn't have an ancestor
        if ( !n.isRoot() )
        {
            recursivelyFlagNodeDirty( n.getParent() );
        }
        
        // set the flag
        dirty_nodes[index] = true;
        
        // if we previously haven't touched this node, then we need to change the active likelihood pointer
        if ( changed_nodes[index] == false )
        {
            active_likelihood[index] = (active_likelihood[index] == 0 ? 1 : 0);
            changed_nodes[index] = true;
        }
        
    }
    
}


void PhyloMultivariateBrownianProcessREML::resetValue( void )
{
    
    // check if the vectors need to be resized
    partial_likelihoods = std::vector<std::vector<double> >(2, std::vector<double>(this->num_nodes, 0) );
    contrasts = std::vector<std::vector<std::vector<double> > >(2, std::vector<std::vector<double> >(this->num_nodes, std::vector<double>(this->num_sites, 0) ) );
    contrast_uncertainty = std::vector<std::vector<MatrixReal> >(2, std::vector<MatrixReal>(this->num_nodes, MatrixReal(this->num_sites, this->num_sites, 0.0) ) );
    observed_dims = std::vector<std::vector<bool> >(this->num_nodes, std::vector<bool>(this->num_sites, false) );
    is_pristine = std::vector<bool>(this->num_nodes, true);
    pristine_time = std::vector<std::vector<double> >(2, std::vector<double>(this->num_nodes, 0.0) );

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

    std::vector<TopologyNode*> nodes = this->tau->getValue().getNodes();
    for (size_t site = 0; site < this->num_sites; ++site)
    {

        for (std::vector<TopologyNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it)
        {
            if ( (*it)->isTip() )
            {
                ContinuousTaxonData& taxon = this->value->getTaxonData( (*it)->getName() );
                double &c = taxon.getCharacter(site_indices[site]);

                size_t node_index = (*it)->getIndex();
                bool is_observed = RbMath::isFinite(c);
                observed_dims[node_index][site] = is_observed;

                double stored_value = is_observed ? c : 0.0;
                contrasts[0][node_index][site] = stored_value;
                contrasts[1][node_index][site] = stored_value;
            }
        }
    }

    // finally we set all the flags for recomputation
    for (std::vector<bool>::iterator it = dirty_nodes.begin(); it != dirty_nodes.end(); ++it)
    {
        (*it) = true;
    }
    
    // flip the active likelihood pointers
    for (size_t index = 0; index < changed_nodes.size(); ++index)
    {
        active_likelihood[index] = 0;
        changed_nodes[index] = true;
    }
    
}


void PhyloMultivariateBrownianProcessREML::restoreSpecialization( const DagNode* affecter )
{
    
    // reset the precision matrix if necessary
    if ( affecter == rate_matrix )
    {
        active_matrix = (active_matrix == 0 ? 1 : 0);
    }
    
    // reset the flags
    for (std::vector<bool>::iterator it = dirty_nodes.begin(); it != dirty_nodes.end(); ++it)
    {
        (*it) = false;
    }
    
    // restore the active likelihoods vector
    for (size_t index = 0; index < changed_nodes.size(); ++index)
    {
        // we have to restore, that means if we have changed the active likelihood vector
        // then we need to revert this change
        if ( changed_nodes[index] == true )
        {
            active_likelihood[index] = (active_likelihood[index] == 0 ? 1 : 0);
        }
        
        // set all flags to false
        changed_nodes[index] = false;
    }
    
}


void PhyloMultivariateBrownianProcessREML::simulateRecursively( const TopologyNode &node, std::vector< ContinuousTaxonData > &taxa)
{
    
    // get the children of the node
    const std::vector<TopologyNode*>& children = node.getChildren();
    
    // get the sequence of this node
    size_t node_index = node.getIndex();
    const ContinuousTaxonData &parent = taxa[ node_index ];
    
    std::vector<double> parent_state(num_sites, 0.0);
    for (size_t i = 0; i < num_sites; ++i)
    {
        parent_state[i] = parent.getCharacter(i);
    }
    
    // simulate the sequence for each child
    RandomNumberGenerator* rng = GLOBAL_RNG;
    for (std::vector< TopologyNode* >::const_iterator it = children.begin(); it != children.end(); ++it)
    {
        const TopologyNode &child = *(*it);
        
        // get the branch length for this child
        size_t child_index   = child.getIndex();
        double branch_length = this->computeBranchTime(child_index, child.getBranchLength());
        
        ContinuousTaxonData &taxon = taxa[ child.getIndex() ];
//        std::vector<double> c = RbStatistics::MultivariateNormal::rvPrecision(parent_state, precision_matrices[active_matrix], *rng, branch_length);
        std::vector<double> c = RbStatistics::MultivariateNormal::rvCovariance(parent_state, rate_matrix->getValue(), *rng, branch_length);
        
        for ( size_t i = 0; i < num_sites; ++i )
        {
            // add the character to the sequence
            taxon.addCharacter( c[i] );
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


std::vector<double> PhyloMultivariateBrownianProcessREML::simulateRootCharacters(size_t n)
{
    
    std::vector<double> chars = std::vector<double>(num_sites, 0);
    for (size_t i=0; i<num_sites; ++i)
    {
        chars[i] = 0.0;
    }
    
    return chars;
}


void PhyloMultivariateBrownianProcessREML::touchSpecialization( const DagNode* affecter, bool touchAll )
{
 
    // if the topology wasn't the culprit for the touch, then we just flag everything as dirty
    if ( affecter == this->heterogeneous_clock_rates )
    {
        
        const std::set<size_t> &indices = this->heterogeneous_clock_rates->getTouchedElementIndices();
        
        // maybe all of them have been touched or the flags haven't been set properly
        if ( indices.size() == 0 )
        {
            // just flag everything for recomputation
            touchAll = true;
        }
        else
        {
            const std::vector<TopologyNode *> &nodes = this->tau->getValue().getNodes();
            // flag recomputation only for the nodes
            for (std::set<size_t>::iterator it = indices.begin(); it != indices.end(); ++it)
            {
                this->recursivelyFlagNodeDirty( *nodes[*it] );
            }
        }
    }
    else if ( affecter == rate_matrix )
    {
        // compute the inverse variance-covariance matrix (the precision matrix)
        active_matrix = (active_matrix == 0 ? 1 : 0);
        precision_matrices[active_matrix] = rate_matrix->getValue().computeInverse();
        precision_matrices[active_matrix].setCholesky(true);
        
        // we need to recompute the likelihood
        touchAll = true;
    }
    else if ( affecter != this->tau ) // if the topology wasn't the culprit for the touch, then we just flag everything as dirty
    {
        touchAll = true;
    }
    
    if ( touchAll )
    {
        for (std::vector<bool>::iterator it = dirty_nodes.begin(); it != dirty_nodes.end(); ++it)
        {
            (*it) = true;
        }
        
        // flip the active likelihood pointers
        for (size_t index = 0; index < changed_nodes.size(); ++index)
        {
            if ( changed_nodes[index] == false )
            {
                active_likelihood[index] = (active_likelihood[index] == 0 ? 1 : 0);
                changed_nodes[index] = true;
            }
        }
        
    }
    
}


/** Swap a parameter of the distribution */
void PhyloMultivariateBrownianProcessREML::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == this->tau)
    {
        this->tau->getValue().getTreeChangeEventHandler().removeListener( this );
        AbstractPhyloBrownianProcess::swapParameterInternal(oldP, newP);
        this->tau->getValue().getTreeChangeEventHandler().addListener( this );
    }
    if (oldP == this->rate_matrix)
    {
        rate_matrix = static_cast<const TypedDagNode< MatrixReal >* >( newP );
        active_matrix = 0;
        precision_matrices[0] = rate_matrix->getValue().computeInverse();
        precision_matrices[0].setCholesky(true);
        precision_matrices[1] = rate_matrix->getValue().computeInverse();
        precision_matrices[1].setCholesky(true);
    }
    else
    {
        AbstractPhyloBrownianProcess::swapParameterInternal(oldP, newP);
    }
    
}

