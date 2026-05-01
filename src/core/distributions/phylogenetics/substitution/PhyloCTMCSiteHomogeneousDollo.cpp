#include "PhyloCTMCSiteHomogeneousDollo.h"

#include <climits>
#include <cmath>
#include <algorithm>
#include <bitset>
#include <string>
#include <iosfwd>

#include "RbMathFunctions.h"
#include "AbstractHomologousDiscreteCharacterData.h"
#include "DiscreteCharacterState.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RateGenerator.h"
#include "RateMatrix.h"
#include "RbBitSet.h"
#include "RbConstants.h"
#include "RbMathLogic.h"
#include "RbSettings.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "StochasticNode.h"
#include "Taxon.h"
#include "TransitionProbabilityMatrix.h"
#include "Tree.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

RevBayesCore::PhyloCTMCSiteHomogeneousDollo::PhyloCTMCSiteHomogeneousDollo(const TypedDagNode<Tree> *t, size_t nc, bool c, size_t nSites, bool amb, DolloAscertainmentBias::Coding ty, bool norm) :
    PhyloCTMCSiteHomogeneousConditional<StandardState>(  t, nc + 1, c, nSites, amb, AscertainmentBias::Coding(ty)), dim(nc), integrationFactors(0), normalize(norm)
{
    massNodeOffset = this->num_site_mixtures*numCorrectionMasks;
    activeMassOffset = this->num_nodes*massNodeOffset;
    perMaskMixtureCorrections = std::vector<double>(2*activeMassOffset, 0.0);

    death_rate = NULL;
    this->addParameter(death_rate);

    survival = std::vector<double>(this->num_site_mixtures,0.0);
    integrationFactors = std::vector<double>(this->num_site_mixtures,0.0);
    maskNodeObservationCounts = std::vector<std::vector<size_t> >(numCorrectionMasks, std::vector<size_t>(num_nodes, 0) );

    // set the offsets for easier iteration through the likelihood vector
    // we use an extra site for the integrated likelihood
    siteOffset                  =  dim + 2;
    mixtureOffset               =  pattern_block_size*siteOffset;
    nodeOffset                  =  num_site_mixtures*mixtureOffset;
    activeLikelihoodOffset      =  num_nodes*nodeOffset;

}

RevBayesCore::PhyloCTMCSiteHomogeneousDollo::PhyloCTMCSiteHomogeneousDollo(const PhyloCTMCSiteHomogeneousDollo& n) :
        PhyloCTMCSiteHomogeneousConditional<StandardState>(n),
        dim(n.dim),
        integrationFactors(n.integrationFactors),
        maskNodeObservationCounts(n.maskNodeObservationCounts),
        survival(n.survival),
        activeMassOffset(n.activeMassOffset),
        massNodeOffset(n.massNodeOffset),
        normalize(n.normalize),
        death_rate(n.death_rate)
{

}

RevBayesCore::PhyloCTMCSiteHomogeneousDollo* RevBayesCore::PhyloCTMCSiteHomogeneousDollo::clone( void ) const {

    return new PhyloCTMCSiteHomogeneousDollo( *this );
}

void RevBayesCore::PhyloCTMCSiteHomogeneousDollo::resizeLikelihoodVectors( void )
{
    RevBayesCore::PhyloCTMCSiteHomogeneousConditional<StandardState>::resizeLikelihoodVectors();

    massNodeOffset = num_site_mixtures*numCorrectionMasks;
    activeMassOffset = num_nodes*massNodeOffset;
    perMaskMixtureCorrections = std::vector<double>(2*activeMassOffset, 0.0);

    maskNodeObservationCounts = std::vector<std::vector<size_t> >(numCorrectionMasks, std::vector<size_t>(num_nodes, 0) );

    integrationFactors = std::vector<double>(num_site_mixtures, 0.0);
    survival = std::vector<double>(this->num_site_mixtures, 0.0);

    // set the offsets for easier iteration through the likelihood vector
    siteOffset                  =  dim + 2;
    mixtureOffset               =  pattern_block_size*siteOffset;
    nodeOffset                  =  num_site_mixtures*mixtureOffset;
    activeLikelihoodOffset      =  num_nodes*nodeOffset;

    // only do this if we are in MCMC mode. This will safe memory
    if ( in_mcmc_mode == true )
    {
        // we resize the partial likelihood vectors to the new dimensions
        delete [] partialLikelihoods;

        partialLikelihoods = new double[2*activeLikelihoodOffset];

        // reinitialize likelihood vectors
        for (size_t i = 0; i < 2*activeLikelihoodOffset; i++)
        {
            partialLikelihoods[i] = 0.0;
        }
    }
}

void RevBayesCore::PhyloCTMCSiteHomogeneousDollo::updateTransitionProbabilities(size_t node_idx)
{
    // first, get the clock rate for the branch
    double rate = 1.0;
    if ( this->branch_heterogeneous_clock_rates == true )
    {
        rate = this->heterogeneous_clock_rates->getValue()[node_idx];
    }
    else if (this->homogeneous_clock_rate != NULL)
    {
        rate = this->homogeneous_clock_rate->getValue();
    }

    // get the death rate
    double death = 1.0;
    if (death_rate != NULL)
    {
        death = death_rate->getValue();
    }

    // and compute the per site rate transition probability matrix
    const TopologyNode* node = tau->getValue().getNodes()[node_idx];

    double endAge = node->getAge();

    // if the tree is not a time tree, then the age will be not a number
    if ( RbMath::isFinite(endAge) == false )
    {
        // we assume by default that the end is at time 0
        endAge = 0.0;
    }
    double startAge = endAge + node->getBranchLength();

    double mu = death * rate;
    double beta = rate;

    // first, get the rate matrix for this branch
    const RateGenerator *rg = NULL;

    size_t numMatrices = node_idx+1;
    size_t startMatrix = node_idx;
    if (this->branch_heterogeneous_substitution_matrices == false )
    {
        numMatrices = this->num_matrices;
        startMatrix = 0;
    }

    for (size_t matrix = startMatrix; matrix < numMatrices; ++matrix)
    {
        if ( this->heterogeneous_rate_matrices != NULL )
        {
            rg = &this->heterogeneous_rate_matrices->getValue()[matrix];
        }
        else if ( this->homogeneous_rate_matrix != NULL )
        {
            rg = &this->homogeneous_rate_matrix->getValue();
        }

        const RateMatrix *rm = dynamic_cast<const RateMatrix *>(rg);
        if ( rm == NULL && rg != NULL)
        {
            throw RbException("Dollo model cannot use RateGenerators that are not RateMatrices");
        }

        if (rm != NULL && normalize)
        {
            double avg = rm->averageRate();

            beta *= avg / (avg + death);
            mu /= (avg + death);
        }

        for (size_t i = 0; i < this->num_site_rates; ++i)
        {
            size_t matrix_index = this->branch_heterogeneous_substitution_matrices ? i : i*this->num_matrices + matrix;

            double r = 1.0;
            if ( this->rate_variation_across_sites == true )
                r = this->site_rates->getValue()[i];

            double scaled_mu = r * mu;

            if (node->isRoot())
            {
                integrationFactors[i] = 1.0/scaled_mu;
            }
            else
            {
                survival[i] = exp( - scaled_mu * (startAge - endAge) );
                integrationFactors[i] = (1.0 - survival[i])/scaled_mu;

                if (rm != NULL)
                {
                    rm->calculateTransitionProbabilities( startAge, endAge,  r * beta, this->transition_prob_matrices[matrix_index] );

                    for (size_t a = 0; a < dim; a++)
                    {
                        for (size_t b = 0; b < dim; b++)
                        {
                            this->transition_prob_matrices[matrix_index][a][b] *= survival[i];
                        }

                        this->transition_prob_matrices[matrix_index][a][dim] = 1.0 - survival[i];
                        this->transition_prob_matrices[matrix_index][dim][a] = 0.0;
                    }

                    this->transition_prob_matrices[matrix_index][dim][dim] = 1.0;
                }
                else
                {
                    this->transition_prob_matrices[matrix_index][1][1] = 1.0;
                    this->transition_prob_matrices[matrix_index][1][0] = 0.0;
                    this->transition_prob_matrices[matrix_index][0][1] = 1.0 - survival[i];
                    this->transition_prob_matrices[matrix_index][0][0] = survival[i];
                }
            }
        }
    }
    TransitionProbabilityMatrix m1 = transition_prob_matrices[0];
    TransitionProbabilityMatrix m2 = transition_prob_matrices[1];
    TransitionProbabilityMatrix m3 = transition_prob_matrices[2];
    TransitionProbabilityMatrix m4 = transition_prob_matrices[3];

    return;
}


void RevBayesCore::PhyloCTMCSiteHomogeneousDollo::updateTransitionProbabilityMatrices( void )
{
    // doing nothing here as this function has only been implemented for the classes that
    // do not have its own updateTransitionProbabilities function.
    // making this an empty function for this class for now so no redundant computation would be incurred when
    // updateTransitionProbabilityMatrices gets called in computeLnProbability in AbstractPhyloCTMCSiteHomogeneous
}


void RevBayesCore::PhyloCTMCSiteHomogeneousDollo::getStationaryFrequencies( std::vector<std::vector<double> >& rf ) const
{
    if ( this->branch_heterogeneous_substitution_matrices == false )
    {
        return getRootFrequencies(rf);
    }

    // first, get the rate matrix for this branch
    const RateGenerator *rg = NULL;
    if ( this->heterogeneous_rate_matrices != NULL )
    {
        for (size_t matrix = 0; matrix < num_matrices; matrix++)
        {
            rg = &this->heterogeneous_rate_matrices->getValue()[matrix];

            const RateMatrix *rm = dynamic_cast<const RateMatrix *>(rg);
            if ( rm != NULL )
            {
                rf.push_back( rm->getStationaryFrequencies() );
            }
            else
            {
                throw RbException("Dollo model cannot use RateGenerators that are not RateMatrices");
            }
        }
    }
    else if (this->homogeneous_rate_matrix != NULL)
    {
        rg = &this->homogeneous_rate_matrix->getValue();

        const RateMatrix *rm = dynamic_cast<const RateMatrix *>(rg);
        if ( rm != NULL )
        {
            rf.push_back( rm->getStationaryFrequencies() );
        }
        else
        {
            throw RbException("Dollo model cannot use RateGenerators that are not RateMatrices");
        }
    }
    else
    {
        rf.push_back( std::vector<double>(1,1.0) );
    }

}

void RevBayesCore::PhyloCTMCSiteHomogeneousDollo::setDeathRate(const TypedDagNode< double > *r)
{

    // remove the old parameter first
    if ( death_rate != NULL )
    {
        this->removeParameter( death_rate );
        death_rate = NULL;
    }

    // set the value
    death_rate = r;

    // add the new parameter
    this->addParameter( death_rate );

    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }

}

void RevBayesCore::PhyloCTMCSiteHomogeneousDollo::computeRootLikelihood( size_t root, size_t left, size_t right)
{
    // compute the transition probability matrix
    updateTransitionProbabilities( root );

    // get the root frequencies
    std::vector<std::vector<double> > ff;
    this->getStationaryFrequencies(ff);

    // get the pointers to the partial likelihoods of the left and right subtree
          double* p        = partialLikelihoods + activeLikelihood[root]  * activeLikelihoodOffset + root  * nodeOffset;
    const double* p_left   = partialLikelihoods + activeLikelihood[left]  * activeLikelihoodOffset + left  * nodeOffset;
    const double* p_right  = partialLikelihoods + activeLikelihood[right] * activeLikelihoodOffset + right * nodeOffset;

    // get pointers the likelihood for both subtrees
          double*   p_mixture          = p;
    const double*   p_mixture_left     = p_left;
    const double*   p_mixture_right    = p_right;

    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < num_site_mixtures; ++mixture)
    {
        // get the root frequencies
        const std::vector<double> &f = branch_heterogeneous_substitution_matrices ? ff[root] : ff[mixture % ff.size()];

        // get pointers to the likelihood for this mixture category
              double*   p_site_mixture          = p_mixture;
        const double*   p_site_mixture_left     = p_mixture_left;
        const double*   p_site_mixture_right    = p_mixture_right;
        // iterate over all sites
        for (size_t site = 0; site < pattern_block_size; ++site)
        {
            p_site_mixture[dim + 1] = 0.0;

            p_site_mixture[dim] = p_site_mixture_left[dim] * p_site_mixture_right[dim];

            if (p_site_mixture[dim] == 0)
            {
                for (size_t c = 0; c < dim; c++)
                {
                    p_site_mixture[dim + 1] += p_site_mixture_left[c] * p_site_mixture_right[c] * f[c];
                }

                p_site_mixture[dim + 1] *= integrationFactors[mixture];
            }

            if ( !RbSettings::userSettings().getUseScaling() )
            {
                p_site_mixture[dim + 1] += ( p_site_mixture_left[dim] > 0 ) * p_site_mixture_right[dim + 1] + ( p_site_mixture_right[dim] > 0 ) * p_site_mixture_left[dim + 1];
            }

            // increment the pointers to the next site
            p_site_mixture+=siteOffset; p_site_mixture_left+=siteOffset; p_site_mixture_right+=siteOffset;

        } // end-for over all sites (=patterns)

        // increment the pointers to the next mixture category
        p_mixture+=mixtureOffset; p_mixture_left+=mixtureOffset; p_mixture_right+=mixtureOffset;

    } // end-for over all mixtures (=rate categories)

    computeRootCorrection(root, left, right);
}

void RevBayesCore::PhyloCTMCSiteHomogeneousDollo::computeRootLikelihood( size_t root, size_t left, size_t right, size_t middle)
{
    // compute the transition probability matrix
    updateTransitionProbabilities( root );

    // get the root frequencies
    std::vector<std::vector<double> > ff;
    this->getRootFrequencies(ff);

    // get the pointers to the partial likelihoods of the left and right subtree
          double* p        = partialLikelihoods + activeLikelihood[root]   * activeLikelihoodOffset + root   * nodeOffset;
    const double* p_left   = partialLikelihoods + activeLikelihood[left]   * activeLikelihoodOffset + left   * nodeOffset;
    const double* p_right  = partialLikelihoods + activeLikelihood[right]  * activeLikelihoodOffset + right  * nodeOffset;
    const double* p_middle = partialLikelihoods + activeLikelihood[middle] * activeLikelihoodOffset + middle * nodeOffset;

    // get pointers the likelihood for both subtrees
          double*   p_mixture          = p;
    const double*   p_mixture_left     = p_left;
    const double*   p_mixture_right    = p_right;
    const double*   p_mixture_middle   = p_middle;

    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < num_site_mixtures; ++mixture)
    {
        // get the root frequencies
        const std::vector<double> &f = ff[mixture % ff.size()];

        // get pointers to the likelihood for this mixture category
              double*   p_site_mixture          = p_mixture;
        const double*   p_site_mixture_left     = p_mixture_left;
        const double*   p_site_mixture_right    = p_mixture_right;
        const double*   p_site_mixture_middle   = p_mixture_middle;
        // iterate over all sites
        for (size_t site = 0; site < pattern_block_size; ++site)
        {
            p_site_mixture[dim + 1] = 0.0;

            p_site_mixture[dim] = p_site_mixture_left[dim] * p_site_mixture_right[dim] * p_site_mixture_middle[dim];

            if (p_site_mixture[dim] == 0)
            {
                for (size_t c = 0; c < dim; c++)
                {
                    p_site_mixture[dim + 1] += p_site_mixture_left[c] * p_site_mixture_right[c] * p_site_mixture_middle[c] * f[c];
                }

                p_site_mixture[dim + 1] *= integrationFactors[mixture];
            }

            if ( !RbSettings::userSettings().getUseScaling() )
            {
                p_site_mixture[dim + 1] += ( p_site_mixture_left[dim] > 0 )   * ( p_site_mixture_middle[dim] > 0 ) * p_site_mixture_right[dim + 1]
                                       + ( p_site_mixture_left[dim] > 0 )   * ( p_site_mixture_right[dim] > 0 )  * p_site_mixture_middle[dim + 1]
                                       + ( p_site_mixture_middle[dim] > 0 ) * ( p_site_mixture_right[dim] > 0 )  * p_site_mixture_left[dim + 1];
            }

            // increment the pointers to the next site
            p_site_mixture+=siteOffset; p_site_mixture_left+=siteOffset; p_site_mixture_right+=siteOffset; p_site_mixture_middle+=siteOffset;

        } // end-for over all sites (=patterns)

        // increment the pointers to the next mixture category
        p_mixture+=mixtureOffset; p_mixture_left+=mixtureOffset; p_mixture_right+=mixtureOffset; p_mixture_middle+=mixtureOffset;

    } // end-for over all mixtures (=rate categories)

    computeRootCorrection(root, left, right, middle);
}

void RevBayesCore::PhyloCTMCSiteHomogeneousDollo::computeInternalNodeLikelihood(const TopologyNode &node, size_t node_index, size_t left, size_t right)
{
    // compute the transition probability matrix
    updateTransitionProbabilities( node_index );

    // get the root frequencies
    std::vector<std::vector<double> > ff;
    getStationaryFrequencies(ff);

    // get the pointers to the partial likelihoods for this node and the two descendant subtrees
    const double*   p_left  = partialLikelihoods + activeLikelihood[left]*activeLikelihoodOffset + left*nodeOffset;
    const double*   p_right = partialLikelihoods + activeLikelihood[right]*activeLikelihoodOffset + right*nodeOffset;
    double*         p_node  = partialLikelihoods + activeLikelihood[node_index]*activeLikelihoodOffset + node_index*nodeOffset;

    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < num_site_mixtures; ++mixture)
    {
        // get the root frequencies
        const std::vector<double> &f = branch_heterogeneous_substitution_matrices ? ff[node_index % ff.size()] : ff[mixture % ff.size()];

        // the transition probability matrix for this mixture category
        const TransitionProbabilityMatrix&    pij = this->transition_prob_matrices[mixture];

        // get the pointers to the likelihood for this mixture category
        double*          p_site_mixture          = p_node;
        const double*    p_site_mixture_left     = p_left;
        const double*    p_site_mixture_right    = p_right;
        // compute the per site probabilities
        for (size_t site = 0; site < pattern_block_size ; ++site)
        {
            p_site_mixture[dim] = p_site_mixture_left[dim] * p_site_mixture_right[dim];

            for (size_t from = 0; from < dim; from++)
            {
                p_site_mixture[from] = 0.0;

                for (size_t to = 0; to < dim + 1; to++)
                {
                    p_site_mixture[from] += p_site_mixture_left[to] * p_site_mixture_right[to] * pij[from][to] ;
                }
            }

            p_site_mixture[dim + 1] = 0.0;

            if (p_site_mixture[dim] == 0)
            {
                for (size_t c = 0; c < dim; c++)
                {
                    p_site_mixture[dim + 1] += p_site_mixture_left[c] * p_site_mixture_right[c] * f[c];
                }

                p_site_mixture[dim + 1] *= integrationFactors[mixture];
            }

            if ( !RbSettings::userSettings().getUseScaling() )
            {
                p_site_mixture[dim + 1] += ( p_site_mixture_left[dim] > 0 ) * p_site_mixture_right[dim + 1] + ( p_site_mixture_right[dim] > 0 ) * p_site_mixture_left[dim + 1];
            }

            // increment the pointers to the next site
            p_site_mixture_left+=siteOffset; p_site_mixture_right+=siteOffset; p_site_mixture+=siteOffset;

        } // end-for over all sites (=patterns)

        p_left+=mixtureOffset; p_right+=mixtureOffset; p_node+=mixtureOffset;
    } // end-for over all mixtures (=rate-categories)

    computeInternalNodeCorrection(node, node_index, left, right);
}

void RevBayesCore::PhyloCTMCSiteHomogeneousDollo::computeInternalNodeLikelihood(const TopologyNode &node, size_t node_index, size_t left, size_t right, size_t middle)
{

    // compute the transition probability matrix
    updateTransitionProbabilities( node_index );

    // get the root frequencies
    std::vector<std::vector<double> > ff;
    getStationaryFrequencies(ff);

    // get the pointers to the partial likelihoods for this node and the two descendant subtrees
    const double*   p_left      = partialLikelihoods + activeLikelihood[left]*activeLikelihoodOffset + left*nodeOffset;
    const double*   p_middle    = partialLikelihoods + activeLikelihood[middle]*activeLikelihoodOffset + middle*nodeOffset;
    const double*   p_right     = partialLikelihoods + activeLikelihood[right]*activeLikelihoodOffset + right*nodeOffset;
    double*         p_node      = partialLikelihoods + activeLikelihood[node_index]*activeLikelihoodOffset + node_index*nodeOffset;

    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < num_site_mixtures; ++mixture)
    {
        // get the root frequencies
        const std::vector<double> &f = branch_heterogeneous_substitution_matrices ? ff[node_index % ff.size()] : ff[mixture % ff.size()];

        // the transition probability matrix for this mixture category
        const TransitionProbabilityMatrix&    pij = this->transition_prob_matrices[mixture];

        // get the pointers to the likelihood for this mixture category
        double*          p_site_mixture          = p_node;
        const double*    p_site_mixture_left     = p_left;
        const double*    p_site_mixture_middle   = p_middle;
        const double*    p_site_mixture_right    = p_right;
        // compute the per site probabilities
        for (size_t site = 0; site < pattern_block_size ; ++site)
        {
            for (size_t from = 0; from < dim; from++)
            {
                p_site_mixture[from] = 0.0;

                for (size_t to = 0; to < dim + 1; to++)
                {
                    p_site_mixture[from] += p_site_mixture_left[to] * p_site_mixture_right[to] * p_site_mixture_middle[to] * pij[from][to];
                }
            }

            p_site_mixture[dim] = p_site_mixture_left[dim] * p_site_mixture_right[dim] * p_site_mixture_middle[dim];
            p_site_mixture[dim + 1] = 0.0;

            if (p_site_mixture[dim] == 0)
            {
                for (size_t c = 0; c < dim; c++)
                {
                    p_site_mixture[dim + 1] += p_site_mixture_left[c] * p_site_mixture_right[c] * p_site_mixture_middle[c] * f[c];
                }

                p_site_mixture[dim + 1] *= integrationFactors[mixture];
            }

            if ( !RbSettings::userSettings().getUseScaling() )
            {
                p_site_mixture[dim + 1] += ( p_site_mixture_left[dim] > 0 )   * ( p_site_mixture_middle[dim] > 0 ) * p_site_mixture_right[dim + 1]
                                       + ( p_site_mixture_left[dim] > 0 )   * ( p_site_mixture_right[dim] > 0 )  * p_site_mixture_middle[dim + 1]
                                       + ( p_site_mixture_middle[dim] > 0 ) * ( p_site_mixture_right[dim] > 0 )  * p_site_mixture_left[dim + 1];
            }

            // increment the pointers to the next site
            p_site_mixture_left+=siteOffset; p_site_mixture_middle+=siteOffset; p_site_mixture_right+=siteOffset; p_site_mixture+=siteOffset;

        } // end-for over all sites (=patterns)

        p_left+=mixtureOffset; p_middle+=mixtureOffset; p_right+=mixtureOffset; p_node+=mixtureOffset;

    } // end-for over all mixtures (=rate-categories)

    computeInternalNodeCorrection(node, node_index, left, right, middle);
}

void RevBayesCore::PhyloCTMCSiteHomogeneousDollo::computeTipLikelihood(const TopologyNode &node, size_t node_index)
{

    double* p_node = partialLikelihoods + activeLikelihood[node_index]*activeLikelihoodOffset + node_index*nodeOffset;

    
    size_t data_tip_index = this->taxon_name_2_tip_index_map[ node.getName() ];
    const std::vector<bool> &gap_node = this->gap_matrix[data_tip_index];
    const std::vector<unsigned long> &char_node = this->char_matrix[data_tip_index];
    const std::vector<RbBitSet> &amb_char_node = this->ambiguous_char_matrix[data_tip_index];
    
    // compute the transition probabilities
    updateTransitionProbabilities( node_index );

    // get the root frequencies
    std::vector<std::vector<double> > ff;
    getStationaryFrequencies(ff);

    double*   p_mixture      = p_node;

    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < num_site_mixtures; ++mixture)
    {
        // get the root frequencies
        const std::vector<double> &f = branch_heterogeneous_substitution_matrices ? ff[node_index % ff.size()] : ff[mixture % ff.size()];

        // the transition probability matrix for this mixture category
        const TransitionProbabilityMatrix&    pij = this->transition_prob_matrices[mixture];

        // get the pointer to the likelihoods for this site and mixture category
        double*     p_site_mixture      = p_mixture;

        // iterate over all sites
        for (size_t site = 0; site != pattern_block_size; ++site)
        {

            // is this site a gap?
            if ( gap_node[site] )
            {
                for (size_t c = 0; c < dim + 1; c++)
                {
                    p_site_mixture[c] = 1.0;
                }

                p_site_mixture[dim + 1] = integrationFactors[mixture];
            }
            else // we have observed a character
            {

                if ( this->using_ambiguous_characters == true )
                {
                    // compute the likelihood that we had a transition from state c1 to the observed state org_val
                    // note, the observed state could be ambiguous!
                    const RbBitSet &val = amb_char_node[site];

                    for (size_t c = 0; c < dim + 1; c++)
                    {
                        double tmp = 0.0;

                        for ( size_t i=0; i<val.size(); ++i )
                        {
                            // check whether we observed this state
                            if ( val.test(i) == true )
                            {
                                tmp += pij[c][i];
                            }
                        } // end-while over all observed states for this character

                        // store the likelihood
                        p_site_mixture[c] = tmp;
                    }

                    p_site_mixture[dim + 1] = 0.0;

                    for ( size_t i=1; i<val.size(); ++i )
                    {
                        if ( val.test(i) == true )
                        {
                            // check whether we observed this state
                            p_site_mixture[dim + 1] += f[i];
                        }
                    } // end-while over all observed states for this character

                    p_site_mixture[dim + 1] *= integrationFactors[mixture];

                }
                else
                {
                    // get the original character
                    // shift the characters so that a zero state is the (n-1)th state
                    unsigned long org_val = char_node[site] == 0 ? dim : char_node[site] - 1;

                    // store the branch likelihoods and integrated node likelihood
                    for (size_t c = 0; c < dim + 1; c++)
                    {
                        p_site_mixture[c] = pij[c][org_val];
                    }

                    p_site_mixture[dim + 1] = (org_val == dim) ? 0.0 : f[org_val] * integrationFactors[mixture];
                }


            } // end-if a gap state

            // increment the pointers to next site
            p_site_mixture+=siteOffset;

        } // end-for over all sites/patterns in the sequence

        // increment the pointers to next mixture category
        p_mixture+=mixtureOffset;

    } // end-for over all mixture categories

    computeTipCorrection(node, node_index);
}

void RevBayesCore::PhyloCTMCSiteHomogeneousDollo::computeTipCorrection(const TopologyNode &node, size_t node_index)
{
    std::vector<double>::iterator p_node = correctionLikelihoods.begin() + this->activeLikelihood[node_index]*this->activeCorrectionOffset + node_index*correctionNodeOffset;
    std::vector<double>::iterator c_node = perMaskMixtureCorrections.begin() + this->activeLikelihood[node_index]*activeMassOffset + node_index*massNodeOffset;

    std::vector<std::vector<std::vector<double> > > partialNodeCorrections = std::vector<std::vector<std::vector<double> > >(dim + 1, std::vector<std::vector<double> >(dim + 1, std::vector<double>(numCorrectionPatterns, 0.0)));

    // get the root frequencies
    std::vector<std::vector<double> > ff;
    getStationaryFrequencies(ff);

    // iterate over correction masks
    for (size_t mask = 0; mask < numCorrectionMasks; mask++)
    {
        bool gap = correctionMaskMatrix[mask][node_index];

        maskNodeObservationCounts[mask][node_index] = !gap;

        std::vector<double>::iterator                 c_mask_mixture = c_node   + num_site_mixtures*mask;

        // iterate over all mixture categories
        for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
        {
            const TransitionProbabilityMatrix&    pij = this->transition_prob_matrices[mixture];

            // iterate over ancestral (non-autapomorphic) states
            for (size_t a = 0; a < dim + 1; a++)
            {
                size_t offset = mixture*correctionMixtureOffset + mask*correctionMaskOffset + a*correctionOffset;

                std::vector<double>::iterator                 u = p_node   + offset;

                // iterate over constant/autapomorphic states
                for (size_t i = 0; i < dim + 1; i++)
                {
                    if (i > 0 && numCorrectionPatterns == 1)
                        break;

                    // get bit pattern for this constant/autapomorphic state
                    size_t c = ( i == 0 ? 0 : pow(2.0f, float(i - 1)) );

                    std::vector<double>::iterator         uc = u  + c*(dim + 1);

                    // if we have a gap, then fill with ones for constant patterns
                    // or zero for other patterns
                    if (gap)
                    {
                        std::fill(uc, uc + dim + 1, i == 0);
                    }
                    // if this is a constant pattern, fill with transition prob from ci to a
                    else if (c == 0)
                    {
                        // iterate over initial states
                        for (size_t ci = 0; ci < dim + 1; ci++)
                        {
                            uc[ci] = pij[ci][a];
                        }

                        partialNodeCorrections[a][a][c] = 1.0;
                    }
                    // if this is an autapomorphic state cj, fill with transition prob from ci to cj
                    else
                    {
                        // get autamorphic state
                        size_t cj = (i-1 >= a ? i : i-1);

                        // iterate over initial states
                        for (size_t ci = 0; ci < dim + 1; ci++)
                        {
                            uc[ci] = pij[ci][cj];
                        }

                        partialNodeCorrections[cj][a][c] = 1.0;
                    }
                }
            }

            c_mask_mixture[mixture] = computeIntegratedNodeCorrection(partialNodeCorrections, node_index, mask, mixture, branch_heterogeneous_substitution_matrices ? ff[node_index % ff.size()] : ff[mixture % ff.size()]);
        }
    }
}

void RevBayesCore::PhyloCTMCSiteHomogeneousDollo::computeInternalNodeCorrection(const TopologyNode &node, size_t node_index, size_t left, size_t right, size_t middle)
{
    // get the pointers to the partial likelihoods for this node and the two descendant subtrees
    std::vector<double>::const_iterator   p_left   = correctionLikelihoods.begin() + this->activeLikelihood[left]*activeCorrectionOffset + left*correctionNodeOffset;
    std::vector<double>::const_iterator   p_right  = correctionLikelihoods.begin() + this->activeLikelihood[right]*activeCorrectionOffset + right*correctionNodeOffset;
    std::vector<double>::const_iterator   p_middle = correctionLikelihoods.begin() + this->activeLikelihood[middle]*activeCorrectionOffset + middle*correctionNodeOffset;
    std::vector<double>::iterator         p_node   = correctionLikelihoods.begin() + this->activeLikelihood[node_index]*activeCorrectionOffset + node_index*correctionNodeOffset;

    std::vector<double>::iterator c_node   = perMaskMixtureCorrections.begin() + this->activeLikelihood[node_index]*activeMassOffset + node_index*massNodeOffset;

    std::vector<std::vector<std::vector<double> > > partialNodeCorrections = std::vector<std::vector<std::vector<double> > >(dim + 1, std::vector<std::vector<double> >(dim + 1, std::vector<double>(numCorrectionPatterns, 0.0)));

    // get the root frequencies
    std::vector<std::vector<double> > ff;
    getStationaryFrequencies(ff);

    // iterate over correction masks
    for (size_t mask = 0; mask < numCorrectionMasks; mask++)
    {
        maskNodeObservationCounts[mask][node_index] = maskNodeObservationCounts[mask][left] + maskNodeObservationCounts[mask][right] + maskNodeObservationCounts[mask][middle];

        std::vector<double>::iterator                 c_mask_mixture = c_node   + num_site_mixtures*mask;

        // iterate over all mixture categories
        for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
        {
            const TransitionProbabilityMatrix&    pij = this->transition_prob_matrices[mixture];

            // iterate over ancestral (non-autapomorphic) states
            for (size_t a = 0; a < dim + 1; a++)
            {
                size_t offset = mixture*correctionMixtureOffset + mask*correctionMaskOffset + a*correctionOffset;

                std::vector<double>::iterator                 u = p_node   + offset;
                std::vector<double>::const_iterator         u_l = p_left   + offset;
                std::vector<double>::const_iterator         u_r = p_right  + offset;
                std::vector<double>::const_iterator         u_m = p_middle + offset;

                // iterate over combinations of autapomorphic states
                for (size_t c = 0; c < numCorrectionPatterns; c++)
                {
                    std::vector<double>::iterator         uc = u  + c*(dim + 1);

                    std::fill(uc, uc + dim + 1, 0.0);

                    // iterate over partitions of c
                    for (size_t p1 = 0; p1 <= c; p1++)
                    {
                        if ( (p1 | c) == c)
                        {
                            size_t p_tmp = p1 ^ c;

                            // iterate over partitions of p_tmp
                            for (size_t p2 = 0; p2 <= p_tmp; p2++)
                            {
                                if ( (p2 | p_tmp) == p_tmp)
                                {
                                    size_t p3 = p2 ^ p_tmp;

                                    std::vector<double>::const_iterator         lc = u_l  + p1*(dim + 1);
                                    std::vector<double>::const_iterator         rc = u_r  + p2*(dim + 1);
                                    std::vector<double>::const_iterator         mc = u_m  + p3*(dim + 1);

                                    // iterate over ending states
                                    for (size_t cj = 0; cj < dim + 1; cj++)
                                    {
                                        double pj = lc[cj] * rc[cj] * mc[cj];

                                        partialNodeCorrections[cj][a][c] += pj;

                                        // iterate over initial states
                                        for (size_t ci = 0; ci < dim + 1; ci++)
                                        {
                                            uc[ci] += pij[ci][cj] * pj;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            c_mask_mixture[mixture] = computeIntegratedNodeCorrection(partialNodeCorrections, node_index, mask, mixture, branch_heterogeneous_substitution_matrices ? ff[node_index % ff.size()] : ff[mixture % ff.size()]);
        }
    }
}

void RevBayesCore::PhyloCTMCSiteHomogeneousDollo::computeInternalNodeCorrection(const TopologyNode &node, size_t node_index, size_t left, size_t right)
{
    // get the pointers to the partial likelihoods for this node and the two descendant subtrees
    std::vector<double>::const_iterator   p_left   = correctionLikelihoods.begin() + this->activeLikelihood[left]*activeCorrectionOffset + left*correctionNodeOffset;
    std::vector<double>::const_iterator   p_right  = correctionLikelihoods.begin() + this->activeLikelihood[right]*activeCorrectionOffset + right*correctionNodeOffset;
    std::vector<double>::iterator         p_node   = correctionLikelihoods.begin() + this->activeLikelihood[node_index]*activeCorrectionOffset + node_index*correctionNodeOffset;

    std::vector<double>::iterator c_node   = perMaskMixtureCorrections.begin() + this->activeLikelihood[node_index]*activeMassOffset + node_index*massNodeOffset;

    std::vector<std::vector<std::vector<double> > > partialNodeCorrections = std::vector<std::vector<std::vector<double> > >(dim + 1, std::vector<std::vector<double> >(dim + 1, std::vector<double>(numCorrectionPatterns, 0.0)));

    // get the root frequencies
    std::vector<std::vector<double> > ff;
    getStationaryFrequencies(ff);

    // iterate over correction masks
    for (size_t mask = 0; mask < numCorrectionMasks; mask++)
    {
        maskNodeObservationCounts[mask][node_index] = maskNodeObservationCounts[mask][left] + maskNodeObservationCounts[mask][right];

        std::vector<double>::iterator                 c_mask_mixture = c_node   + num_site_mixtures*mask;

        // iterate over all mixture categories
        for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
        {
            const TransitionProbabilityMatrix&    pij = this->transition_prob_matrices[mixture];

            // iterate over ancestral (non-autapomorphic) states
            for (size_t a = 0; a < dim + 1; a++)
            {
                size_t offset = mixture*correctionMixtureOffset + mask*correctionMaskOffset + a*correctionOffset;

                std::vector<double>::iterator                 u = p_node   + offset;
                std::vector<double>::const_iterator         u_l = p_left   + offset;
                std::vector<double>::const_iterator         u_r = p_right  + offset;

                // iterate over combinations of autapomorphic states
                for (size_t c = 0; c < numCorrectionPatterns; c++)
                {
                    std::vector<double>::iterator         uc = u  + c*(dim + 1);

                    std::fill(uc, uc + dim + 1, 0.0);

                    // iterate over partitions of c
                    for (size_t p1 = 0; p1 <= c; p1++)
                    {
                        if ( (p1 | c) == c)
                        {
                            size_t p2 = p1 ^ c;

                            std::vector<double>::const_iterator         lc = u_l  + p1*(dim + 1);
                            std::vector<double>::const_iterator         rc = u_r  + p2*(dim + 1);

                            // iterate over ending states
                            for (size_t cj = 0; cj < dim + 1; cj++)
                            {
                                double pj = lc[cj] * rc[cj];

                                partialNodeCorrections[cj][a][c] += pj;

                                // iterate over initial states
                                for (size_t ci = 0; ci < dim + 1; ci++)
                                {
                                    uc[ci] += pij[ci][cj] * pj;
                                }
                            }
                        }
                    }
                }
            }

            c_mask_mixture[mixture] = computeIntegratedNodeCorrection(partialNodeCorrections, node_index, mask, mixture, branch_heterogeneous_substitution_matrices ? ff[node_index % ff.size()] : ff[mixture % ff.size()]);
        }
    }
}

void RevBayesCore::PhyloCTMCSiteHomogeneousDollo::computeRootCorrection( size_t root, size_t left, size_t right, size_t middle)
{
    // get the pointers to the partial likelihoods for this node and the two descendant subtrees
    std::vector<double>::iterator         p_node   = correctionLikelihoods.begin() + this->activeLikelihood[root]*activeCorrectionOffset + root*correctionNodeOffset;
    std::vector<double>::const_iterator   p_left   = correctionLikelihoods.begin() + this->activeLikelihood[left]*activeCorrectionOffset + left*correctionNodeOffset;
    std::vector<double>::const_iterator   p_right  = correctionLikelihoods.begin() + this->activeLikelihood[right]*activeCorrectionOffset + right*correctionNodeOffset;
    std::vector<double>::const_iterator   p_middle = correctionLikelihoods.begin() + this->activeLikelihood[middle]*activeCorrectionOffset + middle*correctionNodeOffset;

    std::vector<double>::iterator c_node   = perMaskMixtureCorrections.begin() + this->activeLikelihood[root]*activeMassOffset + root*massNodeOffset;

    std::vector<std::vector<std::vector<double> > > partialNodeCorrections = std::vector<std::vector<std::vector<double> > >(dim + 1, std::vector<std::vector<double> >(dim + 1, std::vector<double>(numCorrectionPatterns, 0.0)));

    // get the root frequencies
    std::vector<std::vector<double> > ff;
    getStationaryFrequencies(ff);

    // iterate over correction masks
    for (size_t mask = 0; mask < numCorrectionMasks; mask++)
    {
        maskNodeObservationCounts[mask][root] = maskNodeObservationCounts[mask][left] + maskNodeObservationCounts[mask][right] + maskNodeObservationCounts[mask][middle];

        std::vector<double>::iterator                 c_mask_mixture = c_node   + num_site_mixtures*mask;

        // iterate over all mixture categories
        for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
        {
            const std::vector<double> &f = branch_heterogeneous_substitution_matrices ? ff[root % ff.size()] : ff[mixture % ff.size()];

            // iterate over ancestral (non-autapomorphic) states
            for (size_t a = 0; a < dim + 1; a++)
            {
                size_t offset = mixture*correctionMixtureOffset + mask*correctionMaskOffset + a*correctionOffset;

                std::vector<double>::iterator                 u = p_node   + offset;
                std::vector<double>::const_iterator         u_l = p_left   + offset;
                std::vector<double>::const_iterator         u_r = p_right  + offset;
                std::vector<double>::const_iterator         u_m = p_middle + offset;

                // iterate over combinations of autapomorphic states
                for (size_t c = 0; c < numCorrectionPatterns; c++)
                {
                    std::vector<double>::iterator         uc = u  + c*(dim + 1);

                    std::fill(uc, uc + dim + 1, 0.0);

                    // iterate over partitions of c
                    for (size_t p1 = 0; p1 <= c; p1++)
                    {
                        if ( (p1 | c) == c)
                        {
                            size_t p_tmp = p1 ^ c;

                            // iterate over partitions of p_tmp
                            for (size_t p2 = 0; p2 <= p_tmp; p2++)
                            {
                                if ( (p2 | p_tmp) == p_tmp)
                                {
                                    size_t p3 = p2 ^ p_tmp;

                                    std::vector<double>::const_iterator         lc = u_l  + p1*(dim + 1);
                                    std::vector<double>::const_iterator         rc = u_r  + p2*(dim + 1);
                                    std::vector<double>::const_iterator         mc = u_m  + p3*(dim + 1);

                                    // iterate over ending states
                                    for (size_t cj = 0; cj < dim; cj++)
                                    {
                                        double pj = lc[cj] * rc[cj] * mc[cj];

                                        partialNodeCorrections[cj][a][c] += pj;

                                        uc[cj] += f[cj] * pj;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            c_mask_mixture[mixture] = computeIntegratedNodeCorrection(partialNodeCorrections, root, mask, mixture, f);
        }
    }
}

void RevBayesCore::PhyloCTMCSiteHomogeneousDollo::computeRootCorrection( size_t root, size_t left, size_t right)
{
    // get the pointers to the partial likelihoods for this node and the two descendant subtrees
    std::vector<double>::iterator         p_node   = correctionLikelihoods.begin() + this->activeLikelihood[root]*activeCorrectionOffset + root*correctionNodeOffset;
    std::vector<double>::const_iterator   p_left   = correctionLikelihoods.begin() + this->activeLikelihood[left]*activeCorrectionOffset + left*correctionNodeOffset;
    std::vector<double>::const_iterator   p_right  = correctionLikelihoods.begin() + this->activeLikelihood[right]*activeCorrectionOffset + right*correctionNodeOffset;

    std::vector<double>::iterator c_node   = perMaskMixtureCorrections.begin() + this->activeLikelihood[root]*activeMassOffset + root*massNodeOffset;

    std::vector<std::vector<std::vector<double> > > partialNodeCorrections = std::vector<std::vector<std::vector<double> > >(dim + 1, std::vector<std::vector<double> >(dim + 1, std::vector<double>(numCorrectionPatterns, 0.0)));

    // get the root frequencies
    std::vector<std::vector<double> > ff;
    getStationaryFrequencies(ff);

    // iterate over correction masks
    for (size_t mask = 0; mask < numCorrectionMasks; mask++)
    {
        maskNodeObservationCounts[mask][root] = maskNodeObservationCounts[mask][left] + maskNodeObservationCounts[mask][right];

        std::vector<double>::iterator                 c_mask_mixture = c_node   + num_site_mixtures*mask;

        // iterate over all mixture categories
        for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
        {
            const std::vector<double> &f = branch_heterogeneous_substitution_matrices ? ff[root % ff.size()] : ff[mixture % ff.size()];

            // iterate over ancestral (non-autapomorphic) states
            for (size_t a = 0; a < dim + 1; a++)
            {
                size_t offset = mixture*correctionMixtureOffset + mask*correctionMaskOffset + a*correctionOffset;

                std::vector<double>::iterator                 u = p_node   + offset;
                std::vector<double>::const_iterator         u_l = p_left   + offset;
                std::vector<double>::const_iterator         u_r = p_right  + offset;

                // iterate over combinations of autapomorphic states
                for (size_t c = 0; c < numCorrectionPatterns; c++)
                {
                    std::vector<double>::iterator         uc = u  + c*(dim + 1);

                    std::fill(uc, uc + dim + 1, 0.0);

                    // iterate over partitions of c
                    for (size_t p1 = 0; p1 <= c; p1++)
                    {
                        if ( (p1 | c) == c)
                        {
                            size_t p2 = p1 ^ c;

                            // iterate over partitions of p_tmp
                            std::vector<double>::const_iterator         lc = u_l  + p1*(dim + 1);
                            std::vector<double>::const_iterator         rc = u_r  + p2*(dim + 1);

                            // iterate over ending states
                            for (size_t cj = 0; cj < dim; cj++)
                            {
                                double pj = lc[cj] * rc[cj];

                                partialNodeCorrections[cj][a][c] += pj;

                                uc[cj] += f[cj] * pj;
                            }
                        }
                    }
                }
            }

            c_mask_mixture[mixture] = computeIntegratedNodeCorrection(partialNodeCorrections, root, mask, mixture, f);
        }
    }
}


double RevBayesCore::PhyloCTMCSiteHomogeneousDollo::computeIntegratedNodeCorrection(const std::vector<std::vector<std::vector<double> > >& partials, size_t node_index, size_t mask, size_t mixture, const std::vector<double> &f)
{
    /*
     * partials[a][b][c] contains the probability of observing:
     * autapomorphy pattern c
     * on background state b
     * conditioned on node state a
     */

    double prob = 0.0;

    // iterate over background (non-autapomorphic) states
    for (size_t a = 0; a < dim + 1; a++)
    {
        // iterate over combinations of autapomorphic states
        for (size_t c = 0; c < numCorrectionPatterns; c++)
        {
            if (c == 0 && (coding == DolloAscertainmentBias::VARIABLE || (a == dim && coding == DolloAscertainmentBias::NOABSENCESITES) ))
            {
                // iterate over initial states
                for (size_t ci = 0; ci < dim; ci++)
                {
                    prob += partials[ci][a][c] * f[ci];
                }
            }
            else if (coding == DolloAscertainmentBias::INFORMATIVE)
            {
                // get number of autapomorphic states
                size_t popcount = std::bitset<sizeof(size_t)*CHAR_BIT>(c).count();

               /* if the number of observations is <= the number of autapomorphies then
                * 1. don't count patterns with zero prob (num obs < num auto)
                * 2. don't double-count patterns (num obs == num auto)
                *
                * if num obs == num auto + 1
                * then don't double count the all-autapomorphies pattern
                */
                if ( (maskNodeObservationCounts[mask][node_index] == maskObservationCounts[mask] &&
                    ((maskNodeObservationCounts[mask][node_index] == popcount + 1 && a == 0) ||
                     maskNodeObservationCounts[mask][node_index] > popcount + 1))
                        ||
                    (maskNodeObservationCounts[mask][node_index] != maskObservationCounts[mask] &&
                    ((maskNodeObservationCounts[mask][node_index] == popcount + 1 && a == 0) ||
                     (maskNodeObservationCounts[mask][node_index] > popcount + 1 && a == dim)))
                     )
                {
                    // iterate over initial states
                    for (size_t ci = 0; ci < dim; ci++)
                    {
                        prob += partials[ci][a][c] * f[ci];
                    }
                }
            }
        }
    }

    // impose a per-mixture boundary
    if (prob < 0.0 || prob > 1.0)
    {
        prob = RbConstants::Double::nan;
    }

    if (prob == 0.0 && !this->tau->getValue().getNode(node_index).isTip())
    {
        //prob = RbConstants::Double::nan;
    }

    prob = integrationFactors[mixture]*(1.0 - prob);

    return prob;
}

double RevBayesCore::PhyloCTMCSiteHomogeneousDollo::sumRootLikelihood( void )
{
    // get the root node
    const TopologyNode &root = tau->getValue().getRoot();

    // get the index of the root node
    size_t root_index = root.getIndex();

    const double*   p_root  = this->partialLikelihoods + this->activeLikelihood[root_index] * this->activeLikelihoodOffset + root_index*nodeOffset;

    // create a vector for the per mixture likelihoods
    // we need this vector to sum over the different mixture likelihoods
    std::vector<double> per_mixture_Likelihoods = std::vector<double>(pattern_block_size,0.0);

    const double*   p_site_root = p_root;

    // iterate over all mixture categories
    for (size_t site = 0; site < pattern_block_size; ++site)
    {
        if ( RbSettings::userSettings().getUseScaling() )
        {
            std::vector<double> integrated_likelihoods;

            double max = getScaledNodeWeights(root, site, integrated_likelihoods);

            per_mixture_Likelihoods[site] = RbMath::log_sum_exp(integrated_likelihoods, max);
        }
        else
        {
            const double*   p_site_mixture_root = p_site_root;

            for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
            {
                per_mixture_Likelihoods[site] += p_site_mixture_root[dim + 1];

                p_site_mixture_root += mixtureOffset;
            } // end-for over all sites (=patterns)

            p_site_root += siteOffset;
        }
    } // end-for over all mixtures (=rate categories)

    // sum the log-likelihoods for all sites together
    double sumPartialProbs = 0.0;

    std::vector< size_t >::const_iterator patterns = this->pattern_counts.begin();
    for (size_t site = 0; site < pattern_block_size; ++site, ++patterns)
    {
        if ( RbSettings::userSettings().getUseScaling() )
        {
            sumPartialProbs += ( per_mixture_Likelihoods[site] - log(this->num_site_mixtures ) ) * *patterns;
        }
        else
        {
            sumPartialProbs += log( per_mixture_Likelihoods[site] / this->num_site_mixtures ) * *patterns;
        }

    }


#ifdef RB_MPI

    // we only need to send message if there is more than one process
    if ( num_processes > 1 )
    {

        // send the likelihood from the helpers to the master
        if ( process_active == false )
        {
            // send from the workers the log-likelihood to the master
            MPI_Send(&sumPartialProbs, 1, MPI_DOUBLE, int(active_PID), 0, MPI_COMM_WORLD);
        }

        // receive the likelihoods from the helpers
        if ( process_active == true )
        {
            for (size_t i=active_PID+1; i<active_PID+num_processes; ++i)
            {
                double tmp = 0;
                MPI_Status status;
                MPI_Recv(&tmp, 1, MPI_DOUBLE, int(i), 0, MPI_COMM_WORLD, &status);
                
                sumPartialProbs += tmp;
            }
        }

        // now send back the combined likelihood to the helpers
        if ( process_active == true )
        {
            for (size_t i=active_PID+1; i<active_PID+num_processes; ++i)
            {
                MPI_Send(&sumPartialProbs, 1, MPI_DOUBLE, int(i), 0, MPI_COMM_WORLD);
            }
        }
        else
        {
            MPI_Status status;
            MPI_Recv(&sumPartialProbs, 1, MPI_DOUBLE, int(active_PID), 0, MPI_COMM_WORLD, &status);
        }
    
    }

#endif

    std::vector<double> perMaskCorrections = std::vector<double>(numCorrectionMasks, 0.0);

    // iterate over each correction mask
    for (size_t mask = 0; mask < numCorrectionMasks; mask++)
    {
        for (size_t node = 0; node < num_nodes; ++node)
        {
            std::vector<double>::iterator c_node_mixture   = perMaskMixtureCorrections.begin() + this->activeLikelihood[node]*activeMassOffset + node*massNodeOffset;

            // iterate over all mixture categories
            for (size_t mixture = 0; mixture < num_site_mixtures; ++mixture)
            {
                perMaskCorrections[mask] += c_node_mixture[mixture];
            }
        }
        
        if (perMaskCorrections[mask] <= 0.0)
            perMaskCorrections[mask] = RbConstants::Double::nan;

        perMaskCorrections[mask] = log(perMaskCorrections[mask]) - log(num_site_mixtures);

        // apply the correction for this correction mask
        sumPartialProbs -= perMaskCorrections[mask]*correctionMaskCounts[mask];
    }

    // num_sites is poisson
    sumPartialProbs -= log(num_sites);

    return sumPartialProbs;
}

double RevBayesCore::PhyloCTMCSiteHomogeneousDollo::getScaledNodeWeights(const TopologyNode &node, size_t pattern, std::vector<double>& weights)
{
    //reset the per site likelihood
    double max = RbConstants::Double::neginf;

    size_t node_index = node.getIndex();

    const double* p_node  = partialLikelihoods + activeLikelihood[node_index] * activeLikelihoodOffset  + node_index*nodeOffset + pattern*siteOffset;

    double logScalingFactor = perNodeSiteLogScalingFactors[activeLikelihood[node_index]][node_index][pattern];

    //otherwise, it is an ancestral node so we add the integrated likelihood
    for (size_t mixture = 0; mixture < num_site_mixtures; ++mixture)
    {
        double prob = log(p_node[dim + 1]) - logScalingFactor;

        max = std::max(prob, max);

        weights.push_back( prob );

        p_node += mixtureOffset;
    }

    if (node.isTip())
        return max;

    std::vector<TopologyNode*> children = node.getChildren();

    TopologyNode* child = NULL;
    size_t num_with_descendants = 0;

    for (size_t i = 0; i < children.size(); i++)
    {
        size_t child_index = children[i]->getIndex();
        const double* p_child  = partialLikelihoods + activeLikelihood[child_index]  * activeLikelihoodOffset + child_index*nodeOffset  + pattern*siteOffset;

        // does this child have descendants?
        if (p_child[dim] == 0)
        {
            child = children[i];
            num_with_descendants++;
        }
    }

    //if more than one child has observed descendants, then this is the last ancestral node
    if (num_with_descendants == 1)
    {
        max = std::max(max, getScaledNodeWeights(*child, pattern, weights));
    }

    return max;
}

void RevBayesCore::PhyloCTMCSiteHomogeneousDollo::scale( size_t node_index)
{
    double* p_node = this->partialLikelihoods + this->activeLikelihood[node_index]*this->activeLikelihoodOffset + node_index*this->nodeOffset;

    if ( RbSettings::userSettings().getUseScaling() == true && node_index % RbSettings::userSettings().getScalingDensity() == 0 )
    {
        // iterate over all mixture categories
        for (size_t site = 0; site < this->pattern_block_size ; ++site)
        {

            // the max probability
            double max = 0.0;

            // compute the per site probabilities
            for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
            {
                // get the pointers to the likelihood for this mixture category
                size_t offset = mixture*this->mixtureOffset + site*this->siteOffset;

                double*          p_site_mixture          = p_node + offset;

                for ( size_t i=0; i<dim; ++i)
                {
                    if ( p_site_mixture[i] > max )
                    {
                        max = p_site_mixture[i];
                    }
                }

            }

            this->perNodeSiteLogScalingFactors[this->activeLikelihood[node_index]][node_index][site] = -log(max);


            // compute the per site probabilities
            for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
            {
                // get the pointers to the likelihood for this mixture category
                size_t offset = mixture*this->mixtureOffset + site*this->siteOffset;

                double*          p_site_mixture          = p_node + offset;

                for ( size_t i=0; i<=dim + 1; ++i)
                {
                    p_site_mixture[i] /= max;
                }

            }

        }
    }
    else if ( RbSettings::userSettings().getUseScaling() == true )
    {
        // iterate over all mixture categories
        for (size_t site = 0; site < this->pattern_block_size ; ++site)
        {
            this->perNodeSiteLogScalingFactors[this->activeLikelihood[node_index]][node_index][site] = 0;
        }

    }
}


void RevBayesCore::PhyloCTMCSiteHomogeneousDollo::scale( size_t node_index, size_t left, size_t right )
{
    double* p_node = this->partialLikelihoods + this->activeLikelihood[node_index]*this->activeLikelihoodOffset + node_index*this->nodeOffset;

    if ( RbSettings::userSettings().getUseScaling() == true && node_index % RbSettings::userSettings().getScalingDensity() == 0 && node_index < num_nodes -1)
    {
        // iterate over all mixture categories
        for (size_t site = 0; site < this->pattern_block_size ; ++site)
        {

            // the max probability
            double max = 0.0;

            // compute the per site probabilities
            for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
            {
                // get the pointers to the likelihood for this mixture category
                size_t offset = mixture*this->mixtureOffset + site*this->siteOffset;

                double*          p_site_mixture          = p_node + offset;

                for ( size_t i=0; i<dim; ++i)
                {
                    if ( p_site_mixture[i] > max )
                    {
                        max = p_site_mixture[i];
                    }
                }

            }

            this->perNodeSiteLogScalingFactors[this->activeLikelihood[node_index]][node_index][site] = this->perNodeSiteLogScalingFactors[this->activeLikelihood[left]][left][site] + this->perNodeSiteLogScalingFactors[this->activeLikelihood[right]][right][site] - log(max);


            // compute the per site probabilities
            for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
            {
                // get the pointers to the likelihood for this mixture category
                size_t offset = mixture*this->mixtureOffset + site*this->siteOffset;

                double*          p_site_mixture          = p_node + offset;

                for ( size_t i=0; i<=dim + 1; ++i)
                {
                    p_site_mixture[i] /= max;
                }

            }

        }
    }
    else if ( RbSettings::userSettings().getUseScaling() == true )
    {
        // iterate over all mixture categories
        for (size_t site = 0; site < this->pattern_block_size ; ++site)
        {
            this->perNodeSiteLogScalingFactors[this->activeLikelihood[node_index]][node_index][site] = this->perNodeSiteLogScalingFactors[this->activeLikelihood[left]][left][site] + this->perNodeSiteLogScalingFactors[this->activeLikelihood[right]][right][site];
        }

    }
}

void RevBayesCore::PhyloCTMCSiteHomogeneousDollo::scale( size_t node_index, size_t left, size_t right, size_t middle )
{
    double* p_node   = this->partialLikelihoods + this->activeLikelihood[node_index]*this->activeLikelihoodOffset + node_index*this->nodeOffset;

    if ( RbSettings::userSettings().getUseScaling() == true && node_index % RbSettings::userSettings().getScalingDensity() == 0 && node_index < num_nodes -1)
    {
        // iterate over all mixture categories
        for (size_t site = 0; site < this->pattern_block_size ; ++site)
        {

            // the max probability
            double max = 0.0;

            // compute the per site probabilities
            for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
            {
                // get the pointers to the likelihood for this mixture category
                size_t offset = mixture*this->mixtureOffset + site*this->siteOffset;

                double*          p_site_mixture          = p_node + offset;

                for ( size_t i=0; i<dim; ++i)
                {
                    if ( p_site_mixture[i] > max )
                    {
                        max = p_site_mixture[i];
                    }
                }

            }

            this->perNodeSiteLogScalingFactors[this->activeLikelihood[node_index]][node_index][site] = this->perNodeSiteLogScalingFactors[this->activeLikelihood[left]][left][site] + this->perNodeSiteLogScalingFactors[this->activeLikelihood[right]][right][site] + this->perNodeSiteLogScalingFactors[this->activeLikelihood[middle]][middle][site] - log(max);


            // compute the per site probabilities
            for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
            {
                // get the pointers to the likelihood for this mixture category
                size_t offset = mixture*this->mixtureOffset + site*this->siteOffset;

                double*          p_site_mixture          = p_node + offset;

                for ( size_t i=0; i<=dim + 1; ++i)
                {
                    p_site_mixture[i] /= max;
                }

            }

        }
    }
    else if ( RbSettings::userSettings().getUseScaling() == true )
    {
        // iterate over all mixture categories
        for (size_t site = 0; site < this->pattern_block_size ; ++site)
        {
            this->perNodeSiteLogScalingFactors[this->activeLikelihood[node_index]][node_index][site] = this->perNodeSiteLogScalingFactors[this->activeLikelihood[left]][left][site] + this->perNodeSiteLogScalingFactors[this->activeLikelihood[right]][right][site] + this->perNodeSiteLogScalingFactors[this->activeLikelihood[middle]][middle][site];
        }

    }
}

/** Swap a parameter of the distribution */
void RevBayesCore::PhyloCTMCSiteHomogeneousDollo::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{

    if (oldP == death_rate)
    {
        death_rate = static_cast<const TypedDagNode< double >* >( newP );
    }
    else
    {
        PhyloCTMCSiteHomogeneousConditional<StandardState>::swapParameterInternal(oldP,newP);
    }

}

void RevBayesCore::PhyloCTMCSiteHomogeneousDollo::redrawValue( void ) {

    if ( num_sites == 0 )
        return;

    // delete the old value first
    delete this->value;

    // create a new character data object
    this->value = new HomologousDiscreteCharacterData<StandardState>();

    size_t numTips = this->tau->getValue().getNumberOfTips();
    size_t num_nodes = this->tau->getValue().getNumberOfNodes();

    RandomNumberGenerator* rng = GLOBAL_RNG;

    const TopologyNode &root = this->tau->getValue().getRoot();
    size_t rootIndex = this->tau->getValue().getRoot().getIndex();

    updateCorrections(root, rootIndex);

    std::vector< DiscreteTaxonData<StandardState> > taxa = std::vector< DiscreteTaxonData<StandardState> >(numTips, DiscreteTaxonData<StandardState>( Taxon("") ) );

    std::vector<size_t> rateIndices;
    std::vector<size_t> birthNodes;

    double total = 0.0;
    for ( size_t i = 0; i < this->num_nodes; ++i )
        for ( size_t j = 0; j < this->num_site_mixtures; ++j )
            total += perMaskMixtureCorrections[i*massNodeOffset + j];

    for ( size_t i = 0; i < num_sites; i++ )
    {
        double u = rng->uniform01()*total;
        double tmp = 0.0;

        // simulate a birth and rate index for this character
        // by sampling nodes in proportion to the per node survival probs
        size_t birthNode = 0;
        size_t rateIndex = 0;
        while(tmp < u)
        {
            tmp += perMaskMixtureCorrections[birthNode*massNodeOffset + rateIndex];

            if (tmp < u)
            {
                rateIndex++;
                if (rateIndex == num_site_mixtures)
                {
                    rateIndex = 0;
                    birthNode++;
                }
            }
        }

        rateIndices.push_back(rateIndex);
        birthNodes.push_back(birthNode);
    }

    // get root frequencies
    std::vector< std::vector< double > > ff;
    this->getStationaryFrequencies(ff);

    // then sample site-patterns using rejection sampling,
    // rejecting those that match the unobservable ones.
    for ( size_t i = 0; i < this->num_sites; i++ )
    {
        size_t rateIndex = rateIndices[i];
        size_t birthIndex = birthNodes[i];

        const std::vector< double > &freqs = branch_heterogeneous_substitution_matrices ? ff[birthIndex % ff.size()] : ff[rateIndex % ff.size()];

        std::vector<StandardState> siteData(num_nodes, StandardState(dim));

        // create the character
        StandardState &c = siteData[birthIndex];
        c.setStateByIndex(1);
        // draw the state
        double u = rng->uniform01();
        size_t stateIndex = 0;
        while ( true )
        {
            u -= freqs[stateIndex];
            ++stateIndex;

            if ( u > 0.0 && stateIndex < dim)
            {
                ++c;
            }
            else
            {
                break;
            }

        }

        const TopologyNode& birthNode = tau->getValue().getNode(birthIndex);

        // recursively simulate the sequences
        std::map<size_t, size_t> charCounts;
        simulate( birthNode, siteData, rateIndex, charCounts);

        if ( !isSitePatternCompatible(charCounts) )
        {
            i--;
            continue;
        }

        // add the taxon data to the character data
        for (size_t t = 0; t < numTips; ++t)
        {
            taxa[t].addCharacter(siteData[t]);
        }
    }

    // add the taxon data to the character data
    for (size_t i = 0; i < this->tau->getValue().getNumberOfTips(); ++i)
    {
        taxa[i].setTaxon( this->tau->getValue().getNode(i).getTaxon() );
        this->value->addTaxonData( taxa[i] );
    }

    // compress the data and initialize internal variables
    this->reInitialized();

    for (std::vector<bool>::iterator it = this->dirty_nodes.begin(); it != this->dirty_nodes.end(); ++it)
    {
        (*it) = true;
    }

    // flip the active likelihood pointers
    for (size_t index = 0; index < this->changed_nodes.size(); ++index)
    {
        if ( this->changed_nodes[index] == false )
        {
            this->activeLikelihood[index] = (this->activeLikelihood[index] == 0 ? 1 : 0);
            this->changed_nodes[index] = true;
        }
    }

}

void RevBayesCore::PhyloCTMCSiteHomogeneousDollo::simulate( const TopologyNode &node, std::vector<StandardState> &data, size_t rateIndex, std::map<size_t, size_t>& charCounts) {

    // get the children of the node
    const std::vector<TopologyNode*>& children = node.getChildren();

    // get the sequence of this node
    size_t node_index = node.getIndex();
    StandardState &parentState = data[ node_index ];

    // simulate the sequence for each child
    RandomNumberGenerator* rng = GLOBAL_RNG;
    for (std::vector< TopologyNode* >::const_iterator it = children.begin(); it != children.end(); ++it)
    {
        const TopologyNode &child = *(*it);

        // update the transition probability matrix
        this->updateTransitionProbabilities( child.getIndex() );

        double u = rng->uniform01();

        if (u < survival[rateIndex])
        {
            unsigned long cp = parentState.getStateIndex() - 1;

            double *freqs = this->transition_prob_matrices[ rateIndex ][ cp ];

            // create the character
            StandardState &c = data[ child.getIndex() ];
            c.setStateByIndex(1);
            // draw the state
            double u = rng->uniform01();
            size_t stateIndex = 0;
            while ( true )
            {
                u -= *freqs;
                ++stateIndex;

                if ( u > 0.0 && stateIndex < dim)
                {
                    ++c;
                    ++freqs;
                }
                else
                {
                    break;
                }

            }

            if (child.isTip())
                charCounts[c.getStateIndex()]++;
            else
                simulate( child, data, rateIndex, charCounts);
        }
    }

}
