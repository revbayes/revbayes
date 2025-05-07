#ifndef PhyloCTMCSiteHomogeneous_H
#define PhyloCTMCSiteHomogeneous_H

#include <cassert>
#include "AbstractPhyloCTMCSiteHomogeneous.h"
#include "DnaState.h"
#include "RateMatrix.h"
#include "RbVector.h"
#include "TopologyNode.h"
#include "TransitionProbabilityMatrix.h"
#include "TypedDistribution.h"

namespace RevBayesCore {

    template<class charType>
    class PhyloCTMCSiteHomogeneous : public AbstractPhyloCTMCSiteHomogeneous<charType> {

    public:
        PhyloCTMCSiteHomogeneous(const TypedDagNode< Tree > *t, size_t nChars, bool c, size_t nSites, bool amb, bool internal, bool gapmatch);
        virtual                                            ~PhyloCTMCSiteHomogeneous(void);                                                                   //!< Virtual destructor

        // public member functions
        PhyloCTMCSiteHomogeneous*                           clone(void) const;                                                                          //!< Create an independent clone


    protected:

        virtual void                                        computeInternalNodeLikelihoodBranchWise(const TopologyNode &n, size_t nIdx, size_t l, size_t r);
        virtual void                                        computeRootLikelihood(size_t root, size_t l, size_t r);
        virtual void                                        computeRootLikelihood(size_t root, size_t l, size_t r, size_t m);
        virtual void                                        computeTipLikelihood(const TopologyNode &node, size_t nIdx);

        virtual void                                        computeInternalNodeLikelihoodNodeWise(const TopologyNode &n, size_t nIdx, size_t l, size_t r);
        virtual void                                        computeRootLikelihoodNode(size_t root, size_t l, size_t r);
        virtual void                                        computeRootLikelihoodNode(size_t root, size_t l, size_t r, size_t m);

        virtual void                                        computeInternalNodeLikelihoodBranchNodeWise(const TopologyNode &n, size_t nIdx, size_t l, size_t r);
        virtual void                                        computeRootLikelihoodBranchNode(size_t root, size_t l, size_t r);
        virtual void                                        computeRootLikelihoodBranchNode(size_t root, size_t l, size_t r, size_t m);


    private:



    };

}


#include "ConstantNode.h"
#include "DiscreteCharacterState.h"
#include "RateMatrix_JC.h"
#include "RandomNumberFactory.h"

#include <cmath>
#include <cstring>

template<class charType>
RevBayesCore::PhyloCTMCSiteHomogeneous<charType>::PhyloCTMCSiteHomogeneous(const TypedDagNode<Tree> *t, size_t nChars, bool c, size_t nSites, bool amb, bool internal, bool gapmatch) : AbstractPhyloCTMCSiteHomogeneous<charType>(  t, nChars, 1, c, nSites, amb, internal, gapmatch )
{

}


template<class charType>
RevBayesCore::PhyloCTMCSiteHomogeneous<charType>::~PhyloCTMCSiteHomogeneous( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!

}


template<class charType>
RevBayesCore::PhyloCTMCSiteHomogeneous<charType>* RevBayesCore::PhyloCTMCSiteHomogeneous<charType>::clone( void ) const
{

    return new PhyloCTMCSiteHomogeneous<charType>( *this );
}


template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneous<charType>::computeInternalNodeLikelihoodBranchNodeWise(const TopologyNode &node, size_t node_index, size_t left, size_t right)
{

    
    // compute the transition probability matrix
    size_t pmat_offset_left  = this->active_pmatrices[left]  * this->active_P_matrix_offset + left  * this->pmat_node_offset;
    size_t pmat_offset_right = this->active_pmatrices[right] * this->active_P_matrix_offset + right * this->pmat_node_offset;

    // get the pointers to the partial likelihoods for this node and the two descendant subtrees
    double*         p_node         = this->partial_node_likelihoods   + this->active_node_likelihood[node_index]   * this->active_node_likelihood_offset   + (node_index-this->num_tips) * this->node_offset;
    double*         p_branch_left  = this->partial_branch_likelihoods + this->active_branch_likelihood[left]       * this->active_branch_likelihood_offset + left                        * this->node_offset;
    double*         p_branch_right = this->partial_branch_likelihoods + this->active_branch_likelihood[right]      * this->active_branch_likelihood_offset + right                       * this->node_offset;
    const double*   p_left         = NULL;
    const double*   p_right        = NULL;

    // check whether the branches are dirty and need recomputing
    bool left_branch_dirty  = this->dirty_branches[left];
    bool right_branch_dirty = this->dirty_branches[right];
    
    
    bool left_is_tip  = left  < this->num_tips;
    bool right_is_tip = right < this->num_tips;
    if ( left_is_tip   )
    {
        p_left   = this->tip_likelihoods + left   * this->tip_offset;
    }
    else
    {
        p_left  = this->partial_node_likelihoods + this->active_node_likelihood[left]       * this->active_node_likelihood_offset + (left-this->num_tips)       * this->node_offset;
    }
    if ( right_is_tip  )
    {
        p_right  = this->tip_likelihoods + right  * this->tip_offset;
    }
    else
    {
        p_right = this->partial_node_likelihoods + this->active_node_likelihood[right]      * this->active_node_likelihood_offset + (right-this->num_tips)      * this->node_offset;
    }
    
    bool left_use_tip_state  = left_is_tip  && this->using_ambiguous_characters == false && this->using_weighted_characters == false;
    bool right_use_tip_state = right_is_tip && this->using_ambiguous_characters == false && this->using_weighted_characters == false;

    const std::vector<bool>&            left_gap_node   = this->gap_matrix [(left_is_tip  ? left  : 0)];
    const std::vector<bool>&            right_gap_node  = this->gap_matrix [(right_is_tip ? right : 0)];
    const std::vector<std::uint64_t>&   left_char_node  = this->char_matrix[(left_is_tip  ? left  : 0)];
    const std::vector<std::uint64_t>&   right_char_node = this->char_matrix[(right_is_tip ? right : 0)];
    
    if ( left_branch_dirty == true && right_branch_dirty == true )
    {
        // iterate over all mixture categories
        for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
        {
            // the transition probability matrix for this mixture category
            const double* tp_begin_left  = this->pmatrices[pmat_offset_left  + mixture].theMatrix;
            const double* tp_begin_right = this->pmatrices[pmat_offset_right + mixture].theMatrix;

            // get the pointers to the likelihood for this mixture category
            size_t offset = mixture*this->mixture_offset;
            double*          p_node_site_mixture          = p_node          + offset;
            double*          p_branch_site_mixture_left   = p_branch_left   + offset;
            double*          p_branch_site_mixture_right  = p_branch_right  + offset;
            const double*    p_node_site_mixture_left     = p_left          + (left_is_tip  ? 0 : offset);
            const double*    p_node_site_mixture_right    = p_right         + (right_is_tip ? 0 : offset);
            // compute the per site probabilities
            for (size_t site = 0; site < this->pattern_block_size ; ++site)
            {

                // get the pointers for this mixture category and this site
                const double*       tp_a_left    = tp_begin_left;
                const double*       tp_a_right   = tp_begin_right;
                // iterate over the possible starting states
                for (size_t c1 = 0; c1 < this->num_states; ++c1)
                {
                    
                    // initialize the probabilities
                    double sum_left  = 0.0;
                    double sum_right = 0.0;
                    if ( left_use_tip_state == true )
                    {
                        if ( left_gap_node[site] == true )
                        {
                            sum_left = 1.0;
                        }
                        else
                        {
                            sum_left = tp_a_left[left_char_node[site]];
                        }
                    }
                    else
                    {
                        // iterate over all possible terminal states
                        for (size_t c2 = 0; c2 < this->num_states; ++c2 )
                        {
                            sum_left  += p_node_site_mixture_left [c2] * tp_a_left [c2];
                        } // end-for over all distination character
                    }
                    (*p_branch_site_mixture_left)  = sum_left;
                    
                    if ( right_use_tip_state == true )
                    {
                        if ( right_gap_node[site] == true )
                        {
                            sum_right = 1.0;
                        }
                        else
                        {
                            sum_right = tp_a_right[right_char_node[site]];
                        }
                    }
                    else
                    {
                        // iterate over all possible terminal states
                        for (size_t c2 = 0; c2 < this->num_states; ++c2 )
                        {
                            sum_right += p_node_site_mixture_right[c2] * tp_a_right[c2];
                        } // end-for over all distination character
                        
                    }
                    (*p_branch_site_mixture_right)  = sum_right;

                    // store the likelihood for this starting state
                    (*p_node_site_mixture) = sum_left * sum_right;

                    assert(isnan(*p_node_site_mixture) || (0 <= *p_node_site_mixture and *p_node_site_mixture <= 1.00000000001) || *p_node_site_mixture <= 0);

                    // increment the pointers to the next starting state
                    tp_a_left  += this->num_states;
                    tp_a_right += this->num_states;
                    
                    ++p_branch_site_mixture_left;
                    ++p_branch_site_mixture_right;
                    ++p_node_site_mixture;

                } // end-for over all initial characters

                // increment the pointers to the next site
                p_node_site_mixture_left    += this->site_offset;
                p_node_site_mixture_right   += this->site_offset;

            } // end-for over all sites (=patterns)
        
        } // end-for over all mixture categories

    }
    else if ( left_branch_dirty == true )
    {
        
        // iterate over all mixture categories
        for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
        {
            // the transition probability matrix for this mixture category
            const double* tp_begin_left  = this->pmatrices[pmat_offset_left  + mixture].theMatrix;

            // get the pointers to the likelihood for this mixture category
            size_t offset = mixture*this->mixture_offset;
            double*          p_node_site_mixture          = p_node          + offset;
            double*          p_branch_site_mixture_left   = p_branch_left   + offset;
            double*          p_branch_site_mixture_right  = p_branch_right  + offset;
            const double*    p_node_site_mixture_left     = p_left          + (left_is_tip  ? 0 : offset);
            // compute the per site probabilities
            for (size_t site = 0; site < this->pattern_block_size ; ++site)
            {

                // get the pointers for this mixture category and this site
                const double*       tp_a_left    = tp_begin_left;
                // iterate over the possible starting states
                for (size_t c1 = 0; c1 < this->num_states; ++c1)
                {
                        
                    double sum_left  = 0.0;
                    if ( left_use_tip_state == true )
                    {
                        if ( left_gap_node[site] == true )
                        {
                            sum_left = 1.0;
                        }
                        else
                        {
                            sum_left = tp_a_left[left_char_node[site]];
                        }
                    }
                    else
                    {
                        // iterate over all possible terminal states
                        for (size_t c2 = 0; c2 < this->num_states; ++c2 )
                        {
                            sum_left  += p_node_site_mixture_left [c2] * tp_a_left [c2];
                        } // end-for over all distination character
                    }
                    (*p_branch_site_mixture_left)  = sum_left;

                    // store the likelihood for this starting state
                    (*p_node_site_mixture) = sum_left * *p_branch_site_mixture_right;

                    assert(isnan(*p_node_site_mixture) || (0 <= *p_node_site_mixture and *p_node_site_mixture <= 1.00000000001) || *p_node_site_mixture <= 0);

                    // increment the pointers to the next starting state
                    tp_a_left  += this->num_states;
                    
                    ++p_branch_site_mixture_left;
                    ++p_branch_site_mixture_right;
                    ++p_node_site_mixture;
                    
                } // end-for over all initial characters

                // increment the pointers to the next site
                p_node_site_mixture_left    += this->site_offset;
                
            } // end-for over all sites (=patterns)
            
        } // end-for over all mixture categories
    
    }
    else if ( right_branch_dirty == true )
    {
        // iterate over all mixture categories
        for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
        {
            // the transition probability matrix for this mixture category
            const double* tp_begin_right = this->pmatrices[pmat_offset_right + mixture].theMatrix;

            // get the pointers to the likelihood for this mixture category
            size_t offset = mixture*this->mixture_offset;
            double*          p_node_site_mixture          = p_node          + offset;
            double*          p_branch_site_mixture_left   = p_branch_left   + offset;
            double*          p_branch_site_mixture_right  = p_branch_right  + offset;
            const double*    p_node_site_mixture_right    = p_right         + (right_is_tip ? 0 : offset);
            // compute the per site probabilities
            for (size_t site = 0; site < this->pattern_block_size ; ++site)
            {

                // get the pointers for this mixture category and this site
                const double*       tp_a_right   = tp_begin_right;
                // iterate over the possible starting states
                for (size_t c1 = 0; c1 < this->num_states; ++c1)
                {
                    // initialize the probabilities
                    double sum_right = 0.0;
                    if ( right_use_tip_state == true )
                    {
                        if ( right_gap_node[site] == true )
                        {
                            sum_right = 1.0;
                        }
                        else
                        {
                            sum_right = tp_a_right[right_char_node[site]];
                        }
                    }
                    else
                    {
                        // iterate over all possible terminal states
                        for (size_t c2 = 0; c2 < this->num_states; ++c2 )
                        {
                            sum_right += p_node_site_mixture_right[c2] * tp_a_right[c2];
                        } // end-for over all distination character
                        
                    }
                    (*p_branch_site_mixture_right)  = sum_right;

                    // store the likelihood for this starting state
                    (*p_node_site_mixture) = *p_branch_site_mixture_left * sum_right;

                    assert(isnan(*p_node_site_mixture) || (0 <= *p_node_site_mixture and *p_node_site_mixture <= 1.00000000001) || *p_node_site_mixture <= 0);

                    // increment the pointers to the next starting state
                    tp_a_right += this->num_states;
                    
                    ++p_branch_site_mixture_left;
                    ++p_branch_site_mixture_right;
                    ++p_node_site_mixture;

                } // end-for over all initial characters

                // increment the pointers to the next site
                p_node_site_mixture_right   += this->site_offset;

            } // end-for over all sites (=patterns)
        
        } // end-for over all mixture categories

    }
    else
    {
        // iterate over all mixture categories
        for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
        {

            // get the pointers to the likelihood for this mixture category
            size_t offset = mixture*this->mixture_offset;
            double*          p_node_site_mixture          = p_node          + offset;
            double*          p_branch_site_mixture_left   = p_branch_left   + offset;
            double*          p_branch_site_mixture_right  = p_branch_right  + offset;
            // compute the per site probabilities
            for (size_t site = 0; site < this->pattern_block_size ; ++site)
            {

                // iterate over the possible starting states
                for (size_t c1 = 0; c1 < this->num_states; ++c1)
                {

                    // store the likelihood for this starting state
                    *p_node_site_mixture = *p_branch_site_mixture_left * *p_branch_site_mixture_right;

                    assert(isnan(*p_node_site_mixture) || (0 <= *p_node_site_mixture and *p_node_site_mixture <= 1.00000000001) || *p_node_site_mixture <= 0);
                    
                    ++p_branch_site_mixture_left;
                    ++p_branch_site_mixture_right;
                    ++p_node_site_mixture;
                    
                } // end-for over all initial characters

            } // end-for over all sites (=patterns)
        
        } // end-for over all mixture categories

    }

    this->dirty_branches[left]  = false;
    this->dirty_branches[right] = false;
}



template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneous<charType>::computeInternalNodeLikelihoodBranchWise(const TopologyNode &node, size_t node_index, size_t left, size_t right)
{

    // compute the transition probability matrix
    size_t pmat_offset = this->active_pmatrices[node_index] * this->active_P_matrix_offset + node_index * this->pmat_node_offset;

    // get the pointers to the partial likelihoods for this node and the two descendant subtrees
    const double*   p_left  = this->partial_branch_likelihoods + this->active_branch_likelihood[left]       * this->active_branch_likelihood_offset + left       * this->node_offset;
    const double*   p_right = this->partial_branch_likelihoods + this->active_branch_likelihood[right]      * this->active_branch_likelihood_offset + right      * this->node_offset;
    double*         p_node  = this->partial_branch_likelihoods + this->active_branch_likelihood[node_index] * this->active_branch_likelihood_offset + node_index * this->node_offset;
    
    bool test_underflow  = RbSettings::userSettings().getUseScaling() == true;
    bool test_this_node  = ((node_index+1) % RbSettings::userSettings().getScalingDensity() == 0);
    bool scale_threshold = RbSettings::userSettings().getScalingMethod() == "threshold";
    bool scale_per_mixture  = RbSettings::userSettings().getScalingPerMixture();

    test_underflow = test_underflow && scale_per_mixture;
    
    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
    {
        // the transition probability matrix for this mixture category
        const double* tp_begin = this->pmatrices[pmat_offset + mixture].theMatrix;

        // get the pointers to the likelihood for this mixture category
        size_t offset = mixture*this->mixture_offset;
        double*          p_site_mixture          = p_node + offset;
        const double*    p_site_mixture_left     = p_left + offset;
        const double*    p_site_mixture_right    = p_right + offset;
        // compute the per site probabilities
        for (size_t site = 0; site < this->pattern_block_size ; ++site)
        {

            // get the pointers for this mixture category and this site
            const double*       tp_a    = tp_begin;
            // iterate over the possible starting states
            for (size_t c1 = 0; c1 < this->num_states; ++c1)
            {
                // temporary variable
                double sum = 0.0;

                // iterate over all possible terminal states
                for (size_t c2 = 0; c2 < this->num_states; ++c2 )
                {
                    sum += p_site_mixture_left[c2] * p_site_mixture_right[c2] * tp_a[c2];
                } // end-for over all distination character

                // store the likelihood for this starting state
                p_site_mixture[c1] = sum;

                assert(isnan(sum) || (0 <= sum and sum <= 1.00000000001));

                // increment the pointers to the next starting state
                tp_a+=this->num_states;

            } // end-for over all initial characters

            if ( test_underflow )
            {
                if ( test_this_node == false )
                {
                    this->per_node_site_mixture_log_scaling_factors[this->active_branch_likelihood[node_index]][node_index][mixture][site] = this->per_node_site_mixture_log_scaling_factors[this->active_branch_likelihood[left]][left][mixture][site] + this->per_node_site_mixture_log_scaling_factors[this->active_branch_likelihood[right]][right][mixture][site];
                }
                else
                {
                    
                    // iterate over all possible terminal states
                    double max = p_site_mixture[0];
                    for (size_t c2 = 1; c2 < this->num_states; ++c2 )
                    {
                        max = ( p_site_mixture[c2] > max ? p_site_mixture[c2] : max );

                    } // end-for over all distination character
                    
                    if ( scale_threshold == false )
                    {
                        if ( max > 0 )
                        {
                            this->per_node_site_mixture_log_scaling_factors[this->active_branch_likelihood[node_index]][node_index][mixture][site] = this->per_node_site_mixture_log_scaling_factors[this->active_branch_likelihood[left]][left][mixture][site] + this->per_node_site_mixture_log_scaling_factors[this->active_branch_likelihood[right]][right][mixture][site] - log(max);
                            for (size_t c2 = 0; c2 < this->num_states; ++c2 )
                            {
                                p_site_mixture[c2] /= max;
                            } // end-for over all distination character
                        }
                        else
                        {
                            this->per_node_site_mixture_log_scaling_factors[this->active_branch_likelihood[node_index]][node_index][mixture][site] = this->per_node_site_mixture_log_scaling_factors[this->active_branch_likelihood[left]][left][mixture][site] + this->per_node_site_mixture_log_scaling_factors[this->active_branch_likelihood[right]][right][mixture][site];
                        }
                    }
                    else if ( max < RbConstants::SCALING_THRESHOLD )
                    {
                        this->per_node_site_mixture_log_scaling_factors[this->active_branch_likelihood[node_index]][node_index][mixture][site] = this->per_node_site_mixture_log_scaling_factors[this->active_branch_likelihood[left]][left][mixture][site] + this->per_node_site_mixture_log_scaling_factors[this->active_branch_likelihood[right]][right][mixture][site] + 1;
                        for (size_t c2 = 0; c2 < this->num_states; ++c2 )
                        {
                            p_site_mixture[c2] /= RbConstants::SCALING_THRESHOLD;

                        } // end-for over all distination character
                    }
                    else
                    {
                        this->per_node_site_mixture_log_scaling_factors[this->active_branch_likelihood[node_index]][node_index][mixture][site] = this->per_node_site_mixture_log_scaling_factors[this->active_branch_likelihood[left]][left][mixture][site] + this->per_node_site_mixture_log_scaling_factors[this->active_branch_likelihood[right]][right][mixture][site];
                    }
                } // end else if we don't test this node
            } // end-if test for underflow
            
            // increment the pointers to the next site
            p_site_mixture_left  += this->site_offset;
            p_site_mixture_right += this->site_offset;
            p_site_mixture       += this->site_offset;

        } // end-for over all sites (=patterns)

    } // end-for over all mixtures (=rate-categories)

}


template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneous<charType>::computeInternalNodeLikelihoodNodeWise(const TopologyNode &node, size_t node_index, size_t left, size_t right)
{

    // compute the transition probability matrix
    size_t pmat_offset_left  = this->active_pmatrices[left]  * this->active_P_matrix_offset + left  * this->pmat_node_offset;
    size_t pmat_offset_right = this->active_pmatrices[right] * this->active_P_matrix_offset + right * this->pmat_node_offset;

    // get the pointers to the partial likelihoods for this node and the two descendant subtrees
    double*         p_node  = this->partial_node_likelihoods + this->active_node_likelihood[node_index] * this->active_node_likelihood_offset + (node_index-this->num_tips) * this->node_offset;
    const double*   p_left  = NULL;
    const double*   p_right = NULL;

    bool left_is_tip  = left  < this->num_tips;
    bool right_is_tip = right < this->num_tips;
    if ( left_is_tip   )
    {
        p_left   = this->tip_likelihoods + left   * this->tip_offset;
    }
    else
    {
        p_left  = this->partial_node_likelihoods + this->active_node_likelihood[left]       * this->active_node_likelihood_offset + (left-this->num_tips)       * this->node_offset;
    }
    if ( right_is_tip  )
    {
        p_right  = this->tip_likelihoods + right  * this->tip_offset;
    }
    else
    {
        p_right = this->partial_node_likelihoods + this->active_node_likelihood[right]      * this->active_node_likelihood_offset + (right-this->num_tips)      * this->node_offset;
    }
    
    bool left_use_tip_state  = left_is_tip  && this->using_ambiguous_characters == false && this->using_weighted_characters == false;
    bool right_use_tip_state = right_is_tip && this->using_ambiguous_characters == false && this->using_weighted_characters == false;

    const std::vector<bool>&            left_gap_node   = this->gap_matrix [(left_is_tip  ? left  : 0)];
    const std::vector<bool>&            right_gap_node  = this->gap_matrix [(right_is_tip ? right : 0)];
    const std::vector<std::uint64_t>&   left_char_node  = this->char_matrix[(left_is_tip  ? left  : 0)];
    const std::vector<std::uint64_t>&   right_char_node = this->char_matrix[(right_is_tip ? right : 0)];
    
    
    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
    {
        // the transition probability matrix for this mixture category
        const double* tp_begin_left  = this->pmatrices[pmat_offset_left  + mixture].theMatrix;
        const double* tp_begin_right = this->pmatrices[pmat_offset_right + mixture].theMatrix;

        // get the pointers to the likelihood for this mixture category
        size_t offset = mixture*this->mixture_offset;
        double*          p_node_site_mixture          = p_node  + offset;
        const double*    p_node_site_mixture_left     = p_left  + (left_is_tip  ? 0 : offset);
        const double*    p_node_site_mixture_right    = p_right + (right_is_tip ? 0 : offset);
        // compute the per site probabilities
        for (size_t site = 0; site < this->pattern_block_size ; ++site)
        {

            // get the pointers for this mixture category and this site
            const double*       tp_a_left    = tp_begin_left;
            const double*       tp_a_right   = tp_begin_right;
            // iterate over the possible starting states
            for (size_t c1 = 0; c1 < this->num_states; ++c1)
            {
                // temporary variable
                double sum_left  = 0.0;
                double sum_right = 0.0;
                
                if ( left_use_tip_state == true )
                {
                    if ( left_gap_node[site] == true )
                    {
                        sum_left = 1.0;
                    }
                    else
                    {
                        sum_left = tp_a_left[left_char_node[site]];
                    }
                }
                else
                {
                    // iterate over all possible terminal states
                    for (size_t c2 = 0; c2 < this->num_states; ++c2 )
                    {
                        sum_left  += p_node_site_mixture_left [c2] * tp_a_left [c2];
                    } // end-for over all distination character
                }
                
                if ( right_use_tip_state == true )
                {
                    if ( right_gap_node[site] == true )
                    {
                        sum_right = 1.0;
                    }
                    else
                    {
                        sum_right = tp_a_right[right_char_node[site]];
                    }
                }
                else
                {
                    // iterate over all possible terminal states
                    for (size_t c2 = 0; c2 < this->num_states; ++c2 )
                    {
                        sum_right  += p_node_site_mixture_right[c2] * tp_a_right[c2];
                    } // end-for over all distination character
                }

                // store the likelihood for this starting state
                p_node_site_mixture[c1] = sum_left * sum_right;

                assert(isnan(p_node_site_mixture[c1]) || (0 <= p_node_site_mixture[c1] and p_node_site_mixture[c1] <= 1.00000000001));

                // increment the pointers to the next starting state
                tp_a_left  += this->num_states;
                tp_a_right += this->num_states;

            } // end-for over all initial characters

            // increment the pointers to the next site
            p_node_site_mixture_left  += this->site_offset;
            p_node_site_mixture_right += this->site_offset;
            p_node_site_mixture       += this->site_offset;

        } // end-for over all sites (=patterns)

    } // end-for over all mixtures (=rate-categories)

}


template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneous<charType>::computeRootLikelihoodBranchNode( size_t root, size_t left, size_t right)
{
    
    // compute the transition probability matrix
    size_t pmat_offset_left   = this->active_pmatrices[left]   * this->active_P_matrix_offset + left   * this->pmat_node_offset;
    size_t pmat_offset_right  = this->active_pmatrices[right]  * this->active_P_matrix_offset + right  * this->pmat_node_offset;

    // get the pointers to the partial likelihoods of the left and right subtree
    double*         p               = this->partial_node_likelihoods   + this->active_node_likelihood[root]         * this->active_node_likelihood_offset   + (root-this->num_tips)       * this->node_offset;
    double*         p_branch_left   = this->partial_branch_likelihoods + this->active_branch_likelihood[left]       * this->active_branch_likelihood_offset + left                        * this->node_offset;
    double*         p_branch_right  = this->partial_branch_likelihoods + this->active_branch_likelihood[right]      * this->active_branch_likelihood_offset + right                       * this->node_offset;
    const double*   p_left          = NULL;
    const double*   p_right         = NULL;
    
    // check whether the branches are dirty and need recomputing
    bool left_branch_dirty   = this->dirty_branches[left];
    bool right_branch_dirty  = this->dirty_branches[right];

    bool left_is_tip   = left   < this->num_tips;
    bool right_is_tip  = right  < this->num_tips;
    if ( left_is_tip   )
    {
        p_left   = this->tip_likelihoods + left   * this->tip_offset;
    }
    else
    {
        p_left   = this->partial_node_likelihoods + this->active_node_likelihood[left]       * this->active_node_likelihood_offset + (left-this->num_tips)       * this->node_offset;
    }
    if ( right_is_tip  )
    {
        p_right  = this->tip_likelihoods + right  * this->tip_offset;
    }
    else
    {
        p_right  = this->partial_node_likelihoods + this->active_node_likelihood[right]      * this->active_node_likelihood_offset + (right-this->num_tips)      * this->node_offset;
    }
    
    bool left_use_tip_state   = left_is_tip   && this->using_ambiguous_characters == false && this->using_weighted_characters == false;
    bool right_use_tip_state  = right_is_tip  && this->using_ambiguous_characters == false && this->using_weighted_characters == false;

    const std::vector<bool>&            left_gap_node    = this->gap_matrix [(left_is_tip   ? left   : 0)];
    const std::vector<bool>&            right_gap_node   = this->gap_matrix [(right_is_tip  ? right  : 0)];
    const std::vector<std::uint64_t>&   left_char_node   = this->char_matrix[(left_is_tip   ? left   : 0)];
    const std::vector<std::uint64_t>&   right_char_node  = this->char_matrix[(right_is_tip  ? right  : 0)];

    // create a vector for the per mixture likelihoods
    // we need this vector to sum over the different mixture likelihoods
    std::vector<double> per_mixture_Likelihoods = std::vector<double>(this->num_patterns,0.0);

    // get the root frequencies
    std::vector<std::vector<double> >   base_frequencies_vectors;
    this->getRootFrequencies(base_frequencies_vectors);

    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
    {
        
        // get the root frequencies
        const std::vector<double> &base_freqs = base_frequencies_vectors[mixture % base_frequencies_vectors.size()];
        assert(base_freqs.size() == this->num_states);
        
        // the transition probability matrix for this mixture category
        const double* tp_begin_left   = this->pmatrices[pmat_offset_left   + mixture].theMatrix;
        const double* tp_begin_right  = this->pmatrices[pmat_offset_right  + mixture].theMatrix;

        // get the pointers to the likelihood for this mixture category
        size_t offset = mixture*this->mixture_offset;
        double*          p_node_site_mixture            = p                 + offset;
        double*          p_branch_site_mixture_left     = p_branch_left     + offset;
        double*          p_branch_site_mixture_right    = p_branch_right    + offset;
        const double*    p_node_site_mixture_left       = p_left            + (left_is_tip   ? 0 : offset);
        const double*    p_node_site_mixture_right      = p_right           + (right_is_tip  ? 0 : offset);
        
        // iterate over all sites
        for (size_t site = 0; site < this->pattern_block_size; ++site)
        {
            
            // get the pointers for this mixture category and this site
            const double*       tp_a_left    = tp_begin_left;
            const double*       tp_a_right   = tp_begin_right;
            
            // iterate over the possible starting states
            for (size_t c1 = 0; c1 < this->num_states; ++c1)
            {
                // temporary variable
                double sum_left   = 0.0;
                double sum_right  = 0.0;
                
                if ( left_branch_dirty == true )
                {
 
                     if ( left_use_tip_state == true )
                     {
                         if ( left_gap_node[site] == true )
                         {
                             sum_left = 1.0;
                         }
                         else
                         {
                             sum_left = tp_a_left[left_char_node[site]];
                         }
                     }
                     else
                     {
                         // iterate over all possible terminal states
                         for (size_t c2 = 0; c2 < this->num_states; ++c2 )
                         {
                             sum_left  += p_node_site_mixture_left [c2] * tp_a_left [c2];
                         } // end-for over all distination character
                     }
                    (*p_branch_site_mixture_left)  = sum_left;
                }
                else
                {
                    sum_left = (*p_branch_site_mixture_left);
                }
                
                if ( right_branch_dirty == true )
                {
                    if ( right_use_tip_state == true )
                    {
                        if ( right_gap_node[site] == true )
                        {
                            sum_right = 1.0;
                        }
                        else
                        {
                            sum_right = tp_a_right[right_char_node[site]];
                       }
                    }
                    else
                    {
                        // iterate over all possible terminal states
                        for (size_t c2 = 0; c2 < this->num_states; ++c2 )
                        {
                            sum_right += p_node_site_mixture_right[c2] * tp_a_right[c2];
                        } // end-for over all distination character
                    }

                    (*p_branch_site_mixture_right)  = sum_right;
                }
                else
                {
                    sum_right = (*p_branch_site_mixture_right);
                }

                // store the likelihood for this starting state
                (*p_node_site_mixture) = sum_left * sum_right * base_freqs[c1];

                assert(isnan(*p_node_site_mixture) || (0 <= *p_node_site_mixture and *p_node_site_mixture <= 1.00000000001));

                // increment the pointers to the next starting state
                tp_a_left   += this->num_states;
                tp_a_right  += this->num_states;

                ++p_branch_site_mixture_left;
                ++p_branch_site_mixture_right;
                ++p_node_site_mixture;

            } // end-for over all initial characters

            // increment the pointers to the next site
            p_node_site_mixture_left   += this->site_offset;
            p_node_site_mixture_right  += this->site_offset;

        } // end-for over all sites (=patterns)

    } // end-for over all mixtures (=rate categories)
    
    
    this->dirty_branches[left]   = false;
    this->dirty_branches[right]  = false;

}


template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneous<charType>::computeRootLikelihoodBranchNode( size_t root, size_t left, size_t right, size_t middle)
{

    // compute the transition probability matrix
    size_t pmat_offset_left   = this->active_pmatrices[left]   * this->active_P_matrix_offset + left   * this->pmat_node_offset;
    size_t pmat_offset_right  = this->active_pmatrices[right]  * this->active_P_matrix_offset + right  * this->pmat_node_offset;
    size_t pmat_offset_middle = this->active_pmatrices[middle] * this->active_P_matrix_offset + middle * this->pmat_node_offset;

    // get the pointers to the partial likelihoods of the left and right subtree
    double*         p               = this->partial_node_likelihoods   + this->active_node_likelihood[root]         * this->active_node_likelihood_offset   + (root-this->num_tips)       * this->node_offset;
    double*         p_branch_left   = this->partial_branch_likelihoods + this->active_branch_likelihood[left]       * this->active_branch_likelihood_offset + left                        * this->node_offset;
    double*         p_branch_right  = this->partial_branch_likelihoods + this->active_branch_likelihood[right]      * this->active_branch_likelihood_offset + right                       * this->node_offset;
    double*         p_branch_middle = this->partial_branch_likelihoods + this->active_branch_likelihood[middle]     * this->active_branch_likelihood_offset + middle                      * this->node_offset;
    const double*   p_left          = NULL;
    const double*   p_right         = NULL;
    const double*   p_middle        = NULL;
    
    // check whether the branches are dirty and need recomputing
    bool left_branch_dirty   = this->dirty_branches[left];
    bool right_branch_dirty  = this->dirty_branches[right];
    bool middle_branch_dirty = this->dirty_branches[middle];

    bool left_is_tip   = left   < this->num_tips;
    bool right_is_tip  = right  < this->num_tips;
    bool middle_is_tip = middle < this->num_tips;
    if ( left_is_tip   )
    {
        p_left   = this->tip_likelihoods + left   * this->tip_offset;
    }
    else
    {
        p_left   = this->partial_node_likelihoods + this->active_node_likelihood[left]       * this->active_node_likelihood_offset + (left-this->num_tips)       * this->node_offset;
    }
    if ( right_is_tip  )
    {
        p_right  = this->tip_likelihoods + right  * this->tip_offset;
    }
    else
    {
        p_right  = this->partial_node_likelihoods + this->active_node_likelihood[right]      * this->active_node_likelihood_offset + (right-this->num_tips)      * this->node_offset;
    }
    if ( middle_is_tip  )
    {
        p_middle  = this->tip_likelihoods + middle  * this->tip_offset;
    }
    else
    {
        p_middle = this->partial_node_likelihoods + this->active_node_likelihood[middle]     * this->active_node_likelihood_offset + (middle-this->num_tips)      * this->node_offset;
    }
    
    bool left_use_tip_state   = left_is_tip   && this->using_ambiguous_characters == false && this->using_weighted_characters == false;
    bool right_use_tip_state  = right_is_tip  && this->using_ambiguous_characters == false && this->using_weighted_characters == false;
    bool middle_use_tip_state = middle_is_tip && this->using_ambiguous_characters == false && this->using_weighted_characters == false;

    const std::vector<bool>&            left_gap_node    = this->gap_matrix [(left_is_tip   ? left   : 0)];
    const std::vector<bool>&            right_gap_node   = this->gap_matrix [(right_is_tip  ? right  : 0)];
    const std::vector<bool>&            middle_gap_node  = this->gap_matrix [(middle_is_tip ? middle : 0)];
    const std::vector<std::uint64_t>&   left_char_node   = this->char_matrix[(left_is_tip   ? left   : 0)];
    const std::vector<std::uint64_t>&   right_char_node  = this->char_matrix[(right_is_tip  ? right  : 0)];
    const std::vector<std::uint64_t>&   middle_char_node = this->char_matrix[(middle_is_tip ? middle : 0)];

    // create a vector for the per mixture likelihoods
    // we need this vector to sum over the different mixture likelihoods
    std::vector<double> per_mixture_Likelihoods = std::vector<double>(this->num_patterns,0.0);

    // get the root frequencies
    std::vector<std::vector<double> >   base_frequencies_vectors;
    this->getRootFrequencies(base_frequencies_vectors);

    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
    {
        
        // get the root frequencies
        const std::vector<double> &base_freqs = base_frequencies_vectors[mixture % base_frequencies_vectors.size()];
        assert(base_freqs.size() == this->num_states);
        
        // the transition probability matrix for this mixture category
        const double* tp_begin_left   = this->pmatrices[pmat_offset_left   + mixture].theMatrix;
        const double* tp_begin_right  = this->pmatrices[pmat_offset_right  + mixture].theMatrix;
        const double* tp_begin_middle = this->pmatrices[pmat_offset_middle + mixture].theMatrix;

        // get the pointers to the likelihood for this mixture category
        size_t offset = mixture*this->mixture_offset;
        double*          p_node_site_mixture            = p                 + offset;
        double*          p_branch_site_mixture_left     = p_branch_left     + offset;
        double*          p_branch_site_mixture_right    = p_branch_right    + offset;
        double*          p_branch_site_mixture_middle   = p_branch_middle   + offset;
        const double*    p_node_site_mixture_left       = p_left            + (left_is_tip   ? 0 : offset);
        const double*    p_node_site_mixture_right      = p_right           + (right_is_tip  ? 0 : offset);
        const double*    p_node_site_mixture_middle     = p_middle          + (middle_is_tip ? 0 : offset);
        
        // iterate over all sites
        for (size_t site = 0; site < this->pattern_block_size; ++site)
        {
            
            // get the pointers for this mixture category and this site
            const double*       tp_a_left    = tp_begin_left;
            const double*       tp_a_right   = tp_begin_right;
            const double*       tp_a_middle  = tp_begin_middle;
            
            // iterate over the possible starting states
            for (size_t c1 = 0; c1 < this->num_states; ++c1)
            {
                // temporary variable
                double sum_left   = 0.0;
                double sum_right  = 0.0;
                double sum_middle = 0.0;
                
                if ( left_branch_dirty == true )
                {
 
                     if ( left_use_tip_state == true )
                     {
                         if ( left_gap_node[site] == true )
                         {
                             sum_left = 1.0;
                         }
                         else
                         {
                             sum_left = tp_a_left[left_char_node[site]];
                         }
                     }
                     else
                     {
                         // iterate over all possible terminal states
                         for (size_t c2 = 0; c2 < this->num_states; ++c2 )
                         {
                             sum_left  += p_node_site_mixture_left [c2] * tp_a_left [c2];
                         } // end-for over all distination character
                     }
                    (*p_branch_site_mixture_left)  = sum_left;
                }
                else
                {
                    sum_left = (*p_branch_site_mixture_left);
                }
                
                if ( right_branch_dirty == true )
                {
                    if ( right_use_tip_state == true )
                    {
                        if ( right_gap_node[site] == true )
                        {
                            sum_right = 1.0;
                        }
                        else
                        {
                            sum_right = tp_a_right[right_char_node[site]];
                       }
                    }
                    else
                    {
                        // iterate over all possible terminal states
                        for (size_t c2 = 0; c2 < this->num_states; ++c2 )
                        {
                            sum_right += p_node_site_mixture_right[c2] * tp_a_right[c2];
                        } // end-for over all distination character
                    }

                    (*p_branch_site_mixture_right)  = sum_right;
                }
                else
                {
                    sum_right = (*p_branch_site_mixture_right);
                }
                
                if ( middle_branch_dirty == true )
                {
                    if ( middle_use_tip_state == true )
                    {
                        if ( middle_gap_node[site] == true )
                        {
                            sum_middle = 1.0;
                        }
                        else
                        {
                            sum_middle = tp_a_middle[middle_char_node[site]];
                       }
                    }
                    else
                    {
                        // iterate over all possible terminal states
                        for (size_t c2 = 0; c2 < this->num_states; ++c2 )
                        {
                            sum_middle += p_node_site_mixture_middle[c2] * tp_a_middle[c2];
                        } // end-for over all distination character
                    }

                    (*p_branch_site_mixture_middle)  = sum_middle;
                }
                else
                {
                    sum_middle = (*p_branch_site_mixture_middle);
                }

                // store the likelihood for this starting state
                (*p_node_site_mixture) = sum_left * sum_right * sum_middle * base_freqs[c1];

                assert(isnan(*p_node_site_mixture) || (0 <= *p_node_site_mixture and *p_node_site_mixture <= 1.00000000001));

                // increment the pointers to the next starting state
                tp_a_left   += this->num_states;
                tp_a_right  += this->num_states;
                tp_a_middle += this->num_states;

                ++p_branch_site_mixture_left;
                ++p_branch_site_mixture_right;
                ++p_branch_site_mixture_middle;
                ++p_node_site_mixture;

            } // end-for over all initial characters

            // increment the pointers to the next site
            p_node_site_mixture_left   += this->site_offset;
            p_node_site_mixture_right  += this->site_offset;
            p_node_site_mixture_middle += this->site_offset;

        } // end-for over all sites (=patterns)

    } // end-for over all mixtures (=rate categories)    
    
    
    this->dirty_branches[left]   = false;
    this->dirty_branches[right]  = false;
    this->dirty_branches[middle] = false;
    
}




template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneous<charType>::computeRootLikelihood( size_t root, size_t left, size_t right)
{

    // get the pointers to the partial likelihoods of the left and right subtree
          double* p        = this->partial_branch_likelihoods + this->active_branch_likelihood[root]  * this->active_branch_likelihood_offset + root  * this->node_offset;
    const double* p_left   = this->partial_branch_likelihoods + this->active_branch_likelihood[left]  * this->active_branch_likelihood_offset + left  * this->node_offset;
    const double* p_right  = this->partial_branch_likelihoods + this->active_branch_likelihood[right] * this->active_branch_likelihood_offset + right * this->node_offset;

    // create a vector for the per mixture likelihoods
    // we need this vector to sum over the different mixture likelihoods
    std::vector<double> per_mixture_Likelihoods = std::vector<double>(this->num_patterns,0.0);

    // get pointers the likelihood for both subtrees
          double*   p_mixture          = p;
    const double*   p_mixture_left     = p_left;
    const double*   p_mixture_right    = p_right;

    // get the root frequencies
    std::vector<std::vector<double> >   ff;
    this->getRootFrequencies(ff);

    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
    {
        // get the root frequencies
        const std::vector<double> &f                    = ff[mixture % ff.size()];
        assert(f.size() == this->num_states);
        std::vector<double>::const_iterator f_end       = f.end();
        std::vector<double>::const_iterator f_begin     = f.begin();

        // get pointers to the likelihood for this mixture category
              double*   p_site_mixture          = p_mixture;
        const double*   p_site_mixture_left     = p_mixture_left;
        const double*   p_site_mixture_right    = p_mixture_right;
        // iterate over all sites
        for (size_t site = 0; site < this->pattern_block_size; ++site)
        {
            // get the pointer to the stationary frequencies
            std::vector<double>::const_iterator f_j             = f_begin;
            // get the pointers to the likelihoods for this site and mixture category
                  double* p_site_j        = p_site_mixture;
            const double* p_site_left_j   = p_site_mixture_left;
            const double* p_site_right_j  = p_site_mixture_right;
            // iterate over all starting states
            for (; f_j != f_end; ++f_j)
            {
                // add the probability of starting from this state
                *p_site_j = *p_site_left_j * *p_site_right_j * *f_j;

                assert(isnan(*p_site_j) || (0.0 <= *p_site_j and *p_site_j <= 1.00000000001));

                // increment pointers
                ++p_site_j; ++p_site_left_j; ++p_site_right_j;
            }

            // increment the pointers to the next site
            p_site_mixture+=this->site_offset; p_site_mixture_left+=this->site_offset; p_site_mixture_right+=this->site_offset;

        } // end-for over all sites (=patterns)

        // increment the pointers to the next mixture category
        p_mixture+=this->mixture_offset; p_mixture_left+=this->mixture_offset; p_mixture_right+=this->mixture_offset;

    } // end-for over all mixtures (=rate categories)


}


template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneous<charType>::computeRootLikelihood( size_t root, size_t left, size_t right, size_t middle)
{
    
    
    bool test_underflow     = RbSettings::userSettings().getUseScaling() == true;
    bool test_this_node     = ((root+1) % RbSettings::userSettings().getScalingDensity() == 0);
    bool scale_threshold    = RbSettings::userSettings().getScalingMethod() == "threshold";
    bool scale_per_mixture  = RbSettings::userSettings().getScalingPerMixture();

    test_underflow = test_underflow && scale_per_mixture;
    
    // get the pointers to the partial likelihoods of the left and right subtree
          double* p        = this->partial_branch_likelihoods + this->active_branch_likelihood[root]   * this->active_branch_likelihood_offset + root   * this->node_offset;
    const double* p_left   = this->partial_branch_likelihoods + this->active_branch_likelihood[left]   * this->active_branch_likelihood_offset + left   * this->node_offset;
    const double* p_right  = this->partial_branch_likelihoods + this->active_branch_likelihood[right]  * this->active_branch_likelihood_offset + right  * this->node_offset;
    const double* p_middle = this->partial_branch_likelihoods + this->active_branch_likelihood[middle] * this->active_branch_likelihood_offset + middle * this->node_offset;

    // get pointers the likelihood for both subtrees
          double*   p_mixture          = p;
    const double*   p_mixture_left     = p_left;
    const double*   p_mixture_right    = p_right;
    const double*   p_mixture_middle   = p_middle;

    // get the root frequencies
    std::vector<std::vector<double> >   base_frequencies_vector;
    this->getRootFrequencies(base_frequencies_vector);

    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
    {
        
        // get the root frequencies
        const std::vector<double> &base_frequencies     = base_frequencies_vector[mixture % base_frequencies_vector.size()];
        assert(base_frequencies.size() == this->num_states);
        std::vector<double>::const_iterator f_end       = base_frequencies.end();
        std::vector<double>::const_iterator f_begin     = base_frequencies.begin();

        // get pointers to the likelihood for this mixture category
              double*   p_site_mixture          = p_mixture;
        const double*   p_site_mixture_left     = p_mixture_left;
        const double*   p_site_mixture_right    = p_mixture_right;
        const double*   p_site_mixture_middle   = p_mixture_middle;
        // iterate over all sites
        for (size_t site = 0; site < this->pattern_block_size; ++site)
        {

            // get the pointer to the stationary frequencies
            std::vector<double>::const_iterator f_j = f_begin;
            // get the pointers to the likelihoods for this site and mixture category
                  double* p_site_j        = p_site_mixture;
            const double* p_site_left_j   = p_site_mixture_left;
            const double* p_site_right_j  = p_site_mixture_right;
            const double* p_site_middle_j = p_site_mixture_middle;
            // iterate over all starting states
            for (; f_j != f_end; ++f_j)
            {
                // add the probability of starting from this state
                *p_site_j = *p_site_left_j * *p_site_right_j * *p_site_middle_j * *f_j;

                assert(isnan(*p_site_j) || (0.0 <= *p_site_j and *p_site_j <= 1.00000000001));

                // increment pointers
                ++p_site_j;
                ++p_site_left_j;
                ++p_site_right_j;
                ++p_site_middle_j;
            }
            
            if ( test_underflow )
            {
                if ( test_this_node == false )
                {
                    this->per_node_site_mixture_log_scaling_factors[this->active_branch_likelihood[root]][root][mixture][site] = this->per_node_site_mixture_log_scaling_factors[this->active_branch_likelihood[left]][left][mixture][site] + this->per_node_site_mixture_log_scaling_factors[this->active_branch_likelihood[right]][right][mixture][site] + this->per_node_site_mixture_log_scaling_factors[this->active_branch_likelihood[middle]][middle][mixture][site];
                }
                else
                {
                    
                    // iterate over all possible terminal states
                    double max = p_site_mixture[0];
                    for (size_t c2 = 1; c2 < this->num_states; ++c2 )
                    {
                        max = ( p_site_mixture[c2] > max ? p_site_mixture[c2] : max );

                    } // end-for over all distination character
                    
                    if ( scale_threshold == false )
                    {
                        if ( max > 0 )
                        {
                            this->per_node_site_mixture_log_scaling_factors[this->active_branch_likelihood[root]][root][mixture][site] = this->per_node_site_mixture_log_scaling_factors[this->active_branch_likelihood[left]][left][mixture][site] + this->per_node_site_mixture_log_scaling_factors[this->active_branch_likelihood[right]][right][mixture][site] + this->per_node_site_mixture_log_scaling_factors[this->active_branch_likelihood[middle]][middle][mixture][site] - log(max);
                            for (size_t c2 = 0; c2 < this->num_states; ++c2 )
                            {
                                p_site_mixture[c2] /= max;
                            } // end-for over all distination character
                        }
                        else
                        {
                            this->per_node_site_mixture_log_scaling_factors[this->active_branch_likelihood[root]][root][mixture][site] = this->per_node_site_mixture_log_scaling_factors[this->active_branch_likelihood[left]][left][mixture][site] + this->per_node_site_mixture_log_scaling_factors[this->active_branch_likelihood[right]][right][mixture][site] + this->per_node_site_mixture_log_scaling_factors[this->active_branch_likelihood[middle]][middle][mixture][site];
                        }
                    }
                    else if ( max < RbConstants::SCALING_THRESHOLD )
                    {
                        this->per_node_site_mixture_log_scaling_factors[this->active_branch_likelihood[root]][root][mixture][site] = this->per_node_site_mixture_log_scaling_factors[this->active_branch_likelihood[left]][left][mixture][site] + this->per_node_site_mixture_log_scaling_factors[this->active_branch_likelihood[right]][right][mixture][site] + this->per_node_site_mixture_log_scaling_factors[this->active_branch_likelihood[middle]][middle][mixture][site] + 1;
                        for (size_t c2 = 0; c2 < this->num_states; ++c2 )
                        {
                            p_site_mixture[c2] /= RbConstants::SCALING_THRESHOLD;

                        } // end-for over all distination character
                    }
                    else
                    {
                        this->per_node_site_mixture_log_scaling_factors[this->active_branch_likelihood[root]][root][mixture][site] = this->per_node_site_mixture_log_scaling_factors[this->active_branch_likelihood[left]][left][mixture][site] + this->per_node_site_mixture_log_scaling_factors[this->active_branch_likelihood[right]][right][mixture][site] + this->per_node_site_mixture_log_scaling_factors[this->active_branch_likelihood[middle]][middle][mixture][site];
                    }
                }
            }

            // increment the pointers to the next site
            p_site_mixture+=this->site_offset; p_site_mixture_left+=this->site_offset; p_site_mixture_right+=this->site_offset; p_site_mixture_middle+=this->site_offset;

        } // end-for over all sites (=patterns)

        // increment the pointers to the next mixture category
        p_mixture+=this->mixture_offset; p_mixture_left+=this->mixture_offset; p_mixture_right+=this->mixture_offset; p_mixture_middle+=this->mixture_offset;

    } // end-for over all mixtures (=rate categories)

}



template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneous<charType>::computeRootLikelihoodNode( size_t root, size_t left, size_t right)
{
    
    // compute the transition probability matrix
    size_t pmat_offset_left   = this->active_pmatrices[left]   * this->active_P_matrix_offset + left   * this->pmat_node_offset;
    size_t pmat_offset_right  = this->active_pmatrices[right]  * this->active_P_matrix_offset + right  * this->pmat_node_offset;

    // get the pointers to the partial likelihoods of the left and right subtree
          double*   p        = this->partial_node_likelihoods + this->active_node_likelihood[root]   * this->active_node_likelihood_offset + (root-this->num_tips)   * this->node_offset;
    const double*   p_left   = NULL;
    const double*   p_right  = NULL;
    
    bool left_is_tip   = left   < this->num_tips;
    bool right_is_tip  = right  < this->num_tips;
    if ( left_is_tip   )
    {
        p_left   = this->tip_likelihoods + left   * this->tip_offset;
    }
    else
    {
        p_left   = this->partial_node_likelihoods + this->active_node_likelihood[left]       * this->active_node_likelihood_offset + (left-this->num_tips)       * this->node_offset;
    }
    if ( right_is_tip  )
    {
        p_right  = this->tip_likelihoods + right  * this->tip_offset;
    }
    else
    {
        p_right  = this->partial_node_likelihoods + this->active_node_likelihood[right]      * this->active_node_likelihood_offset + (right-this->num_tips)      * this->node_offset;
    }

    // create a vector for the per mixture likelihoods
    // we need this vector to sum over the different mixture likelihoods
    std::vector<double> per_mixture_Likelihoods = std::vector<double>(this->num_patterns,0.0);

    // get the root frequencies
    std::vector<std::vector<double> >   base_frequencies_vectors;
    this->getRootFrequencies(base_frequencies_vectors);

    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
    {
        
        // get the root frequencies
        const std::vector<double> &base_freqs = base_frequencies_vectors[mixture % base_frequencies_vectors.size()];
        assert(base_freqs.size() == this->num_states);
        
        // the transition probability matrix for this mixture category
        const double* tp_begin_left   = this->pmatrices[pmat_offset_left   + mixture].theMatrix;
        const double* tp_begin_right  = this->pmatrices[pmat_offset_right  + mixture].theMatrix;

        // get the pointers to the likelihood for this mixture category
        size_t offset = mixture*this->mixture_offset;
        double*          p_site_mixture          = p        + offset;
        const double*    p_site_mixture_left     = p_left   + (left_is_tip   ? 0 : offset);
        const double*    p_site_mixture_right    = p_right  + (right_is_tip  ? 0 : offset);
        
        // iterate over all sites
        for (size_t site = 0; site < this->pattern_block_size; ++site)
        {
            
            // get the pointers for this mixture category and this site
            const double*       tp_a_left    = tp_begin_left;
            const double*       tp_a_right   = tp_begin_right;
            
            // iterate over the possible starting states
            for (size_t c1 = 0; c1 < this->num_states; ++c1)
            {
                // temporary variable
                double sum_left   = 0.0;
                double sum_right  = 0.0;

                // iterate over all possible terminal states
                for (size_t c2 = 0; c2 < this->num_states; ++c2 )
                {
                    sum_left   += p_site_mixture_left  [c2] * tp_a_left  [c2];
                    sum_right  += p_site_mixture_right [c2] * tp_a_right [c2];
                } // end-for over all distination character

                // store the likelihood for this starting state
                p_site_mixture[c1] = sum_left * sum_right * base_freqs[c1];

                assert(isnan(p_site_mixture[c1]) || (0 <= p_site_mixture[c1] and p_site_mixture[c1] <= 1.00000000001));

                // increment the pointers to the next starting state
                tp_a_left   += this->num_states;
                tp_a_right  += this->num_states;

            } // end-for over all initial characters

            // increment the pointers to the next site
            p_site_mixture_left   += this->site_offset;
            p_site_mixture_right  += this->site_offset;
            p_site_mixture        += this->site_offset;

        } // end-for over all sites (=patterns)

    } // end-for over all mixtures (=rate categories)

}


template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneous<charType>::computeRootLikelihoodNode( size_t root, size_t left, size_t right, size_t middle)
{

    // compute the transition probability matrix
    size_t pmat_offset_left   = this->active_pmatrices[left]   * this->active_P_matrix_offset + left   * this->pmat_node_offset;
    size_t pmat_offset_right  = this->active_pmatrices[right]  * this->active_P_matrix_offset + right  * this->pmat_node_offset;
    size_t pmat_offset_middle = this->active_pmatrices[middle] * this->active_P_matrix_offset + middle * this->pmat_node_offset;

    // get the pointers to the partial likelihoods of the left and right subtree
          double*   p        = this->partial_node_likelihoods + this->active_node_likelihood[root]   * this->active_node_likelihood_offset + (root-this->num_tips)   * this->node_offset;
    const double*   p_left   = NULL;
    const double*   p_right  = NULL;
    const double*   p_middle = NULL;
    
    bool left_is_tip   = left   < this->num_tips;
    bool right_is_tip  = right  < this->num_tips;
    bool middle_is_tip = middle < this->num_tips;
    if ( left_is_tip   )
    {
        p_left   = this->tip_likelihoods + left   * this->tip_offset;
    }
    else
    {
        p_left   = this->partial_node_likelihoods + this->active_node_likelihood[left]       * this->active_node_likelihood_offset + (left-this->num_tips)       * this->node_offset;
    }
    if ( right_is_tip  )
    {
        p_right  = this->tip_likelihoods + right  * this->tip_offset;
    }
    else
    {
        p_right  = this->partial_node_likelihoods + this->active_node_likelihood[right]      * this->active_node_likelihood_offset + (right-this->num_tips)      * this->node_offset;
    }
    if ( middle_is_tip  )
    {
        p_middle  = this->tip_likelihoods + middle  * this->tip_offset;
    }
    else
    {
        p_middle = this->partial_node_likelihoods + this->active_node_likelihood[middle]     * this->active_node_likelihood_offset + (middle-this->num_tips)      * this->node_offset;
    }

    // create a vector for the per mixture likelihoods
    // we need this vector to sum over the different mixture likelihoods
    std::vector<double> per_mixture_Likelihoods = std::vector<double>(this->num_patterns,0.0);

    // get the root frequencies
    std::vector<std::vector<double> >   base_frequencies_vectors;
    this->getRootFrequencies(base_frequencies_vectors);

    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
    {
        
        // get the root frequencies
        const std::vector<double> &base_freqs = base_frequencies_vectors[mixture % base_frequencies_vectors.size()];
        assert(base_freqs.size() == this->num_states);
        
        // the transition probability matrix for this mixture category
        const double* tp_begin_left   = this->pmatrices[pmat_offset_left   + mixture].theMatrix;
        const double* tp_begin_right  = this->pmatrices[pmat_offset_right  + mixture].theMatrix;
        const double* tp_begin_middle = this->pmatrices[pmat_offset_middle + mixture].theMatrix;

        // get the pointers to the likelihood for this mixture category
        size_t offset = mixture*this->mixture_offset;
        double*          p_site_mixture          = p        + offset;
        const double*    p_site_mixture_left     = p_left   + (left_is_tip   ? 0 : offset);
        const double*    p_site_mixture_right    = p_right  + (right_is_tip  ? 0 : offset);
        const double*    p_site_mixture_middle   = p_middle + (middle_is_tip ? 0 : offset);
        
        // iterate over all sites
        for (size_t site = 0; site < this->pattern_block_size; ++site)
        {
            
            // get the pointers for this mixture category and this site
            const double*       tp_a_left    = tp_begin_left;
            const double*       tp_a_right   = tp_begin_right;
            const double*       tp_a_middle  = tp_begin_middle;
            
            // iterate over the possible starting states
            for (size_t c1 = 0; c1 < this->num_states; ++c1)
            {
                // temporary variable
                double sum_left   = 0.0;
                double sum_right  = 0.0;
                double sum_middle = 0.0;

                // iterate over all possible terminal states
                for (size_t c2 = 0; c2 < this->num_states; ++c2 )
                {
                    sum_left   += p_site_mixture_left  [c2] * tp_a_left  [c2];
                    sum_right  += p_site_mixture_right [c2] * tp_a_right [c2];
                    sum_middle += p_site_mixture_middle[c2] * tp_a_middle[c2];
                } // end-for over all distination character

                // store the likelihood for this starting state
                p_site_mixture[c1] = sum_left * sum_right * sum_middle * base_freqs[c1];

                assert(isnan(p_site_mixture[c1]) || (0 <= p_site_mixture[c1] and p_site_mixture[c1] <= 1.00000000001));

                // increment the pointers to the next starting state
                tp_a_left   += this->num_states;
                tp_a_right  += this->num_states;
                tp_a_middle += this->num_states;

            } // end-for over all initial characters

            // increment the pointers to the next site
            p_site_mixture_left   += this->site_offset;
            p_site_mixture_right  += this->site_offset;
            p_site_mixture_middle += this->site_offset;
            p_site_mixture        += this->site_offset;

        } // end-for over all sites (=patterns)

    } // end-for over all mixtures (=rate categories)

}



template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneous<charType>::computeTipLikelihood(const TopologyNode &node, size_t node_index)
{

    double* p_node = this->partial_branch_likelihoods + this->active_branch_likelihood[node_index]*this->active_branch_likelihood_offset + node_index*this->node_offset;
    
    // get the current correct tip index in case the whole tree change (after performing an empiricalTree Proposal)
    size_t data_tip_index = this->taxon_name_2_tip_index_map[ node.getName() ];
    const std::vector<bool>&            gap_node        = this->gap_matrix[data_tip_index];
    const std::vector<std::uint64_t>&   char_node       = this->char_matrix[data_tip_index];
    const std::vector<RbBitSet>&        amb_char_node   = this->ambiguous_char_matrix[data_tip_index];

    size_t char_data_node_index = this->value->indexOfTaxonWithName(node.getName());
    std::vector<size_t> site_indices;
    if ( this->using_weighted_characters == true )
        site_indices = this->getIncludedSiteIndices();
    
    // compute the transition probabilities
    size_t pmat_offset = this->active_pmatrices[node_index] * this->active_P_matrix_offset + node_index * this->pmat_node_offset;

    double* p_mixture = p_node;

    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
    {
        // the transition probability matrix for this mixture category
        const double* tp_begin = this->pmatrices[pmat_offset + mixture].theMatrix;

        // get the pointer to the likelihoods for this site and mixture category
        double* p_site_mixture = p_mixture;

        // iterate over all sites
        for (size_t site = 0; site != this->pattern_block_size; ++site)
        {

            // is this site a gap?
            if ( gap_node[site] )
            {
                // since this is a gap we need to assume that the actual state could have been any state

                // iterate over all initial states for the transitions
                for (size_t c1 = 0; c1 < this->num_states; ++c1)
                {

                    // store the likelihood
                    p_site_mixture[c1] = 1.0;

                }
            }
            else // we have observed a character
            {

                // iterate over all possible initial states
                for (size_t c1 = 0; c1 < this->num_states; ++c1)
                {

                    if ( this->using_ambiguous_characters == true && this->using_weighted_characters == false)
                    {
                        // compute the likelihood that we had a transition from state c1 to the observed state org_val
                        // note, the observed state could be ambiguous!
                        const RbBitSet &val = amb_char_node[site];

                        // get the pointer to the transition probabilities for the terminal states
                        const double* d  = tp_begin+(this->num_states*c1);

                        double tmp = 0.0;

                        for ( size_t i=0; i<this->num_states; ++i )
                        {
                            // check whether we observed this state
                            if ( val.test(i) == true )
                            {
                                // add the probability
                                tmp += *d;
                            }

                            // increment the pointer to the next transition probability
                            ++d;
                        } // end-while over all observed states for this character

                        // store the likelihood
                        p_site_mixture[c1] = tmp;
                        
                    }
                    else if ( this->using_weighted_characters == true )
                    {
                        // compute the likelihood that we had a transition from state c1 to the observed state org_val
                        // note, the observed state could be ambiguous!
//                        const RbBitSet &val = amb_char_node[site];
                        size_t this_site_index = site_indices[site];
                        const RbBitSet &val = this->value->getCharacter(char_data_node_index, this_site_index).getState();

                        // get the pointer to the transition probabilities for the terminal states
                        const double* d = tp_begin+(this->num_states*c1);

                        double tmp = 0.0;
                        const std::vector< double >& weights = this->value->getCharacter(char_data_node_index, this_site_index).getWeights();
                        for ( size_t i=0; i<this->num_states; ++i )
                        {
                            // check whether we observed this state
                            if ( val.test(i) == true )
                            {
                                // add the probability
                                tmp += *d * weights[i] ;
                            }

                            // increment the pointer to the next transition probability
                            ++d;
                        } // end-while over all observed states for this character

                        // store the likelihood
                        p_site_mixture[c1] = tmp;
                        
                    }
                    else // no ambiguous characters in use
                    {
                        std::uint64_t org_val = char_node[site];

                        // store the likelihood
                        p_site_mixture[c1] = tp_begin[c1*this->num_states+org_val];

                    }

                } // end-for over all possible initial character for the branch

            } // end-if a gap state

            // increment the pointers to next site
            p_site_mixture+=this->site_offset;

        } // end-for over all sites/patterns in the sequence

        // increment the pointers to next mixture category
        p_mixture+=this->mixture_offset;

    } // end-for over all mixture categories

}


#endif
