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

        virtual void                                        computeRootLikelihood(size_t root, size_t l, size_t r);
        virtual void                                        computeRootLikelihood(size_t root, size_t l, size_t r, size_t m);
        virtual void                                        computeInternalNodeLikelihood(const TopologyNode &n, size_t nIdx, size_t l, size_t r);
        virtual void                                        computeInternalNodeLikelihood(const TopologyNode &n, size_t nIdx, size_t l, size_t r, size_t m);
        virtual void                                        computeTipLikelihood(const TopologyNode &node, size_t nIdx);


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
void RevBayesCore::PhyloCTMCSiteHomogeneous<charType>::computeRootLikelihood( size_t root, size_t left, size_t right)
{

    // get the pointers to the partial likelihoods of the left and right subtree
          double* p        = this->getCreatePartialLikelihoodsForNode(root);
    const double* p_left   = this->getPartialLikelihoodsForNode(left);
    const double* p_right  = this->getPartialLikelihoodsForNode(right);

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
        assert(f.size() == this->num_chars);
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
            p_site_mixture+=this->siteOffset; p_site_mixture_left+=this->siteOffset; p_site_mixture_right+=this->siteOffset;

        } // end-for over all sites (=patterns)

        // increment the pointers to the next mixture category
        p_mixture+=this->mixtureOffset; p_mixture_left+=this->mixtureOffset; p_mixture_right+=this->mixtureOffset;

    } // end-for over all mixtures (=rate categories)


}


template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneous<charType>::computeRootLikelihood( size_t root, size_t left, size_t right, size_t middle)
{

    // get the pointers to the partial likelihoods of the left and right subtree
          double* p        = this->getCreatePartialLikelihoodsForNode(root);
    const double* p_left   = this->getPartialLikelihoodsForNode(left);
    const double* p_right  = this->getPartialLikelihoodsForNode(right);
    const double* p_middle = this->getPartialLikelihoodsForNode(middle);

    // get pointers the likelihood for both subtrees
          double*   p_mixture          = p;
    const double*   p_mixture_left     = p_left;
    const double*   p_mixture_right    = p_right;
    const double*   p_mixture_middle   = p_middle;

    // get the root frequencies
    std::vector<std::vector<double> >   ff;
    this->getRootFrequencies(ff);

    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
    {
        
        // get the root frequencies
        const std::vector<double> &f                    = ff[mixture % ff.size()];
        assert(f.size() == this->num_chars);
        std::vector<double>::const_iterator f_end       = f.end();
        std::vector<double>::const_iterator f_begin     = f.begin();

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
                ++p_site_j; ++p_site_left_j; ++p_site_right_j; ++p_site_middle_j;
            }

            // increment the pointers to the next site
            p_site_mixture+=this->siteOffset; p_site_mixture_left+=this->siteOffset; p_site_mixture_right+=this->siteOffset; p_site_mixture_middle+=this->siteOffset;

        } // end-for over all sites (=patterns)

        // increment the pointers to the next mixture category
        p_mixture+=this->mixtureOffset; p_mixture_left+=this->mixtureOffset; p_mixture_right+=this->mixtureOffset; p_mixture_middle+=this->mixtureOffset;

    } // end-for over all mixtures (=rate categories)

}


template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneous<charType>::computeInternalNodeLikelihood(const TopologyNode &node, size_t node_index, size_t left, size_t right)
{

    // compute the transition probability matrix
//    this->updateTransitionProbabilities( node_index );
    size_t pmat_offset = this->active_pmatrices[node_index] * this->activePmatrixOffset + node_index * this->pmatNodeOffset;

    // get the pointers to the partial likelihoods for this node and the two descendant subtrees
    const double*   p_left  = this->getPartialLikelihoodsForNode(left);
    const double*   p_right = this->getPartialLikelihoodsForNode(right);
    double*         p_node  = this->getCreatePartialLikelihoodsForNode(node_index);

    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
    {
        // the transition probability matrix for this mixture category
//        const double*    tp_begin                = this->transition_prob_matrices[mixture].theMatrix;
        const double* tp_begin = this->pmatrices[pmat_offset + mixture].theMatrix;

        // get the pointers to the likelihood for this mixture category
        size_t offset = mixture*this->mixtureOffset;
        double*          p_site_mixture          = p_node + offset;
        const double*    p_site_mixture_left     = p_left + offset;
        const double*    p_site_mixture_right    = p_right + offset;
        // compute the per site probabilities
        for (size_t site = 0; site < this->pattern_block_size ; ++site)
        {

            // get the pointers for this mixture category and this site
            const double*       tp_a    = tp_begin;
            // iterate over the possible starting states
            for (size_t c1 = 0; c1 < this->num_chars; ++c1)
            {
                // temporary variable
                double sum = 0.0;

                // iterate over all possible terminal states
                for (size_t c2 = 0; c2 < this->num_chars; ++c2 )
                {
                    sum += p_site_mixture_left[c2] * p_site_mixture_right[c2] * tp_a[c2];

                } // end-for over all distination character

                // store the likelihood for this starting state
                p_site_mixture[c1] = sum;

                assert(isnan(sum) || (0 <= sum and sum <= 1.00000000001));

                // increment the pointers to the next starting state
                tp_a+=this->num_chars;

            } // end-for over all initial characters

            // increment the pointers to the next site
            p_site_mixture_left+=this->siteOffset; p_site_mixture_right+=this->siteOffset; p_site_mixture+=this->siteOffset;

        } // end-for over all sites (=patterns)

    } // end-for over all mixtures (=rate-categories)

}


template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneous<charType>::computeInternalNodeLikelihood(const TopologyNode &node, size_t node_index, size_t left, size_t right, size_t middle)
{

    // compute the transition probability matrix
//    this->updateTransitionProbabilities( node_index );
    size_t pmat_offset = this->active_pmatrices[node_index] * this->activePmatrixOffset + node_index * this->pmatNodeOffset;

    // get the pointers to the partial likelihoods for this node and the two descendant subtrees
    const double*   p_left      = this->getPartialLikelihoodsForNode(left);
    const double*   p_middle    = this->getPartialLikelihoodsForNode(middle);
    const double*   p_right     = this->getPartialLikelihoodsForNode(right);
    double*         p_node      = this->getCreatePartialLikelihoodsForNode(node_index);

    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
    {
        // the transition probability matrix for this mixture category
//        const double*    tp_begin                = this->transition_prob_matrices[mixture].theMatrix;
        const double* tp_begin = this->pmatrices[pmat_offset + mixture].theMatrix;

        // get the pointers to the likelihood for this mixture category
        size_t offset = mixture*this->mixtureOffset;
        double*          p_site_mixture          = p_node + offset;
        const double*    p_site_mixture_left     = p_left + offset;
        const double*    p_site_mixture_middle   = p_middle + offset;
        const double*    p_site_mixture_right    = p_right + offset;
        // compute the per site probabilities
        for (size_t site = 0; site < this->pattern_block_size ; ++site)
        {

            // get the pointers for this mixture category and this site
            const double*       tp_a    = tp_begin;
            // iterate over the possible starting states
            for (size_t c1 = 0; c1 < this->num_chars; ++c1)
            {
                // temporary variable
                double sum = 0.0;

                // iterate over all possible terminal states
                for (size_t c2 = 0; c2 < this->num_chars; ++c2 )
                {
                    sum += p_site_mixture_left[c2] * p_site_mixture_middle[c2] * p_site_mixture_right[c2] * tp_a[c2];

                } // end-for over all distination character

                assert(isnan(sum) || (0 <= sum and sum <= 1.00000000001));

                // store the likelihood for this starting state
                p_site_mixture[c1] = sum;

                // increment the pointers to the next starting state
                tp_a+=this->num_chars;

            } // end-for over all initial characters

            // increment the pointers to the next site
            p_site_mixture_left+=this->siteOffset; p_site_mixture_middle+=this->siteOffset; p_site_mixture_right+=this->siteOffset; p_site_mixture+=this->siteOffset;

        } // end-for over all sites (=patterns)

    } // end-for over all mixtures (=rate-categories)

}




template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneous<charType>::computeTipLikelihood(const TopologyNode &node, size_t node_index)
{
    double* p_node = this->getCreatePartialLikelihoodsForNode(node_index);
    
    // get the current correct tip index in case the whole tree change (after performing an empiricalTree Proposal)
    size_t data_tip_index = this->taxon_name_2_tip_index_map[ node.getName() ];
    const std::vector<bool> &gap_node = this->gap_matrix[data_tip_index];
    const std::vector<std::uint64_t> &char_node = this->char_matrix[data_tip_index];
    const std::vector<RbBitSet> &amb_char_node = this->ambiguous_char_matrix[data_tip_index];

    size_t char_data_node_index = this->value->indexOfTaxonWithName(node.getName());
    std::vector<size_t> site_indices;
    if ( this->using_weighted_characters == true )
        site_indices = this->getIncludedSiteIndices();
    
    // compute the transition probabilities
//    this->updateTransitionProbabilities( node_index );
    size_t pmat_offset = this->active_pmatrices[node_index] * this->activePmatrixOffset + node_index * this->pmatNodeOffset;

    double* p_mixture = p_node;

    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
    {
        // the transition probability matrix for this mixture category
//         const double* tp_begin = this->transition_prob_matrices[mixture].theMatrix;
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
                for (size_t c1 = 0; c1 < this->num_chars; ++c1)
                {

                    // store the likelihood
                    p_site_mixture[c1] = 1.0;

                }
            }
            else // we have observed a character
            {

                // iterate over all possible initial states
                for (size_t c1 = 0; c1 < this->num_chars; ++c1)
                {

                    if ( this->using_ambiguous_characters == true && this->using_weighted_characters == false)
                    {
                        // compute the likelihood that we had a transition from state c1 to the observed state org_val
                        // note, the observed state could be ambiguous!
                        const RbBitSet &val = amb_char_node[site];

                        // get the pointer to the transition probabilities for the terminal states
                        const double* d  = tp_begin+(this->num_chars*c1);

                        double tmp = 0.0;

                        for ( size_t i=0; i<this->num_chars; ++i )
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
                        const double* d = tp_begin+(this->num_chars*c1);

                        double tmp = 0.0;
                        const std::vector< double >& weights = this->value->getCharacter(char_data_node_index, this_site_index).getWeights();
                        for ( size_t i=0; i<this->num_chars; ++i )
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
                        p_site_mixture[c1] = tp_begin[c1*this->num_chars+org_val];

                    }

                } // end-for over all possible initial character for the branch

            } // end-if a gap state

            // increment the pointers to next site
            p_site_mixture+=this->siteOffset;

        } // end-for over all sites/patterns in the sequence

        // increment the pointers to next mixture category
        p_mixture+=this->mixtureOffset;

    } // end-for over all mixture categories

}


#endif
