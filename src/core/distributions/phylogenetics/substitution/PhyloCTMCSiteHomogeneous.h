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
    auto& pl_left  = this->getPartialLikelihoodsForNode(left);
    auto& pl_right = this->getPartialLikelihoodsForNode(right);
    const double* p_left  = pl_left.likelihoods.data();
    const double* p_right = pl_right.likelihoods.data();
    assert(pl_left.dims() == pl_right.dims());

    auto& pl_root = this->createEmptyPartialLikelihoodsForNode(root, pl_left.dims());
    double* p = pl_root.likelihoods.data();

    // create a vector for the per mixture likelihoods
    // we need this vector to sum over the different mixture likelihoods
    std::vector<double> per_mixture_Likelihoods(this->num_patterns, 0.0);

    // get the root frequencies
    std::vector<std::vector<double>> ff;
    this->getRootFrequencies(ff);

    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
    {
        const std::vector<double>& f = ff[mixture % ff.size()];
        assert(f.size() == this->num_chars);

        for (size_t site = 0; site < this->pattern_block_size; ++site)
        {
            const size_t offset = mixture * this->mixtureOffset + site * this->siteOffset;

            for (size_t j = 0; j < this->num_chars; ++j)
            {
                const double v = p_left[offset + j] * p_right[offset + j] * f[j];
                assert(std::isnan(v) || (0.0 <= v && v <= 1.00000000001));
                p[offset + j] = v;
            }
        }
    }

    this->scale(root, left, right);
}


template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneous<charType>::computeRootLikelihood( size_t root, size_t left, size_t right, size_t middle)
{
    // get the pointers to the partial likelihoods of the left, right, and middle subtrees
    auto& pl_left   = this->getPartialLikelihoodsForNode(left);
    auto& pl_right  = this->getPartialLikelihoodsForNode(right);
    auto& pl_middle = this->getPartialLikelihoodsForNode(middle);
    const double* p_left   = pl_left.likelihoods.data();
    const double* p_right  = pl_right.likelihoods.data();
    const double* p_middle = pl_middle.likelihoods.data();
    assert(pl_left.dims() == pl_right.dims());
    assert(pl_left.dims() == pl_middle.dims());

    auto& pl_node = this->createEmptyPartialLikelihoodsForNode(root, pl_left.dims());
    double* p = pl_node.likelihoods.data();

    // get the root frequencies
    std::vector<std::vector<double>> ff;
    this->getRootFrequencies(ff);

    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
    {
        const std::vector<double>& f = ff[mixture % ff.size()];
        assert(f.size() == this->num_chars);

        for (size_t site = 0; site < this->pattern_block_size; ++site)
        {
            const size_t offset = mixture * this->mixtureOffset + site * this->siteOffset;

            for (size_t j = 0; j < this->num_chars; ++j)
            {
                const double v = p_left[offset + j] * p_right[offset + j] * p_middle[offset + j] * f[j];
                assert(std::isnan(v) || (0.0 <= v && v <= 1.00000000001));
                p[offset + j] = v;
            }
        }
    }

    this->scale(root, left, right, middle);
}


template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneous<charType>::computeInternalNodeLikelihood(const TopologyNode &node, size_t node_index, size_t left, size_t right)
{
    // update transition probability matrices
    this->updateTransitionProbabilityMatrix(node_index);

    // get the pointers to the partial likelihoods for this node and the two descendant subtrees
    auto& pl_left  = this->getPartialLikelihoodsForNode(left);
    auto& pl_right = this->getPartialLikelihoodsForNode(right);
    const double* p_left  = pl_left.likelihoods.data();
    const double* p_right = pl_right.likelihoods.data();
    assert(pl_left.dims() == pl_right.dims());

    auto& pl_node = this->createEmptyPartialLikelihoodsForNode(node_index, pl_left.dims());
    double* p_node = pl_node.likelihoods.data();

    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
    {
        // the transition probability matrix for this mixture category
        const double* tp = this->pmatrices[node_index][mixture].theMatrix;

        // compute the per site probabilities
        for (size_t site = 0; site < this->pattern_block_size; ++site)
        {
            const size_t offset = mixture * this->mixtureOffset + site * this->siteOffset;

            // iterate over the possible starting states
            for (size_t c1 = 0; c1 < this->num_chars; ++c1)
            {
                // iterate over all possible terminal states
                double sum = 0.0;
                for (size_t c2 = 0; c2 < this->num_chars; ++c2)
                {
                    // NOTE that the first part here doesn't depend on c1 and could be hoisted out.
                    sum += p_left[offset + c2] * p_right[offset + c2] * tp[c1 * this->num_chars + c2];
                }

                p_node[offset + c1] = sum;
                assert(std::isnan(sum) || (0.0 <= sum && sum <= 1.00000000001));
            }
        }
    }

    this->scale(node_index, left, right);
}


template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneous<charType>::computeTipLikelihood(const TopologyNode &node, size_t node_index)
{
    // update transition probability matrices
    this->updateTransitionProbabilityMatrix(node_index);

    double* p_node = this->createEmptyPartialLikelihoodsForNode(node_index, {this->num_site_mixtures, this->pattern_block_size, this->num_chars}).likelihoods.data();
    
    // get the current correct tip index in case the whole tree change (after performing an empiricalTree Proposal)
    size_t data_tip_index = this->taxon_name_2_tip_index_map[ node.getName() ];
    const std::vector<bool> &gap_node = this->gap_matrix[data_tip_index];
    const std::vector<std::uint64_t> &char_node = this->char_matrix[data_tip_index];
    const std::vector<RbBitSet>& amb_char_node = this->ambiguous_char_matrix[data_tip_index];

    size_t char_data_node_index = this->value->indexOfTaxonWithName(node.getName());
    std::vector<size_t> site_indices;
    if (this->using_weighted_characters == true)
        site_indices = this->getIncludedSiteIndices();

    double obs_error_prob = 0.0;
    // Basanta: Initialize with a dummy simplex; overwritten by the model if enabled.
    Simplex obs_error_freqs = Simplex(this->num_chars);
    if (this->using_observation_error)
    {
        obs_error_prob  = this->observation_error_probability->getValue();
        obs_error_freqs = this->observation_error_frequencies->getValue();
    }

    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
    {
        // the transition probability matrix for this mixture category
        const double* tp = this->pmatrices[node_index][mixture].theMatrix;

        // iterate over all sites
        for (size_t site = 0; site < this->pattern_block_size; ++site)
        {
            const size_t offset = mixture * this->mixtureOffset + site * this->siteOffset;
            double* p_site = p_node + offset;

            // is this site a gap?
            if (gap_node[site])
            {
                // since this is a gap we need to assume that the actual state could have been any state
                for (size_t c1 = 0; c1 < this->num_chars; ++c1)
                {
                    p_site[c1] = 1.0;
                }
            }
            else if (this->using_ambiguous_characters == true && this->using_weighted_characters == false)
            {
                const RbBitSet& val = amb_char_node[site];

                for (size_t c1 = 0; c1 < this->num_chars; ++c1)
                {
                    // transition probabilities from state c1 to each terminal state
                    const double* d = tp + c1 * this->num_chars;

                    double tmp = 0.0;
                    for (size_t i = 0; i < this->num_chars; ++i)
                    {
                        if (this->using_observation_error == false)
                        {
                            if (val.test(i) == true)
                            {
                                tmp += d[i];
                            }
                        }
                        else
                        {
                            double tmp2 = 0.0;
                            for (size_t j = 0; j < this->num_chars; ++j)
                            {
                                if (val.test(j) == true)
                                {
                                    if (i == j)
                                    {
                                        tmp2 += 1.0 - obs_error_prob * (1.0 - obs_error_freqs[i]);
                                    }
                                    else
                                    {
                                        tmp2 += obs_error_prob * obs_error_freqs[j];
                                    }
                                }
                            }
                            tmp += d[i] * tmp2;
                        }
                    }

                    p_site[c1] = tmp;
                }
            }
            else if (this->using_weighted_characters == true)
            {
                const size_t this_site_index = site_indices[site];
                const RbBitSet& val = this->value->getCharacter(char_data_node_index, this_site_index).getState();
                const std::vector<double>& weights = this->value->getCharacter(char_data_node_index, this_site_index).getWeights();

                for (size_t c1 = 0; c1 < this->num_chars; ++c1)
                {
                    const double* d = tp + c1 * this->num_chars;

                    double tmp = 0.0;
                    for (size_t i = 0; i < this->num_chars; ++i)
                    {
                        if (val.test(i) == true)
                        {
                            tmp += d[i] * weights[i];
                        }
                    }

                    p_site[c1] = tmp;
                }
            }
            else // no ambiguous characters in use
            {
                const std::uint64_t org_val = char_node[site];

                if (this->using_observation_error == false)
                {
                    for (size_t c1 = 0; c1 < this->num_chars; ++c1)
                    {
                        p_site[c1] = tp[c1 * this->num_chars + org_val];
                    }
                }
                else
                {
                    for (size_t c1 = 0; c1 < this->num_chars; ++c1)
                    {
                        double tmp = 0.0;
                        for (size_t c2 = 0; c2 < this->num_chars; ++c2)
                        {
                            if (c2 == org_val)
                            {
                                tmp += tp[c1 * this->num_chars + c2] *
                                    (1.0 - obs_error_prob * (1.0 - obs_error_freqs[c2]));
                            }
                            else
                            {
                                tmp += tp[c1 * this->num_chars + c2] *
                                    (obs_error_prob * obs_error_freqs[org_val]);
                            }
                        }
                        p_site[c1] = tmp;
                    }
                }
            }
        }
    }

    this->scale(node_index);
}


#endif
