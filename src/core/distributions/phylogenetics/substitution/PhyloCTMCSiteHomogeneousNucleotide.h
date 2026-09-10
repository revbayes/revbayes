#ifndef PhyloCTMCSiteHomogeneousNucleotide_H
#define PhyloCTMCSiteHomogeneousNucleotide_H

#include "AbstractPhyloCTMCSiteHomogeneous.h"
#include "DnaState.h"
#include "RateMatrix.h"
#include "RbVector.h"
#include "TopologyNode.h"
#include "TransitionProbabilityMatrix.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
    
    template<class charType>
    class PhyloCTMCSiteHomogeneousNucleotide : public AbstractPhyloCTMCSiteHomogeneous<charType> {
        
    public:
        PhyloCTMCSiteHomogeneousNucleotide(const TypedDagNode< Tree > *t, bool c, size_t nSites, bool amb, bool internal, bool gapmatch);
        virtual                                            ~PhyloCTMCSiteHomogeneousNucleotide(void);                                                                   //!< Virtual destructor
        
        // public member functions
        PhyloCTMCSiteHomogeneousNucleotide*                 clone(void) const;                                                                          //!< Create an independent clone
        
    protected:
        
        void                                                computeInternalNodeLikelihood(const TopologyNode &n, size_t nIdx, size_t l, size_t r);
        void                                                computeInternalNodeLikelihood(const TopologyNode &n, size_t nIdx, size_t l, size_t r,  size_t m);
        void                                                computeRootLikelihood( size_t root, size_t left, size_t right);
        void                                                computeRootLikelihood( size_t root, size_t left, size_t right, size_t middle);
        void                                                computeTipLikelihood(const TopologyNode &node, size_t nIdx);
        
        
    private:        
        
    };
    
}


#include "ConstantNode.h"
#include "DiscreteCharacterState.h"
#include "RateMatrix_JC.h"
#include "RandomNumberFactory.h"
#include "RbMathLogic.h"

#include <cmath>
#include <cstring>
#if defined (__AVX__)
#include <emmintrin.h>
#include <pmmintrin.h>
#include <immintrin.h>
#endif

template<class charType>
RevBayesCore::PhyloCTMCSiteHomogeneousNucleotide<charType>::PhyloCTMCSiteHomogeneousNucleotide(const TypedDagNode<Tree> *t, bool c, size_t nSites, bool amb, bool internal, bool gapmatch) : AbstractPhyloCTMCSiteHomogeneous<charType>(  t, 4, 1, c, nSites, amb, internal, gapmatch )
{
    
}

template<class charType>
RevBayesCore::PhyloCTMCSiteHomogeneousNucleotide<charType>::~PhyloCTMCSiteHomogeneousNucleotide( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
    
}


template<class charType>
RevBayesCore::PhyloCTMCSiteHomogeneousNucleotide<charType>* RevBayesCore::PhyloCTMCSiteHomogeneousNucleotide<charType>::clone( void ) const
{
    
    return new PhyloCTMCSiteHomogeneousNucleotide<charType>( *this );
}



template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneousNucleotide<charType>::computeRootLikelihood( size_t root, size_t left, size_t right)
{
    
    // reset the likelihood
    this->lnProb = 0.0;
    
    // get the root frequencies
    std::vector<std::vector<double> > ff;
    this->getRootFrequencies(ff);
    
    // get the pointers to the partial likelihoods of the left and right subtree
    auto& pl_left = this->getPartialLikelihoodsForNode(left);
    auto& pl_right= this->getPartialLikelihoodsForNode(right);
    const double* p_left   = pl_left.likelihoods.data();
    const double* p_right  = pl_right.likelihoods.data();
    assert(pl_left.dims() == pl_right.dims());
    
    double* p        = this->createEmptyPartialLikelihoodsForNode(root, pl_left.dims()).likelihoods.data();
    
    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
    {
        // get the root frequencies
        const std::vector<double> &f = ff[mixture % ff.size()];

        // get pointers to the likelihood for this mixture category
              double*   p_mixture          = p       + mixture * this->mixtureOffset;
        const double*   p_mixture_left     = p_left  + mixture * this->mixtureOffset;
        const double*   p_mixture_right    = p_right + mixture * this->mixtureOffset;

        // iterate over all sites
        for (size_t site = 0; site < this->pattern_block_size; ++site)
        {
            // get pointers to the likelihood for this site and mixture category
                  double*   p_site_mixture          = p_mixture       + site * this->siteOffset;
            const double*   p_site_mixture_left     = p_mixture_left  + site * this->siteOffset;
            const double*   p_site_mixture_right    = p_mixture_right + site * this->siteOffset;
            
            p_site_mixture[0] = p_site_mixture_left[0] * p_site_mixture_right[0] * f[0];
            p_site_mixture[1] = p_site_mixture_left[1] * p_site_mixture_right[1] * f[1];
            p_site_mixture[2] = p_site_mixture_left[2] * p_site_mixture_right[2] * f[2];
            p_site_mixture[3] = p_site_mixture_left[3] * p_site_mixture_right[3] * f[3];
            
        } // end-for over all sites (=patterns)
        
    } // end-for over all mixtures (=rate categories)

    this->scale( root, left, right );
}

template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneousNucleotide<charType>::computeRootLikelihood( size_t root, size_t left, size_t right, size_t middle)
{
    
    // reset the likelihood
    this->lnProb = 0.0;
    
    // get the root frequencies
    std::vector<std::vector<double> > ff;
    this->getRootFrequencies(ff);
    
    // get the pointers to the partial likelihoods of the left and right subtree
    auto& pl_left = this->getPartialLikelihoodsForNode(left);
    auto& pl_right = this->getPartialLikelihoodsForNode(right);
    auto& pl_middle = this->getPartialLikelihoodsForNode(middle);
    const double* p_left   = pl_left.likelihoods.data();
    const double* p_right  = pl_right.likelihoods.data();
    const double* p_middle = pl_middle.likelihoods.data();
    assert(pl_left.dims() == pl_right.dims());
    assert(pl_left.dims() == pl_middle.dims());

    double* p        = this->createEmptyPartialLikelihoodsForNode(root, pl_left.dims()).likelihoods.data();
    
    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
    {
        // get the root frequencies
        const std::vector<double> &f = ff[mixture % ff.size()];

        // get pointers to the likelihood for this mixture category
              double*   p_mixture          = p        + mixture * this->mixtureOffset;
        const double*   p_mixture_left     = p_left   + mixture * this->mixtureOffset;
        const double*   p_mixture_right    = p_right  + mixture * this->mixtureOffset;
        const double*   p_mixture_middle   = p_middle + mixture * this->mixtureOffset;

        // iterate over all sites
        for (size_t site = 0; site < this->pattern_block_size; ++site)
        {
            // get pointers to the likelihood for this site and mixture category
                  double*   p_site_mixture          = p_mixture        + site * this->siteOffset;
            const double*   p_site_mixture_left     = p_mixture_left   + site * this->siteOffset;
            const double*   p_site_mixture_right    = p_mixture_right  + site * this->siteOffset;
            const double*   p_site_mixture_middle   = p_mixture_middle + site * this->siteOffset;

            p_site_mixture[0] = p_site_mixture_left[0] * p_site_mixture_right[0] * p_site_mixture_middle[0] * f[0];
            p_site_mixture[1] = p_site_mixture_left[1] * p_site_mixture_right[1] * p_site_mixture_middle[1] * f[1];
            p_site_mixture[2] = p_site_mixture_left[2] * p_site_mixture_right[2] * p_site_mixture_middle[2] * f[2];
            p_site_mixture[3] = p_site_mixture_left[3] * p_site_mixture_right[3] * p_site_mixture_middle[3] * f[3];
            
        } // end-for over all sites (=patterns)
        
    } // end-for over all mixtures (=rate categories)
    
    this->scale( root, left, right, middle );
}

template <bool do_scaling, bool scale_this_branch>
void computeInternalNodeLikelihood4(PartialLikelihoods& pl_node, const PartialLikelihoods& pl_left, const PartialLikelihoods& pl_right,
                                    const std::vector<RevBayesCore::TransitionProbabilityMatrix>& pmatrices)
{
    constexpr int num_states = 4;
    constexpr int siteOffset = num_states;
    const int pattern_block_size = pl_left.dims().num_patterns;
    const int mixtureOffset = pattern_block_size * num_states;
    const int num_site_mixtures = pmatrices.size();

    std::vector<char> site_needs_scaling;
    if constexpr (do_scaling and scale_this_branch)
    {
        site_needs_scaling.resize(pattern_block_size,1);
    }

    // get the pointers to the partial likelihoods for this node and the two descendant subtrees
    assert(pl_left.dims() == pl_right.dims());

    const double* p_left   = pl_left.likelihoods.data();
    const double* p_right  = pl_right.likelihoods.data();
    auto& scale_left = pl_left.scale;
    auto& scale_right = pl_right.scale;

    double* p_node   = pl_node.likelihoods.data();
    auto& scale_node = pl_node.scale;

#   if defined(__AVX__)
    const __m256d scale_min_v = _mm256_set1_pd(scale_min);
#   endif

    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < num_site_mixtures; ++mixture)
    {
        // the transition probability matrix for this mixture category
        const double* tp_begin = pmatrices[mixture].theMatrix;
        
        // get the offset for this mixture category
        size_t offset = mixture*mixtureOffset;
        
#       if defined(__AVX__)
        
        __m256d tp_a = _mm256_loadu_pd(tp_begin);
        __m256d tp_c = _mm256_loadu_pd(tp_begin+4);
        __m256d tp_g = _mm256_loadu_pd(tp_begin+8);
        __m256d tp_t = _mm256_loadu_pd(tp_begin+12);
        
#       endif

        // compute the per site probabilities
        for (size_t site = 0; site < pattern_block_size ; ++site)
        {
            // get the pointers to the likelihood for this site and mixture category
            size_t site_offset = offset + site*siteOffset;

            double*          p_site_mixture          = p_node  + site_offset;
            const double*    p_site_mixture_left     = p_left  + site_offset;
            const double*    p_site_mixture_right    = p_right + site_offset;
            
#           if defined(__AVX__)
 
            __m256d a = _mm256_loadu_pd(p_site_mixture_left);
            __m256d b = _mm256_loadu_pd(p_site_mixture_right);
            __m256d p = _mm256_mul_pd(a,b);
            
            __m256d a_acgt = _mm256_mul_pd(p, tp_a );
            __m256d c_acgt = _mm256_mul_pd(p, tp_c );
            __m256d g_acgt = _mm256_mul_pd(p, tp_g );
            __m256d t_acgt = _mm256_mul_pd(p, tp_t );
            
            __m256d ac   = _mm256_hadd_pd(a_acgt,c_acgt);
            __m256d gt   = _mm256_hadd_pd(g_acgt,t_acgt);

            __m256d lo  = _mm256_permute2f128_pd(ac, gt, 0x20);
            __m256d hi  = _mm256_permute2f128_pd(ac, gt, 0x31);
            __m256d sum = _mm256_add_pd(lo, hi);

            _mm256_storeu_pd(p_site_mixture, sum);

            if constexpr (do_scaling and scale_this_branch)
            {
                __m256d cmp = _mm256_cmp_pd(sum, scale_min_v, _CMP_GE_OQ);
                int mask = _mm256_movemask_pd(cmp);
                site_needs_scaling[site] &= (mask == 0);
            }

#           else

            double p0 = p_site_mixture_left[0] * p_site_mixture_right[0];
            double p1 = p_site_mixture_left[1] * p_site_mixture_right[1];
            double p2 = p_site_mixture_left[2] * p_site_mixture_right[2];
            double p3 = p_site_mixture_left[3] * p_site_mixture_right[3];
            
            double sum = p0 * tp_begin[0];
            sum += p1 * tp_begin[1];
            sum += p2 * tp_begin[2];
            sum += p3 * tp_begin[3];
            
            p_site_mixture[0] = sum;
            
            sum = p0 * tp_begin[4];
            sum += p1 * tp_begin[5];
            sum += p2 * tp_begin[6];
            sum += p3 * tp_begin[7];
            
            p_site_mixture[1] = sum;
            
            sum = p0 * tp_begin[8];
            sum += p1 * tp_begin[9];
            sum += p2 * tp_begin[10];
            sum += p3 * tp_begin[11];
            
            p_site_mixture[2] = sum;
            
            sum = p0 * tp_begin[12];
            sum += p1 * tp_begin[13];
            sum += p2 * tp_begin[14];
            sum += p3 * tp_begin[15];
            
            p_site_mixture[3] = sum;

            if constexpr (do_scaling and scale_this_branch)
            {
                bool mixture_needs_scaling = (p_site_mixture[0] < scale_min) && (p_site_mixture[1] < scale_min) && (p_site_mixture[2] < scale_min) && (p_site_mixture[3] < scale_min);
                site_needs_scaling[site] &= mixture_needs_scaling;
            }

#           endif
            
        } // end-for over all sites (=patterns)
        
    } // end-for over all mixtures (=rate-categories)

    // Rescale the sites, but only if (i) rescaling is enabled and (ii) they need it.
    if constexpr ( do_scaling and scale_this_branch )
    {
        // iterate over all mixture categories
        for (size_t site = 0; site < pattern_block_size ; ++site)
        {
            scale_node[site] = scale_left[site] + scale_right[site];
            if (site_needs_scaling[site])
                scale_node[site]++;
        }

        for (size_t mixture = 0; mixture < num_site_mixtures; ++mixture)
        {
            for (size_t site = 0; site < pattern_block_size ; ++site)
            {
                if (site_needs_scaling[site])
                {
                    // get the pointers to the likelihood for this mixture category
                    double* p_site_mixture = p_node + mixture*mixtureOffset + site*siteOffset;
                    
                    for ( size_t i=0; i<4; ++i)
                        p_site_mixture[i] *= scale_factor;
                }
            }
        }
    }
    else if constexpr ( do_scaling )
    {
        // iterate over all mixture categories
        for (size_t site = 0; site < pattern_block_size ; ++site)
            scale_node[site] = scale_left[site] + scale_right[site];
    }
    
}


template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneousNucleotide<charType>::computeInternalNodeLikelihood(const TopologyNode &node, size_t node_index, size_t left, size_t right) 
{
    // update the transition probability matrix
    this->updateTransitionProbabilityMatrix( node_index );

    // get the pointers to the partial likelihoods for this node and the two descendant subtrees
    auto& pl_left = this->getPartialLikelihoodsForNode(left);
    auto& pl_right = this->getPartialLikelihoodsForNode(right);
    auto& pl_node = this->createEmptyPartialLikelihoodsForNode(node_index, pl_left.dims());

    auto& pmatrices = this->pmatrices[node_index];

    bool do_scaling = RbSettings::userSettings().getUseScaling();
    bool scale_this_branch = node_index % RbSettings::userSettings().getScalingDensity() == 0;

    if (do_scaling and scale_this_branch)
        computeInternalNodeLikelihood4<true,true>(pl_node, pl_left, pl_right, pmatrices);
    else if (do_scaling and not scale_this_branch)
        computeInternalNodeLikelihood4<true,false>(pl_node, pl_left, pl_right, pmatrices);
    else
        computeInternalNodeLikelihood4<false,false>(pl_node, pl_left, pl_right, pmatrices);
}


template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneousNucleotide<charType>::computeInternalNodeLikelihood(const TopologyNode &node, size_t node_index, size_t left, size_t right, size_t middle)
{
    
    // update the transition probability matrix
    this->updateTransitionProbabilityMatrix( node_index );
    
    // get the pointers to the partial likelihoods for this node and the two descendant subtrees
    auto& pl_left = this->getPartialLikelihoodsForNode(left);
    auto& pl_middle = this->getPartialLikelihoodsForNode(middle);
    auto& pl_right = this->getPartialLikelihoodsForNode(right);
    const double*   p_left      = pl_left.likelihoods.data();
    const double*   p_middle    = pl_middle.likelihoods.data();
    const double*   p_right     = pl_right.likelihoods.data();
    assert(pl_left.dims() == pl_middle.dims());
    assert(pl_left.dims() == pl_right.dims());
    
    double*         p_node      = this->createEmptyPartialLikelihoodsForNode(node_index, pl_left.dims()).likelihoods.data();

    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
    {
        // the transition probability matrix for this mixture category
//         const double* tp_begin = this->transition_prob_matrices[mixture].theMatrix;
        const double* tp_begin = this->pmatrices[node_index][mixture].theMatrix;
        
        // get the offset for this mixture category
        size_t offset = mixture*this->mixtureOffset;
        
#       if defined(__AVX__)
        
        __m256d tp_a = _mm256_loadu_pd(tp_begin);
        __m256d tp_c = _mm256_loadu_pd(tp_begin+4);
        __m256d tp_g = _mm256_loadu_pd(tp_begin+8);
        __m256d tp_t = _mm256_loadu_pd(tp_begin+12);
        
#       endif
        
        // compute the per site probabilities
        for (size_t site = 0; site < this->pattern_block_size ; ++site)
        {
            // get the pointers to the likelihood for this site and mixture category
            size_t site_offset = offset + site*this->siteOffset;

            double*          p_site_mixture          = p_node   + site_offset;
            const double*    p_site_mixture_left     = p_left   + site_offset;
            const double*    p_site_mixture_middle   = p_middle + site_offset;
            const double*    p_site_mixture_right    = p_right  + site_offset;
            
#           if defined (__AVX__)
            
            __m256d a = _mm256_loadu_pd(p_site_mixture_left);
            __m256d b = _mm256_loadu_pd(p_site_mixture_middle);
            __m256d c = _mm256_loadu_pd(p_site_mixture_right);
            __m256d p = _mm256_mul_pd(_mm256_mul_pd(a, b), c);

            __m256d a_acgt = _mm256_mul_pd(p, tp_a);
            __m256d c_acgt = _mm256_mul_pd(p, tp_c);
            __m256d g_acgt = _mm256_mul_pd(p, tp_g);
            __m256d t_acgt = _mm256_mul_pd(p, tp_t);

            // First-level reduction: pair-sum within each 128-bit lane.
            // ac = [a0+a1, c0+c1, a2+a3, c2+c3]
            // gt = [g0+g1, t0+t1, g2+g3, t2+t3]
            __m256d ac = _mm256_hadd_pd(a_acgt, c_acgt);
            __m256d gt = _mm256_hadd_pd(g_acgt, t_acgt);

            // Second-level reduction: cross the 128-bit lane boundary.
            // lo  = [a0+a1, c0+c1, g0+g1, t0+t1]
            // hi  = [a2+a3, c2+c3, g2+g3, t2+t3]
            // sum = [a0..3,  c0..3,  g0..3,  t0..3]
            __m256d lo  = _mm256_permute2f128_pd(ac, gt, 0x20);
            __m256d hi  = _mm256_permute2f128_pd(ac, gt, 0x31);
            __m256d sum = _mm256_add_pd(lo, hi);

            _mm256_storeu_pd(p_site_mixture, sum);
            
#           else
            
            double p0 = p_site_mixture_left[0] * p_site_mixture_middle[0] * p_site_mixture_right[0];
            double p1 = p_site_mixture_left[1] * p_site_mixture_middle[1] * p_site_mixture_right[1];
            double p2 = p_site_mixture_left[2] * p_site_mixture_middle[2] * p_site_mixture_right[2];
            double p3 = p_site_mixture_left[3] * p_site_mixture_middle[3] * p_site_mixture_right[3];
            
            double sum = p0 * tp_begin[0];
            sum += p1 * tp_begin[1];
            sum += p2 * tp_begin[2];
            sum += p3 * tp_begin[3];
            
            p_site_mixture[0] = sum;
            
            sum = p0 * tp_begin[4];
            sum += p1 * tp_begin[5];
            sum += p2 * tp_begin[6];
            sum += p3 * tp_begin[7];
            
            p_site_mixture[1] = sum;
            
            sum = p0 * tp_begin[8];
            sum += p1 * tp_begin[9];
            sum += p2 * tp_begin[10];
            sum += p3 * tp_begin[11];
            
            p_site_mixture[2] = sum;
            
            sum = p0 * tp_begin[12];
            sum += p1 * tp_begin[13];
            sum += p2 * tp_begin[14];
            sum += p3 * tp_begin[15];
            
            p_site_mixture[3] = sum;
            
#           endif
            
        } // end-for over all sites (=patterns)
        
    } // end-for over all mixtures (=rate-categories)
    
}


template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneousNucleotide<charType>::computeTipLikelihood(const TopologyNode &node, size_t node_index) 
{    
    
    double* p_node = this->createEmptyPartialLikelihoodsForNode(node_index, {this->num_site_mixtures, this->pattern_block_size, this->num_chars}).likelihoods.data();
    
    size_t data_tip_index = this->taxon_name_2_tip_index_map[ node.getName() ];
    const std::vector<bool> &gap_node = this->gap_matrix[data_tip_index];
    const std::vector<std::uint64_t> &char_node = this->char_matrix[data_tip_index];
    const std::vector<RbBitSet> &amb_char_node = this->ambiguous_char_matrix[data_tip_index];
    
    // update the transition probabilities
    this->updateTransitionProbabilityMatrix( node_index );
    
    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
    {
        // the transition probability matrix for this mixture category
//         const double*       tp_begin    = this->transition_prob_matrices[mixture].theMatrix;
        const double*       tp_begin    = this->pmatrices[node_index][mixture].theMatrix;
        
        // get the pointer to the likelihoods for this mixture category
        double*     p_mixture      = p_node + mixture * this->mixtureOffset;
        
        // iterate over all sites
        for (size_t site = 0; site < this->pattern_block_size; ++site)
        {
            // get the pointer to the likelihoods for this site and mixture category
            double*     p_site_mixture      = p_mixture + site * this->siteOffset;
            
            // is this site a gap?
            if ( gap_node[site] ) 
            {
                // since this is a gap we need to assume that the actual state could have been any state
                p_site_mixture[0] = 1.0;
                p_site_mixture[1] = 1.0;
                p_site_mixture[2] = 1.0;
                p_site_mixture[3] = 1.0;
                
            } 
            else // we have observed a character
            {
                                    
                if ( this->using_ambiguous_characters == true )
                {
                    // get the original character
                    const RbBitSet &org_val = amb_char_node[site];
                    
                    double p0 = 0.0;
                    double p1 = 0.0;
                    double p2 = 0.0;
                    double p3 = 0.0;
                    
                    if ( org_val.test(0) == true )
                    {
                        p0 = tp_begin[0];
                        p1 = tp_begin[4];
                        p2 = tp_begin[8];
                        p3 = tp_begin[12];
                    }
                    
                    if ( org_val.test(1) == true )
                    {
                        p0 += tp_begin[1];
                        p1 += tp_begin[5];
                        p2 += tp_begin[9];
                        p3 += tp_begin[13];
                    }
                    
                    if ( org_val.test(2) == true )
                    {
                        p0 += tp_begin[2];
                        p1 += tp_begin[6];
                        p2 += tp_begin[10];
                        p3 += tp_begin[14];
                    }
                    
                    if ( org_val.test(3) == true )
                    {
                        p0 += tp_begin[3];
                        p1 += tp_begin[7];
                        p2 += tp_begin[11];
                        p3 += tp_begin[15];
                    }
                    
                    p_site_mixture[0] = p0;
                    p_site_mixture[1] = p1;
                    p_site_mixture[2] = p2;
                    p_site_mixture[3] = p3;
                    
                } 
                else // no ambiguous characters in use
                {
                    
                    // get the original character
                    std::uint64_t org_val = char_node[site];
                    
                    // store the likelihood
                    p_site_mixture[0] = tp_begin[org_val];
                    p_site_mixture[1] = tp_begin[4+org_val];
                    p_site_mixture[2] = tp_begin[8+org_val];
                    p_site_mixture[3] = tp_begin[12+org_val];
                        
                }
                
            } // end-if a gap state
            
        } // end-for over all sites/patterns in the sequence
        
    } // end-for over all mixture categories

    this->scale( node_index );
}


#endif
