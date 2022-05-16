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
        
        void                                                computeInternalNodeLikelihood(const TopologyNode &node, size_t node_index, size_t left, size_t right, double* likelihoods, size_t likelihood_offset, size_t node_offset, size_t pattern_block_size, size_t mixture_offset);
        void                                                computeInternalNodeLikelihood(const TopologyNode &node, size_t node_index, size_t left, size_t right, size_t middle, double* likelihoods, size_t likelihood_offset, size_t node_offset, size_t pattern_block_size, size_t mixture_offset);
        void                                                computeRootLikelihood( size_t root, size_t left, size_t right, double* likelihoods, size_t likelihood_offset, size_t node_offset, size_t pattern_block_size, size_t mixture_offset);
        void                                                computeRootLikelihood( size_t root, size_t left, size_t right, size_t middle, double* likelihoods, size_t likelihood_offset, size_t node_offset, size_t pattern_block_size, size_t mixture_offset);
        void                                                computeTipLikelihood(const TopologyNode &node, size_t node_index, double* likelihoods, size_t likelihood_offset,  size_t node_offset, size_t pattern_block_size, size_t mixture_offset, const std::vector<std::vector<RbBitSet> >& ambiguous_char_matrix, const std::vector<std::vector<unsigned long> >& char_matrix, const std::vector<std::vector<bool> >& gap_matrix);
        
        
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
#if defined (SSE_ENABLED)
#include <xmmintrin.h>
#include <emmintrin.h>
#include <pmmintrin.h>
#elif defined ( AVX_ENABLED )
#include <xmmintrin.h>
#include <emmintrin.h>
#include <pmmintrin.h>
#include <immintrin.h>
#endif

template<class charType>
RevBayesCore::PhyloCTMCSiteHomogeneousNucleotide<charType>::PhyloCTMCSiteHomogeneousNucleotide(const TypedDagNode<Tree> *t, bool c, size_t nSites, bool amb, bool internal, bool gapmatch) : AbstractPhyloCTMCSiteHomogeneous<charType>(  t, 4, 1, c, nSites, amb, false, internal, gapmatch )
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
void RevBayesCore::PhyloCTMCSiteHomogeneousNucleotide<charType>::computeRootLikelihood( size_t root, size_t left, size_t right, double* likelihoods, size_t likelihood_offset, size_t node_offset, size_t pattern_block_size, size_t mixture_offset)
{

    // reset the likelihood
    this->ln_prob = 0.0;
    
    // get the root frequencies
    std::vector<std::vector<double> > ff;
    this->getRootFrequencies(ff);
    
    // get the pointers to the partial likelihoods of the left and right subtree
          double* p        = likelihoods + this->active_likelihood[root]  * likelihood_offset + root   * node_offset;
    const double* p_left   = likelihoods + this->active_likelihood[left]  * likelihood_offset + left   * node_offset;
    const double* p_right  = likelihoods + this->active_likelihood[right] * likelihood_offset + right  * node_offset;
    
    // get pointers the likelihood for both subtrees
          double*   p_mixture          = p;
    const double*   p_mixture_left     = p_left;
    const double*   p_mixture_right    = p_right;
    
    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
    {
        // get the root frequencies
        const std::vector<double> &f = ff[mixture % ff.size()];

        // get pointers to the likelihood for this mixture category
              double*   p_site_mixture          = p_mixture;
        const double*   p_site_mixture_left     = p_mixture_left;
        const double*   p_site_mixture_right    = p_mixture_right;
        // iterate over all sites
        for (size_t site = 0; site < pattern_block_size; ++site)
        {
            
            p_site_mixture[0] = p_site_mixture_left[0] * p_site_mixture_right[0] * f[0];
            p_site_mixture[1] = p_site_mixture_left[1] * p_site_mixture_right[1] * f[1];
            p_site_mixture[2] = p_site_mixture_left[2] * p_site_mixture_right[2] * f[2];
            p_site_mixture[3] = p_site_mixture_left[3] * p_site_mixture_right[3] * f[3];
            
            // increment the pointers to the next site
            p_site_mixture+=this->site_offset; p_site_mixture_left+=this->site_offset; p_site_mixture_right+=this->site_offset;
            
        } // end-for over all sites (=patterns)
        
        // increment the pointers to the next mixture category
        p_mixture+=mixture_offset; p_mixture_left+=mixture_offset; p_mixture_right+=mixture_offset;
        
    } // end-for over all mixtures (=rate categories)
    
}


template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneousNucleotide<charType>::computeRootLikelihood( size_t root, size_t left, size_t right, size_t middle, double* likelihoods, size_t likelihood_offset, size_t node_offset, size_t pattern_block_size, size_t mixture_offset)
{
    
    // reset the likelihood
    this->ln_prob = 0.0;
    
    // get the root frequencies
    std::vector<std::vector<double> > ff;
    this->getRootFrequencies(ff);
    
    // get the pointers to the partial likelihoods of the left and right subtree
          double* p        = likelihoods + this->active_likelihood[root]  * likelihood_offset + root   * node_offset;
    const double* p_left   = likelihoods + this->active_likelihood[left]  * likelihood_offset + left   * node_offset;
    const double* p_right  = likelihoods + this->active_likelihood[right] * likelihood_offset + right  * node_offset;
    const double* p_middle = likelihoods + this->active_likelihood[middle]* likelihood_offset + middle * node_offset;
    
    // get pointers the likelihood for both subtrees
          double*   p_mixture          = p;
    const double*   p_mixture_left     = p_left;
    const double*   p_mixture_right    = p_right;
    const double*   p_mixture_middle   = p_middle;
    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
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
            p_site_mixture[0] = p_site_mixture_left[0] * p_site_mixture_right[0] * p_site_mixture_middle[0] * f[0];
            p_site_mixture[1] = p_site_mixture_left[1] * p_site_mixture_right[1] * p_site_mixture_middle[1] * f[1];
            p_site_mixture[2] = p_site_mixture_left[2] * p_site_mixture_right[2] * p_site_mixture_middle[2] * f[2];
            p_site_mixture[3] = p_site_mixture_left[3] * p_site_mixture_right[3] * p_site_mixture_middle[3] * f[3];
            
            // increment the pointers to the next site
            p_site_mixture+=this->site_offset; p_site_mixture_left+=this->site_offset; p_site_mixture_right+=this->site_offset; p_site_mixture_middle+=this->site_offset;
            
        } // end-for over all sites (=patterns)
        
        // increment the pointers to the next mixture category
        p_mixture+=mixture_offset; p_mixture_left+=mixture_offset; p_mixture_right+=mixture_offset; p_mixture_middle+=mixture_offset;
        
    } // end-for over all mixtures (=rate categories)
    
}



template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneousNucleotide<charType>::computeInternalNodeLikelihood(const TopologyNode &node, size_t node_index, size_t left, size_t right, double* likelihoods, size_t likelihood_offset, size_t node_offset, size_t pattern_block_size, size_t mixture_offset)
{
    
    // compute the transition probability matrix
    this->updateTransitionProbabilities( node_index );
    
#   if defined ( SSE_ENABLED )
    
    double* p_left   = likelihoods + this->active_likelihood[left]       * likelihood_offset + left          * node_offset;
    double* p_right  = likelihoods + this->active_likelihood[right]      * likelihood_offset + right         * node_offset;
    double* p_node   = likelihoods + this->active_likelihood[node_index] * likelihood_offset + node_index    * node_offset;
    //    __m128d* p_left   = (__m128d *) this->partialLikelihoods + this->active_likelihood[left]*this->active_likelihoodOffset + left*this->nodeOffset;
    //    __m128d* p_right  = (__m128d *) this->partialLikelihoods + this->active_likelihood[right]*this->active_likelihoodOffset + right*this->nodeOffset;
    //    __m128d* p_node   = (__m128d *) this->partialLikelihoods + this->active_likelihood[node_index]*this->active_likelihoodOffset + node_index*this->nodeOffset;
    
#   elif defined ( AVX_ENABLED )

    double* p_left   = likelihoods + this->active_likelihood[left]       * likelihood_offset + left          * node_offset;
    double* p_right  = likelihoods + this->active_likelihood[right]      * likelihood_offset + right         * node_offset;
    double* p_node   = likelihoods + this->active_likelihood[node_index] * likelihood_offset + node_index    * node_offset;

    double* tmp_ac = new double[4];
    double* tmp_gt = new double[4];
//    double tmp_ac[4];
//    double tmp_gt[4];
    
    
#   else

    // get the pointers to the partial likelihoods for this node and the two descendant subtrees
    double* p_left   = likelihoods + this->active_likelihood[left]       * likelihood_offset + left          * node_offset;
    double* p_right  = likelihoods + this->active_likelihood[right]      * likelihood_offset + right         * node_offset;
    double* p_node   = likelihoods + this->active_likelihood[node_index] * likelihood_offset + node_index    * node_offset;
    
#   endif
    
    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
    {
        // the transition probability matrix for this mixture category
        const double* tp_begin = this->transition_prob_matrices[mixture].theMatrix;
        
        // get the pointers to the likelihood for this mixture category
        size_t offset = mixture*mixture_offset;
        
#       if defined ( SSE_ENABLED )
        
        double*          p_site_mixture          = p_node  + offset;
        const double*    p_site_mixture_left     = p_left  + offset;
        const double*    p_site_mixture_right    = p_right + offset;
        
        __m128d tp_a_ac = _mm_load_pd(tp_begin);
        __m128d tp_a_gt = _mm_load_pd(tp_begin+2);
        __m128d tp_c_ac = _mm_load_pd(tp_begin+4);
        __m128d tp_c_gt = _mm_load_pd(tp_begin+6);
        __m128d tp_g_ac = _mm_load_pd(tp_begin+8);
        __m128d tp_g_gt = _mm_load_pd(tp_begin+10);
        __m128d tp_t_ac = _mm_load_pd(tp_begin+12);
        __m128d tp_t_gt = _mm_load_pd(tp_begin+14);
        
#       elif defined ( AVX_ENABLED )
        
        double*          p_site_mixture          = p_node  + offset;
        const double*    p_site_mixture_left     = p_left  + offset;
        const double*    p_site_mixture_right    = p_right + offset;
        
        __m256d tp_a = _mm256_load_pd(tp_begin);
        __m256d tp_c = _mm256_load_pd(tp_begin+4);
        __m256d tp_g = _mm256_load_pd(tp_begin+8);
        __m256d tp_t = _mm256_load_pd(tp_begin+12);
        
#       else

        double*          p_site_mixture          = p_node  + offset;
        const double*    p_site_mixture_left     = p_left  + offset;
        const double*    p_site_mixture_right    = p_right + offset;

#       endif

        // compute the per site probabilities
        for (size_t site = 0; site < pattern_block_size ; ++site)
        {
            
#           if defined ( SSE_ENABLED )
            
            __m128d a01 = _mm_load_pd(p_site_mixture_left);
            __m128d a23 = _mm_load_pd(p_site_mixture_left+2);
            
            __m128d b01 = _mm_load_pd(p_site_mixture_right);
            __m128d b23 = _mm_load_pd(p_site_mixture_right+2);
            
            __m128d p01 = _mm_mul_pd(a01,b01);
            __m128d p23 = _mm_mul_pd(a23,b23);
            
            __m128d a_ac = _mm_mul_pd(p01, tp_a_ac   );
            __m128d a_gt = _mm_mul_pd(p23, tp_a_gt );
            __m128d a_acgt = _mm_hadd_pd(a_ac,a_gt);
            
            __m128d c_ac = _mm_mul_pd(p01, tp_c_ac );
            __m128d c_gt = _mm_mul_pd(p23, tp_c_gt );
            __m128d c_acgt = _mm_hadd_pd(c_ac,c_gt);
            
            __m128d ac = _mm_hadd_pd(a_acgt,c_acgt);
            _mm_store_pd(p_site_mixture,ac);
            
            
            __m128d g_ac = _mm_mul_pd(p01, tp_g_ac  );
            __m128d g_gt = _mm_mul_pd(p23, tp_g_gt );
            __m128d g_acgt = _mm_hadd_pd(g_ac,g_gt);
            
            __m128d t_ac = _mm_mul_pd(p01, tp_t_ac );
            __m128d t_gt = _mm_mul_pd(p23, tp_t_gt );
            __m128d t_acgt = _mm_hadd_pd(t_ac,t_gt);
            
            __m128d gt = _mm_hadd_pd(g_acgt,t_acgt);
            _mm_store_pd(p_site_mixture+2,gt);
 
#           elif defined ( AVX_ENABLED )
 
            __m256d a = _mm256_load_pd(p_site_mixture_left);
            __m256d b = _mm256_load_pd(p_site_mixture_right);
            __m256d p = _mm256_mul_pd(a,b);
            
            __m256d a_acgt = _mm256_mul_pd(p, tp_a );
            __m256d c_acgt = _mm256_mul_pd(p, tp_c );
            __m256d g_acgt = _mm256_mul_pd(p, tp_g );
            __m256d t_acgt = _mm256_mul_pd(p, tp_t );
            
            __m256d ac   = _mm256_hadd_pd(a_acgt,c_acgt);
            __m256d gt   = _mm256_hadd_pd(g_acgt,t_acgt);
            
            
            _mm256_store_pd(tmp_ac,ac);
            _mm256_store_pd(tmp_gt,gt);
            
            p_site_mixture[0] = tmp_ac[0] + tmp_ac[2];
            p_site_mixture[1] = tmp_ac[1] + tmp_ac[3];
            p_site_mixture[2] = tmp_gt[0] + tmp_gt[2];
            p_site_mixture[3] = tmp_gt[1] + tmp_gt[3];

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

#           endif
            
            // increment the pointers to the next site
            p_site_mixture_left+=this->site_offset; p_site_mixture_right+=this->site_offset; p_site_mixture+=this->site_offset;

                        
        } // end-for over all sites (=patterns)
        
    } // end-for over all mixtures (=rate-categories)

# if defined ( AVX_ENABLED )
    delete[] tmp_ac;
    delete[] tmp_gt;
# endif
    
}

template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneousNucleotide<charType>::computeInternalNodeLikelihood(const TopologyNode &node, size_t node_index, size_t left, size_t right, size_t middle, double* likelihoods, size_t likelihood_offset, size_t node_offset, size_t pattern_block_size, size_t mixture_offset)
{
    
    // compute the transition probability matrix
    this->updateTransitionProbabilities( node_index );
    
    
    // get the pointers to the partial likelihoods for this node and the two descendant subtrees
    const double*   p_left      = likelihoods + this->active_likelihood[left]        * likelihood_offset + left          * node_offset;
    const double*   p_middle    = likelihoods + this->active_likelihood[middle]      * likelihood_offset + middle        * node_offset;
    const double*   p_right     = likelihoods + this->active_likelihood[right]       * likelihood_offset + right         * node_offset;
    double*         p_node      = likelihoods + this->active_likelihood[node_index]  * likelihood_offset + node_index    * node_offset;
    
    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
    {
        // the transition probability matrix for this mixture category
        const double* tp_begin = this->transition_prob_matrices[mixture].theMatrix;
        
        // get the pointers to the likelihood for this mixture category
        size_t offset = mixture * mixture_offset;
        
#       if defined ( SSE_ENABLED )
        
        double*          p_site_mixture          = p_node   + offset;
        const double*    p_site_mixture_left     = p_left   + offset;
        const double*    p_site_mixture_middle   = p_middle + offset;
        const double*    p_site_mixture_right    = p_right  + offset;
        
        __m128d tp_a_ac = _mm_load_pd(tp_begin);
        __m128d tp_a_gt = _mm_load_pd(tp_begin+2);
        __m128d tp_c_ac = _mm_load_pd(tp_begin+4);
        __m128d tp_c_gt = _mm_load_pd(tp_begin+6);
        __m128d tp_g_ac = _mm_load_pd(tp_begin+8);
        __m128d tp_g_gt = _mm_load_pd(tp_begin+10);
        __m128d tp_t_ac = _mm_load_pd(tp_begin+12);
        __m128d tp_t_gt = _mm_load_pd(tp_begin+14);
        
#       elif defined ( AVX_ENABLED )
        
        double*          p_site_mixture          = p_node   + offset;
        const double*    p_site_mixture_left     = p_left   + offset;
        const double*    p_site_mixture_right    = p_right  + offset;
        
        __m256d tp_a = _mm256_load_pd(tp_begin);
        __m256d tp_c = _mm256_load_pd(tp_begin+4);
        __m256d tp_g = _mm256_load_pd(tp_begin+8);
        __m256d tp_t = _mm256_load_pd(tp_begin+12);
        
#       else
        
        double*          p_site_mixture          = p_node   + offset;
        const double*    p_site_mixture_left     = p_left   + offset;
        const double*    p_site_mixture_middle   = p_middle + offset;
        const double*    p_site_mixture_right    = p_right  + offset;
        
#       endif
        
        // compute the per site probabilities
        for (size_t site = 0; site < pattern_block_size ; ++site)
        {
            
#           if defined ( SSE_ENABLED )
            
            __m128d a01 = _mm_load_pd(p_site_mixture_left);
            __m128d a23 = _mm_load_pd(p_site_mixture_left+2);
            
            __m128d b01 = _mm_load_pd(p_site_mixture_middle);
            __m128d b23 = _mm_load_pd(p_site_mixture_middle+2);
            
            __m128d c01 = _mm_load_pd(p_site_mixture_right);
            __m128d c23 = _mm_load_pd(p_site_mixture_right+2);
            
            __m128d tmp_p01 = _mm_mul_pd(a01,b01);
            __m128d p01 = _mm_mul_pd(tmp_p01,c01);
            __m128d tmp_p23 = _mm_mul_pd(a23,b23);
            __m128d p23 = _mm_mul_pd(tmp_p23,c23);
            
            __m128d a_ac = _mm_mul_pd(p01, tp_a_ac   );
            __m128d a_gt = _mm_mul_pd(p23, tp_a_gt );
            __m128d a_acgt = _mm_hadd_pd(a_ac,a_gt);
            
            __m128d c_ac = _mm_mul_pd(p01, tp_c_ac );
            __m128d c_gt = _mm_mul_pd(p23, tp_c_gt );
            __m128d c_acgt = _mm_hadd_pd(c_ac,c_gt);
            

            //            *p_site_mixture = _mm_hadd_pd(a_acgt,c_acgt);
            __m128d ac = _mm_hadd_pd(a_acgt,c_acgt);
            _mm_store_pd(p_site_mixture,ac);
            
            
            __m128d g_ac = _mm_mul_pd(p01, tp_g_ac  );
            __m128d g_gt = _mm_mul_pd(p23, tp_g_gt );
            __m128d g_acgt = _mm_hadd_pd(g_ac,g_gt);
            
            __m128d t_ac = _mm_mul_pd(p01, tp_t_ac );
            __m128d t_gt = _mm_mul_pd(p23, tp_t_gt );
            __m128d t_acgt = _mm_hadd_pd(t_ac,t_gt);
            
            //            p_site_mixture[2] = _mm_hadd_pd(g_acgt,t_acgt);
            __m128d gt = _mm_hadd_pd(g_acgt,t_acgt);
            _mm_store_pd(p_site_mixture+2,gt);
            
#           elif defined ( AVX_ENABLED )
            
            __m256d a = _mm256_load_pd(p_site_mixture_left);
            __m256d b = _mm256_load_pd(p_site_mixture_right);
            __m256d p = _mm_mul_pd(a,b);
            
            __m256d a_acgt = _mm256_mul_pd(p, tp_a );
            __m256d c_acgt = _mm256_mul_pd(p, tp_c );
            __m256d g_acgt = _mm256_mul_pd(p, tp_g );
            __m256d t_acgt = _mm256_mul_pd(p, tp_t );
            
            __m256d ac   = _mm256_hadd_pd(a_acgt,c_acgt);
            __m256d gt   = _mm256_hadd_pd(g_acgt,t_acgt)
            
            __m256d acgt = _mm256_hadd_pd(ac,gt);
            
            _mm256_store_pd(p_site_mixture,acgt);
            
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
            
            // increment the pointers to the next site
            p_site_mixture_left+=this->site_offset; p_site_mixture_middle+=this->site_offset; p_site_mixture_right+=this->site_offset; p_site_mixture+=this->site_offset;
            
            
        } // end-for over all sites (=patterns)
        
    } // end-for over all mixtures (=rate-categories)
    
}


template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneousNucleotide<charType>::computeTipLikelihood(const TopologyNode &node, size_t node_index, double* likelihoods, size_t likelihood_offset,  size_t node_offset, size_t pattern_block_size, size_t mixture_offset, const std::vector<std::vector<RbBitSet> >& ambiguous_char_matrix, const std::vector<std::vector<unsigned long> >& char_matrix, const std::vector<std::vector<bool> >& gap_matrix)
{    
    
    double* p_node = likelihoods + this->active_likelihood[node_index]*likelihood_offset + node_index*node_offset;
    
    size_t data_tip_index = this->taxon_name_2_tip_index_map[ node.getName() ];
    const std::vector<bool> &gap_node           = gap_matrix[data_tip_index];
    const std::vector<unsigned long> &char_node = char_matrix[data_tip_index];
    const std::vector<RbBitSet> &amb_char_node  = ambiguous_char_matrix[data_tip_index];
    
    // compute the transition probabilities
    this->updateTransitionProbabilities( node_index );
    
    double*   p_mixture      = p_node;
    
    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
    {
        // the transition probability matrix for this mixture category
        const double*       tp_begin    = this->transition_prob_matrices[mixture].theMatrix;
        
        // get the pointer to the likelihoods for this site and mixture category
        double*     p_site_mixture      = p_mixture;
        
        // iterate over all sites
        for (size_t site = 0; site < pattern_block_size; ++site)
        {
            
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
                    
                    if ( org_val.isSet(0) == true )
                    {
                        p0 = tp_begin[0];
                        p1 = tp_begin[4];
                        p2 = tp_begin[8];
                        p3 = tp_begin[12];
                    }
                    
                    if ( org_val.isSet(1) == true )
                    {
                        p0 += tp_begin[1];
                        p1 += tp_begin[5];
                        p2 += tp_begin[9];
                        p3 += tp_begin[13];
                    }
                    
                    if ( org_val.isSet(2) == true )
                    {
                        p0 += tp_begin[2];
                        p1 += tp_begin[6];
                        p2 += tp_begin[10];
                        p3 += tp_begin[14];
                    }
                    
                    if ( org_val.isSet(3) == true )
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
                    unsigned long org_val = char_node[site];
                    
                    // store the likelihood
                    p_site_mixture[0] = tp_begin[org_val];
                    p_site_mixture[1] = tp_begin[4+org_val];
                    p_site_mixture[2] = tp_begin[8+org_val];
                    p_site_mixture[3] = tp_begin[12+org_val];
                        
                }
                
            } // end-if a gap state
            
            
            // increment the pointers to next site
            p_site_mixture+=this->site_offset;
            
        } // end-for over all sites/patterns in the sequence
        
        // increment the pointers to next mixture category
        p_mixture+=mixture_offset;
        
    } // end-for over all mixture categories
    
}


#endif
