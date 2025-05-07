#ifndef PhyloCTMCSiteHomogeneousConditional_H
#define PhyloCTMCSiteHomogeneousConditional_H

#include <bitset>
#include <cmath>
#include "PhyloCTMCSiteHomogeneous.h"
#include "DistributionBinomial.h"
#include "DistributionNegativeBinomial.h"
#include "RlUserInterface.h"

namespace RevBayesCore {

    struct AscertainmentBias {
        
        enum Coding { ALL                 = 0x00,
                      VARIABLE            = 0x03,
                      INFORMATIVE         = 0x0F
                    };
    };

    template<class charType>
    class PhyloCTMCSiteHomogeneousConditional : public PhyloCTMCSiteHomogeneous<charType> {

    public:
        PhyloCTMCSiteHomogeneousConditional(const TypedDagNode< Tree > *t, size_t nChars, bool c, size_t nSites, bool amb, AscertainmentBias::Coding cod = AscertainmentBias::ALL, bool internal = false, bool gapmatch = true);
        PhyloCTMCSiteHomogeneousConditional(const PhyloCTMCSiteHomogeneousConditional &n);
        virtual                                            ~PhyloCTMCSiteHomogeneousConditional(void);                                                                   //!< Virtual destructor

        // public member functions
        PhyloCTMCSiteHomogeneousConditional*                clone(void) const;                                                                        //!< Create an independent clone
        void                                                setValue(AbstractHomologousDiscreteCharacterData *v, bool f=false);
        virtual void                                        redrawValue(void);

    protected:

        virtual void                                        computeRootLikelihood(size_t root, size_t l, size_t r);
        virtual void                                        computeRootLikelihood(size_t root, size_t l, size_t r, size_t m);
        virtual void                                        computeInternalNodeLikelihoodBranchWise(const TopologyNode &n, size_t nIdx, size_t l, size_t r);
        virtual void                                        computeTipLikelihood(const TopologyNode &node, size_t nIdx);

        virtual void                                        computeRootCorrection(size_t root, size_t l, size_t r);
        virtual void                                        computeRootCorrection(size_t root, size_t l, size_t r, size_t m);
        virtual void                                        computeInternalNodeCorrection(const TopologyNode &n, size_t nIdx, size_t l, size_t r);
        virtual void                                        computeTipCorrection(const TopologyNode &node, size_t nIdx);

        virtual void                                        resizeLikelihoodVectors(void);

        bool                                                warned;

        int                                                 coding;
        size_t                                              N;

        size_t                                              numCorrectionMasks;
        size_t                                              activeCorrectionOffset;
        size_t                                              correctionnode_offset;
        size_t                                              correctionmixture_offset;
        size_t                                              correctionMaskOffset;
        size_t                                              correctionOffset;
        size_t                                              numCorrectionPatterns;

        std::vector<double>                                 correctionLikelihoods;

        std::vector<std::vector<bool> >                     correctionMaskMatrix;
        std::vector<size_t>                                 correctionMaskCounts;
        std::vector<size_t>                                 maskObservationCounts;
        std::vector<double>                                 perMaskMixtureCorrections;

        virtual double                                      sumRootLikelihood( void );
        virtual bool                                        isSitePatternCompatible( std::map<size_t, size_t> );
        virtual bool                                        isSitePatternCompatible( std::map<RbBitSet, size_t> );
        std::vector<size_t>                                 getIncludedSiteIndices( void );
        void                                                updateCorrections( const TopologyNode& node, size_t nodeIndex );

    private:
        virtual void                                        simulateConditional( const TopologyNode &node, std::vector<charType> &taxa, size_t rateIndex, std::map<size_t, size_t>& charCounts);
    };

}

#include <climits>

template<class charType>
RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::PhyloCTMCSiteHomogeneousConditional(const TypedDagNode<Tree> *t, size_t nChars, bool c, size_t nSites, bool amb, AscertainmentBias::Coding ty, bool internal, bool gapmatch) :
    PhyloCTMCSiteHomogeneous<charType>(  t, nChars, c, nSites, amb, internal, gapmatch ), warned(false), coding(ty), N(nSites), numCorrectionMasks(0)
{
    if (coding != AscertainmentBias::ALL)
    {
        numCorrectionMasks      = 1;
        correctionMaskMatrix    = std::vector<std::vector<bool> >(1, std::vector<bool>(this->num_nodes,0) );
        
        // number of correction patterns per character state
        if (coding & (AscertainmentBias::INFORMATIVE ^ AscertainmentBias::VARIABLE))
        {
            numCorrectionPatterns = std::pow(2.0f, float(this->num_states - 1));
        }
        else
        {
            numCorrectionPatterns = 1;
        }

        correctionOffset        = this->num_states*numCorrectionPatterns;
        correctionMaskOffset    = correctionOffset*this->num_states;
        correctionmixture_offset = numCorrectionMasks*correctionMaskOffset;
        correctionnode_offset    = this->num_site_mixtures*correctionmixture_offset;
        activeCorrectionOffset  = this->num_nodes*correctionnode_offset;

        correctionLikelihoods = std::vector<double>(activeCorrectionOffset*2, 0.0);

        perMaskMixtureCorrections = std::vector<double>(this->num_site_mixtures*numCorrectionMasks, 0.0);
    }
}


template<class charType>
RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::~PhyloCTMCSiteHomogeneousConditional( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!

}


template<class charType>
RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>* RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::clone( void ) const {

    return new PhyloCTMCSiteHomogeneousConditional<charType>( *this );
}

template<class charType>
RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::PhyloCTMCSiteHomogeneousConditional(const PhyloCTMCSiteHomogeneousConditional &n) :
    PhyloCTMCSiteHomogeneous< charType>( n ),
    warned(n.warned),
    coding(n.coding),
    N(n.N),
    numCorrectionMasks(n.numCorrectionMasks),
    activeCorrectionOffset(n.activeCorrectionOffset),
    correctionnode_offset(n.activeCorrectionOffset),
    correctionmixture_offset(n.activeCorrectionOffset),
    correctionMaskOffset(n.correctionMaskOffset),
    correctionOffset(n.correctionOffset),
    numCorrectionPatterns(n.numCorrectionPatterns),
    correctionLikelihoods(n.correctionLikelihoods),
    correctionMaskMatrix(n.correctionMaskMatrix),
    correctionMaskCounts(n.correctionMaskCounts),
    maskObservationCounts(n.maskObservationCounts),
    perMaskMixtureCorrections(n.perMaskMixtureCorrections)
{

}

template<class charType>
std::vector<size_t> RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::getIncludedSiteIndices( void )
{
    correctionMaskCounts.clear();
    maskObservationCounts.clear();
    correctionMaskMatrix.clear();

    // find the unique site patterns and compute their respective frequencies
    std::vector<TopologyNode*> nodes = this->tau->getValue().getNodes();
    size_t tips = this->tau->getValue().getNumberOfTips();

    // create a vector with the correct site indices
    // some of the sites may have been excluded
    std::vector<size_t> siteIndices;
    size_t siteIndex = 0;
    size_t incompatible = 0;

    std::map<std::string, size_t> maskIndices;

    if (coding != AscertainmentBias::ALL)
    {
        // if we are using a correction, add a correction mask with 0 gaps.
        // it is required when simulating data, but not necessarily used
        // in computing the likelihood (e.g. if all sites have at least one gap)

        std::string gapless(tips, ' ');
        std::vector<bool> mask(tips, false);

        maskIndices[gapless] = 0;

        correctionMaskCounts.push_back(0);
        maskObservationCounts.push_back(tips);
        correctionMaskMatrix.push_back(mask);
    }

    this->num_sites = this->value->getNumberOfIncludedCharacters();

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


        std::map<size_t,size_t> charCounts;
        std::map<RbBitSet,size_t> bitCounts;
        size_t numGap = 0;

        std::string mask = "";
        std::vector<bool> maskData;
        for (std::vector<TopologyNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it)
        {
            if ( (*it)->isTip() )
            {
                AbstractDiscreteTaxonData& taxon = this->value->getTaxonData( (*it)->getName() );
                DiscreteCharacterState &c = taxon.getCharacter(siteIndex);

                bool gap = c.isGapState();
                // if we treat unknown characters as gaps and this is an unknown character then we change it
                // because we might then have a pattern more
                if ( this->treatAmbiguousAsGaps && (c.isAmbiguous() || c.isMissingState()) )
                {
                    gap = true;
                }
                else if ( this->treatUnknownAsGap && (c.getNumberOfStates() == c.getNumberObservedStates() || c.isMissingState()) )
                {
                    gap = true;
                }

                if (gap)
                {
                    if (coding != AscertainmentBias::ALL)
                        mask += "-";
                    numGap++;
                }
                else
                {
                    if (coding != AscertainmentBias::ALL)
                        mask += " ";

                    bitCounts[c.getState()]++;
                }

                if (coding != AscertainmentBias::ALL)
                    maskData.push_back(gap);
            }
        }

        if ( !isSitePatternCompatible(bitCounts) )
        {
            incompatible++;
        }
        else
        {
            siteIndices.push_back(siteIndex);

            if (coding != AscertainmentBias::ALL)
            {
                // increase the count for this mask
                std::map<std::string, size_t>::iterator it = maskIndices.find(mask);
                if (it != maskIndices.end())
                {
                    correctionMaskCounts[it->second]++;
                }
                else
                {
                    maskIndices[mask] = correctionMaskCounts.size();

                    correctionMaskCounts.push_back(1);
                    maskObservationCounts.push_back(tips - numGap);
                    correctionMaskMatrix.push_back(maskData);
                }
            }
        }

        siteIndex++;
    }

    // test if there were additional sites that we did not use
    while ( siteIndex < this->value->getNumberOfCharacters() )
    {
        if ( !this->value->isCharacterExcluded(siteIndex)  )
        {
            throw RbException( "The character matrix cannot set to this variable because it has too many included characters." );
        }
        siteIndex++;
    }

    // warn if we have to exclude some incompatible characters.
    if (incompatible > 0)
    {
        if (!warned)
        {
            std::stringstream ss;
            ss << "WARNING: There are " << incompatible << " characters incompatible with the specified coding bias. These characters will be excluded.";
            RBOUT(ss.str());

            warned = true;
        }

        // resize our datset to account for the newly excluded characters
        this->num_sites = siteIndices.size();
        this->site_pattern.clear();
        this->site_pattern.resize(this->num_sites);
    }

    // readjust the number of correction sites to account for masked sites
    if (coding != AscertainmentBias::ALL)
    {
        numCorrectionMasks = correctionMaskCounts.size();
    }

    return siteIndices;
}


template<class charType>
bool RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::isSitePatternCompatible( std::map<size_t, size_t> charCounts )
{

    size_t num_informative = 0;
    for (std::map<size_t, size_t>::iterator it = charCounts.begin(); it != charCounts.end(); it++)
    {
        if (it->second > 1)
        {
            num_informative++;
        }
    }

    bool informative = num_informative > 1;

    if ( (charCounts.size() == 1 && (coding == AscertainmentBias::VARIABLE)) || (!informative && (coding == AscertainmentBias::INFORMATIVE)) )
    {
        return false;
    }

    return true;
}


// Thanks to MH
template<class charType>
bool RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::isSitePatternCompatible( std::map<RbBitSet, size_t> charCounts )
{
    // if the pattern is constant, then it is incompatible with the variable coding
    if (charCounts.size() == 1) return !(coding & AscertainmentBias::VARIABLE);

    // find the common_states
    std::map<size_t,size_t> stateCounts;
    size_t max = 0;

    for (std::map<RbBitSet, size_t>::iterator it = charCounts.begin(); it != charCounts.end(); it++)
    {
        RbBitSet r = it->first;
        for (size_t i = 0; i < r.size(); i++)
        {
            stateCounts[i] += r.test(i) * it->second;
            if (stateCounts[i] > max)
            {
                max = stateCounts[i];
            }
        }
    }

    std::vector<size_t> common_states;
    for (std::map<size_t, size_t>::iterator it = stateCounts.begin(); it != stateCounts.end(); it++)
    {
        if (it->second == max)
        {
            common_states.push_back(it->first);
        }
    }

    // find characters not intersecting common state
    // then get state counts
    stateCounts.clear();
    for (size_t i = 0; i < common_states.size(); i++)
    {
        for (std::map<RbBitSet, size_t>::iterator it = charCounts.begin(); it != charCounts.end(); it++)
        {
            RbBitSet r = it->first;

            if ( r.test(common_states[i]) ) continue;

            if ( it->second > 1 )
            {
                return true;
            }

            for (size_t i = 0; i < r.size(); i++)
            {
                // if a state is found more than once among characters lacking the common state
                // then this site pattern is parsimony informative
                if (stateCounts.find(i) != stateCounts.end() )
                {
                    return true;
                }
                else
                {
                    stateCounts[i] += r.test(i);
                }
            }
        }
    }

    // if the pattern is uninformative, then it is incompatible with the informative coding
    return (coding != AscertainmentBias::INFORMATIVE);
}


template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::computeRootLikelihood( size_t root, size_t left, size_t right)
{

    PhyloCTMCSiteHomogeneous<charType>::computeRootLikelihood(root, left, right);

    if (coding != AscertainmentBias::ALL)
    {
        computeRootCorrection(root, left, right);
    }
}


template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::computeRootLikelihood( size_t root, size_t left, size_t right, size_t middle)
{

    PhyloCTMCSiteHomogeneous<charType>::computeRootLikelihood(root, left, right, middle);

    if (coding != AscertainmentBias::ALL)
    {
        computeRootCorrection(root, left, right, middle);
    }
}


template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::computeInternalNodeLikelihoodBranchWise(const TopologyNode &node, size_t node_index, size_t left, size_t right)
{

    PhyloCTMCSiteHomogeneous<charType>::computeInternalNodeLikelihoodBranchWise(node, node_index, left, right);

    if (coding != AscertainmentBias::ALL)
    {
        computeInternalNodeCorrection(node, node_index, left, right);
    }
}




template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::computeTipLikelihood(const TopologyNode &node, size_t node_index)
{

    PhyloCTMCSiteHomogeneous<charType>::computeTipLikelihood(node, node_index);

    if (coding != AscertainmentBias::ALL)
    {
        computeTipCorrection(node, node_index);
    }
}

template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::computeTipCorrection(const TopologyNode &node, size_t node_index)
{
    std::vector<double>::iterator p_node = correctionLikelihoods.begin() + this->active_branch_likelihood[node_index]*this->activeCorrectionOffset + node_index*correctionnode_offset;
    
    size_t data_tip_index = this->taxon_name_2_tip_index_map[ node.getName() ];
    
    size_t pmat_offset = this->active_pmatrices[node_index] * this->active_P_matrix_offset + node_index * this->pmat_node_offset;

    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
    {
        // const TransitionProbabilityMatrix&    pij = this->transition_prob_matrices[mixture];
        const TransitionProbabilityMatrix&    pij = this->pmatrices[pmat_offset + mixture];

        // iterate over correction masks
        for (size_t mask = 0; mask < numCorrectionMasks; mask++)
        {
            bool gap = correctionMaskMatrix[mask][data_tip_index];

            // iterate over ancestral (non-autapomorphic) states
            for (size_t a = 0; a < this->num_states; a++)
            {
                size_t offset = mixture*correctionmixture_offset + mask*correctionMaskOffset + a*correctionOffset;

                std::vector<double>::iterator                 u = p_node   + offset;

                // iterate over constant/autapomorphic states
                for (size_t i = 0; i < this->num_states; i++)
                {
                    if (i > 0 && numCorrectionPatterns == 1)
                        break;

                    // get bit pattern for this constant/autapomorphic state
                    size_t c = ( i == 0 ? 0 : std::pow(2.0f, float(i - 1)) );

                    std::vector<double>::iterator         uc = u  + c*this->num_states;

                    // if we have a gap, then fill with ones for constant patterns
                    // or zero for other patterns
                    if (gap)
                    {
                        std::fill(uc, uc + this->num_states, i == 0);
                    }
                    // if this is a constant pattern, fill with transition prob from ci to a
                    else if (c == 0)
                    {
                        // iterate over initial states
                        for (size_t ci = 0; ci < this->num_states; ci++)
                        {
                            uc[ci] = pij[ci][a];
                        }
                    }
                    // if this is an autapomorphic state cj, fill with transition prob from ci to cj
                    else
                    {
                        // get autamorphic state
                        size_t cj = (i-1 >= a ? i : i-1);

                        // iterate over initial states
                        for (size_t ci = 0; ci < this->num_states; ci++)
                        {
                            uc[ci] = pij[ci][cj];
                        }
                    }
                }
            }
        }
    }
}


template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::computeInternalNodeCorrection(const TopologyNode &node, size_t node_index, size_t left, size_t right)
{
    // get the pointers to the partial likelihoods for this node and the two descendant subtrees
    std::vector<double>::const_iterator   p_left  = correctionLikelihoods.begin() + this->active_branch_likelihood[left]*activeCorrectionOffset + left*correctionnode_offset;
    std::vector<double>::const_iterator   p_right = correctionLikelihoods.begin() + this->active_branch_likelihood[right]*activeCorrectionOffset + right*correctionnode_offset;
    std::vector<double>::iterator         p_node  = correctionLikelihoods.begin() + this->active_branch_likelihood[node_index]*activeCorrectionOffset + node_index*correctionnode_offset;
    
    size_t pmat_offset = this->active_pmatrices[node_index] * this->active_P_matrix_offset + node_index * this->pmat_node_offset;

    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
    {
        // const TransitionProbabilityMatrix&    pij = this->transition_prob_matrices[mixture];
        const TransitionProbabilityMatrix&    pij = this->pmatrices[pmat_offset + mixture];

        // iterate over correction masks
        for (size_t mask = 0; mask < numCorrectionMasks; mask++)
        {
            // iterate over ancestral (non-autapomorphic) states
            for (size_t a = 0; a < this->num_states; a++)
            {
                size_t offset = mixture*correctionmixture_offset + mask*correctionMaskOffset + a*correctionOffset;

                std::vector<double>::iterator                 u = p_node   + offset;
                std::vector<double>::const_iterator         u_l = p_left   + offset;
                std::vector<double>::const_iterator         u_r = p_right  + offset;

                // iterate over combinations of autapomorphic states
                for (size_t c = 0; c < numCorrectionPatterns; c++)
                {
                    std::vector<double>::iterator         uc = u  + c*this->num_states;

                    std::fill(uc, uc + this->num_states, 0.0);

                    // iterate over partitions of c
                    for (size_t p1 = 0; p1 <= c; p1++)
                    {
                        if ( (p1 | c) == c)
                        {
                            size_t p2 = p1 ^ c;

                            std::vector<double>::const_iterator         lc = u_l  + p1*this->num_states;
                            std::vector<double>::const_iterator         rc = u_r  + p2*this->num_states;

                            // iterate over initial states
                            for (size_t ci = 0; ci < this->num_states; ci++)
                            {
                                // iterate over ending states
                                for (size_t cj = 0; cj < this->num_states; cj++)
                                {
                                    uc[ci] += pij[ci][cj] * lc[cj] * rc[cj];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::computeRootCorrection( size_t root, size_t left, size_t right, size_t middle)
{
    // get the pointers to the partial likelihoods for this node and the two descendant subtrees
    std::vector<double>::iterator         p_node   = correctionLikelihoods.begin() + this->active_branch_likelihood[root]*activeCorrectionOffset + root*correctionnode_offset;
    std::vector<double>::const_iterator   p_left   = correctionLikelihoods.begin() + this->active_branch_likelihood[left]*activeCorrectionOffset + left*correctionnode_offset;
    std::vector<double>::const_iterator   p_right  = correctionLikelihoods.begin() + this->active_branch_likelihood[right]*activeCorrectionOffset + right*correctionnode_offset;
    std::vector<double>::const_iterator   p_middle = correctionLikelihoods.begin() + this->active_branch_likelihood[middle]*activeCorrectionOffset + middle*correctionnode_offset;

    // get the root frequencies
    std::vector<std::vector<double> > ff;
    this->getRootFrequencies(ff);

    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
    {
        const std::vector<double> &f = ff[mixture % ff.size()];

        // iterate over correction masks
        for (size_t mask = 0; mask < numCorrectionMasks; mask++)
        {
            // iterate over ancestral (non-autapomorphic) states
            for (size_t a = 0; a < this->num_states; a++)
            {
                size_t offset = mixture*correctionmixture_offset + mask*correctionMaskOffset + a*correctionOffset;

                std::vector<double>::iterator                 u = p_node   + offset;
                std::vector<double>::const_iterator         u_l = p_left   + offset;
                std::vector<double>::const_iterator         u_r = p_right  + offset;
                std::vector<double>::const_iterator         u_m = p_middle + offset;

                // iterate over combinations of autapomorphic states
                for (size_t c = 0; c < numCorrectionPatterns; c++)
                {
                    std::vector<double>::iterator         uc = u  + c*this->num_states;

                    std::fill(uc, uc + this->num_states, 0.0);

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

                                    std::vector<double>::const_iterator         lc = u_l  + p1*this->num_states;
                                    std::vector<double>::const_iterator         rc = u_r  + p2*this->num_states;
                                    std::vector<double>::const_iterator         mc = u_m  + p3*this->num_states;

                                    // iterate over initial states
                                    for (size_t ci = 0; ci < this->num_states; ci++)
                                    {
                                        uc[ci] += f[ci] * lc[ci] * rc[ci] * mc[ci];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::computeRootCorrection( size_t root, size_t left, size_t right)
{
    // get the pointers to the partial likelihoods for this node and the two descendant subtrees
    std::vector<double>::iterator         p_node  = correctionLikelihoods.begin() + this->active_branch_likelihood[root]*activeCorrectionOffset + root*correctionnode_offset;
    std::vector<double>::const_iterator   p_left  = correctionLikelihoods.begin() + this->active_branch_likelihood[left]*activeCorrectionOffset + left*correctionnode_offset;
    std::vector<double>::const_iterator   p_right = correctionLikelihoods.begin() + this->active_branch_likelihood[right]*activeCorrectionOffset + right*correctionnode_offset;

    // get the root frequencies
    std::vector<std::vector<double> > ff;
    this->getRootFrequencies(ff);

    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
    {
        const std::vector<double> &f = ff[mixture % ff.size()];

        // iterate over correction masks
        for (size_t mask = 0; mask < numCorrectionMasks; mask++)
        {
            // iterate over ancestral (non-autapomorphic) states
            for (size_t a = 0; a < this->num_states; a++)
            {
                size_t offset = mixture*correctionmixture_offset + mask*correctionMaskOffset + a*correctionOffset;

                std::vector<double>::iterator                 u = p_node   + offset;
                std::vector<double>::const_iterator         u_l = p_left   + offset;
                std::vector<double>::const_iterator         u_r = p_right  + offset;

                // iterate over combinations of autapomorphic states
                for (size_t c = 0; c < numCorrectionPatterns; c++)
                {
                    std::vector<double>::iterator         uc = u  + c*this->num_states;

                    std::fill(uc, uc + this->num_states, 0.0);
                    
                    // iterate over partitions of c
                    for (size_t p1 = 0; p1 <= c; p1++)
                    {
                        if ( (p1 | c) == c)
                        {
                            size_t p2 = p1 ^ c;

                            std::vector<double>::const_iterator         lc = u_l  + p1*this->num_states;
                            std::vector<double>::const_iterator         rc = u_r  + p2*this->num_states;

                            // iterate over initial states
                            for (size_t ci = 0; ci < this->num_states; ci++)
                            {
                                uc[ci] += f[ci] * lc[ci] * rc[ci];
                            }
                        }
                    }
                }
            }
        }
    }
}

template<class charType>
double RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::sumRootLikelihood( void )
{
    double sumPartialProbs = PhyloCTMCSiteHomogeneous<charType>::sumRootLikelihood();
    
    if (coding == AscertainmentBias::ALL)
        return sumPartialProbs;

    // get the root node
    const TopologyNode &root = this->tau->getValue().getRoot();

    // get the index of the root node
    size_t node_index = root.getIndex();
    
    std::vector<double>::const_iterator p_node = correctionLikelihoods.begin() + this->active_branch_likelihood[node_index] * activeCorrectionOffset  + node_index*correctionnode_offset;
    
    std::vector<double> perMaskCorrections = std::vector<double>(numCorrectionMasks, 0.0);
    
    std::vector<double> mixtureProbs = this->getMixtureProbs();
    
    // iterate over each correction mask
    for (size_t mask = 0; mask < numCorrectionMasks; mask++)
    {
        // iterate over all mixture categories
        for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
        {
            double prob = 0.0;
            
            // iterate over ancestral (non-autapomorphic) states
            for (size_t a = 0; a < this->num_states; a++)
            {
                size_t offset = mixture*correctionmixture_offset + mask*correctionMaskOffset + a*correctionOffset;

                std::vector<double>::const_iterator             u = p_node   + offset;

                // iterate over combinations of autapomorphic states
                for (size_t c = 0; c < numCorrectionPatterns; c++)
                {
                    std::vector<double>::const_iterator         uc = u  + c*this->num_states;

                    if (c == 0 && coding == AscertainmentBias::VARIABLE)
                    {
                        // iterate over initial states
                        for (size_t ci = 0; ci < this->num_states; ci++)
                        {
                            prob += uc[ci];
                        }
                    }
                    else if (coding == AscertainmentBias::INFORMATIVE)
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
                        if ( (maskObservationCounts[mask] == popcount + 1 && a == 0) ||
                             maskObservationCounts[mask] > popcount + 1)
                        {
                            // iterate over initial states
                            for (size_t ci = 0; ci < this->num_states; ci++)
                            {
                                prob += uc[ci];
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

//            perMaskCorrections[mask] += prob * 0.25;
//            perMaskMixtureCorrections[mask*this->num_site_mixtures + mixture] = (1.0 - prob) * 0.25;
            perMaskCorrections[mask] += prob * mixtureProbs[mixture];
            perMaskMixtureCorrections[mask*this->num_site_mixtures + mixture] = (1.0 - prob) * mixtureProbs[mixture];
        
            // add corrections for invariant sites
            double prob_invariant = this->getPInv();
            if (prob_invariant > 0.0)
            {
                perMaskCorrections[mask] *= (1.0 - prob_invariant);
                perMaskMixtureCorrections[mask*this->num_site_mixtures + mixture] *= (1.0 - prob_invariant);
            }
            
        }
        
        // add corrections for invariant sites
        double prob_invariant = this->getPInv();
        if (prob_invariant > 0.0)
        {
            if ( coding & AscertainmentBias::VARIABLE)
            {
                perMaskCorrections[mask] += prob_invariant;
            }
        }

        // normalize the log-probability
//        perMaskCorrections[mask] /= this->num_site_mixtures;

        // impose a per-mask boundary
        if (perMaskCorrections[mask] < 0.0 || perMaskCorrections[mask] >= 1.0)
        {
            perMaskCorrections[mask] = RbConstants::Double::nan;
        }

        perMaskCorrections[mask] = log(1.0 - perMaskCorrections[mask]);
        
        // apply the correction for this correction mask
        sumPartialProbs -= perMaskCorrections[mask]*correctionMaskCounts[mask];
    }

    return sumPartialProbs;
}

template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::resizeLikelihoodVectors( void )
{

    RevBayesCore::PhyloCTMCSiteHomogeneous<charType>::resizeLikelihoodVectors();
    
    if (coding != AscertainmentBias::ALL)
    {
        correctionmixture_offset = numCorrectionMasks*correctionMaskOffset;
        correctionnode_offset    = this->num_site_mixtures*correctionmixture_offset;
        activeCorrectionOffset  = this->num_nodes*correctionnode_offset;

        correctionLikelihoods = std::vector<double>(activeCorrectionOffset*2, 0.0);

        perMaskMixtureCorrections = std::vector<double>(this->num_site_mixtures*numCorrectionMasks, 0.0);
    }
}

template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::updateCorrections( const TopologyNode& node, size_t nodeIndex )
{

    // only necessary if the node is actually dirty
    //if ( dirtyNodes[nodeIndex] )
    //{
    //    dirtyNodes[nodeIndex] = false;

        // start by filling the likelihood vector for the children of the root
        if ( node.isTip() )
        {
            // this->updateTransitionProbabilities(nodeIndex );
            computeTipCorrection( node, nodeIndex );
        }
        else if ( node.getNumberOfChildren() == 2 ) // rooted trees have two children for the root
        {
            const TopologyNode &left = node.getChild(0);
            size_t leftIndex = left.getIndex();
            updateCorrections( left, leftIndex );
            const TopologyNode &right = node.getChild(1);
            size_t rightIndex = right.getIndex();
            updateCorrections( right, rightIndex );

            if (node.isRoot())
                computeRootCorrection( nodeIndex, leftIndex, rightIndex );
            else
            {
                // this->updateTransitionProbabilities( nodeIndex );
                computeInternalNodeCorrection( node, nodeIndex, leftIndex, rightIndex );
            }

        }
        else if ( node.getNumberOfChildren() == 3 ) // unrooted trees have three children for the root
        {
            const TopologyNode &left = node.getChild(0);
            size_t leftIndex = left.getIndex();
            updateCorrections( left, leftIndex );
            const TopologyNode &right = node.getChild(1);
            size_t rightIndex = right.getIndex();
            updateCorrections( right, rightIndex );
            const TopologyNode &middle = node.getChild(2);
            size_t middleIndex = middle.getIndex();
            updateCorrections( middle, middleIndex );

            if (node.isRoot())
                computeRootCorrection( nodeIndex, leftIndex, rightIndex, middleIndex );
            else
            {
                // this->updateTransitionProbabilities( nodeIndex );
//                computeInternalNodeCorrection( node, nodeIndex, leftIndex, rightIndex, middleIndex );
            }

        }
        else
        {
            throw RbException("Node has an unexpected number of children. Only 2 or 3 are allowed.");
        }
    //}
}

template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::redrawValue( void ) {

    if (coding == AscertainmentBias::ALL)
    {
        PhyloCTMCSiteHomogeneous<charType>::redrawValue();
        
        return;
    }

    // delete the old value first
    delete this->value;

    // create a new character data object
    this->value = new HomologousDiscreteCharacterData<charType>();

    size_t num_tips = this->tau->getValue().getNumberOfTips();
    size_t num_nodes = this->tau->getValue().getNumberOfNodes();

    RandomNumberGenerator* rng = GLOBAL_RNG;

    const TopologyNode &root = this->tau->getValue().getRoot();
    size_t root_index = this->tau->getValue().getRoot().getIndex();

    // for safety, we resize the likelihood vectors here
    resizeLikelihoodVectors();
    
    updateCorrections(root, root_index);

    std::vector< DiscreteTaxonData<charType> > taxa = std::vector< DiscreteTaxonData<charType> >(num_tips, DiscreteTaxonData<charType>( Taxon("") ) );

    /*// first sample a total number of characters (M) from the marginal posterior:
    // M - N | N ~ NegBinomial(N+1, exp(lnCorrection) )
    double M_minus_N = RbStatistics::NegativeBinomial::rv(N + 1, exp(perMaskCorrections[0]), *rng);

    // then sample the observed number of characters (num_sites) from the likelihood:
    // num_sites | M ~ Binomial(M, exp(lnCorrection) )
    this->num_sites = RbStatistics::Binomial::rv( M_minus_N + N, exp(perMaskCorrections[0]), *rng);*/

    // sample the rate categories in proportion to the total probability (correction) for each mixture.
    double total = 0.0;
    std::vector<double> mixtureProbs = this->getMixtureProbs();
    for ( size_t i = 0; i < this->num_site_mixtures; ++i )
        total += mixtureProbs[i];

    std::vector<size_t> perSiteRates;
    for ( size_t i = 0; i < this->num_sites; ++i )
    {
        // draw the state
        double u = rng->uniform01()*total;
        size_t rateIndex = 0;

        double tmp = 0.0;
        while(tmp < u)
        {
            tmp += mixtureProbs[rateIndex];
            if (tmp < u)
                rateIndex++;
        }

        perSiteRates.push_back( rateIndex );
    }
    
    // get the root frequencies
    std::vector<std::vector<double> > ff;
    this->getRootFrequencies(ff);

    // then sample site-patterns using rejection sampling,
    // rejecting those that match the unobservable ones.
    for ( size_t i = 0; i < this->num_sites; i++ )
    {
        size_t rateIndex = perSiteRates[i];

        const std::vector< double > &freqs = ff[perSiteRates[i] % ff.size()];

        std::vector<charType> siteData(num_nodes, charType(this->num_states));

        // create the character
        charType &c = siteData[root_index];
        c.setToFirstState();
        // draw the state
        double u = rng->uniform01();
        size_t stateIndex = 0;
        while ( true )
        {
            u -= freqs[stateIndex];
            ++stateIndex;

            if ( u > 0.0 && stateIndex < this->num_states)
            {
                ++c;
            }
            else
            {
                break;
            }

        }

        // recursively simulate the sequences
        std::map<size_t, size_t> charCounts;
        simulateConditional( root, siteData, rateIndex, charCounts);

        if ( !isSitePatternCompatible(charCounts) )
        {
            i--;
            continue;
        }

        // add the taxon data to the character data
        for (size_t t = 0; t < num_tips; ++t)
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

    for (std::vector<bool>::iterator it = this->dirty_branches.begin(); it != this->dirty_branches.end(); ++it)
    {
        (*it) = true;
    }

    // flip the active likelihood pointers
    for (size_t index = 0; index < this->changed_nodes.size(); ++index)
    {
        if ( this->changed_nodes[index] == false )
        {
            this->active_branch_likelihood[index] = (this->active_branch_likelihood[index] == 0 ? 1 : 0);
            this->changed_nodes[index] = true;
        }
    }
    
    for (std::vector<bool>::iterator it = this->pmat_dirty_branches.begin(); it != this->pmat_dirty_branches.end(); ++it)
    {
        (*it) = true;
    }
    
    // update transition probability matrices
    this->updateTransitionProbabilityMatrices();
    updateCorrections(root, root_index);

    for (std::vector<bool>::iterator it = this->pmat_dirty_branches.begin(); it != this->pmat_dirty_branches.end(); ++it)
    {
        (*it) = true;
    }
    
    for (size_t index = 0; index < this->pmat_changed_nodes.size(); ++index)
    {
        if ( this->pmat_changed_nodes[index] == false )
        {
            this->active_pmatrices[index] = (this->active_pmatrices[index] == 0 ? 1 : 0);
            this->pmat_changed_nodes[index] = true;
        }
    }

}

template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::simulateConditional( const TopologyNode &node, std::vector<charType> &data, size_t rateIndex, std::map<size_t, size_t>& charCounts) {

    // get the children of the node
    const std::vector<TopologyNode*>& children = node.getChildren();

    // get the sequence of this node
    size_t node_index = node.getIndex();
    charType &parentState = data[ node_index ];

    // simulate the sequence for each child
    RandomNumberGenerator* rng = GLOBAL_RNG;
    for (std::vector< TopologyNode* >::const_iterator it = children.begin(); it != children.end(); ++it)
    {
        const TopologyNode &child = *(*it);

        // update the transition probability matrix
        this->updateTransitionProbabilities( child.getIndex() );

        std::uint64_t cp = parentState.getStateIndex();

        double *freqs = this->transition_prob_matrices[ rateIndex ][ cp ];

        // create the character
        charType &c = data[ child.getIndex() ];
        c.setToFirstState();
        // draw the state
        double u = rng->uniform01();
        size_t stateIndex = 0;
        while ( true )
        {
            u -= *freqs;
            ++stateIndex;

            if ( u > 0.0 && stateIndex < this->num_states)
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
            simulateConditional( child, data, rateIndex, charCounts);
    }

}

template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::setValue(AbstractHomologousDiscreteCharacterData *v, bool force)
{
    warned = false;

    // delegate to the parent class
    PhyloCTMCSiteHomogeneous<charType>::setValue(v, force);
}

#endif
