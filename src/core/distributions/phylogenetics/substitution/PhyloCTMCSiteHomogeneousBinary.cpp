#include "PhyloCTMCSiteHomogeneousBinary.h"

#include <cmath>
#include <cstddef>
#include <utility>

#include "RbBitSet.h"
#include "TopologyNode.h"
#include "Tree.h"
#include "TypedDagNode.h"

RevBayesCore::PhyloCTMCSiteHomogeneousBinary::PhyloCTMCSiteHomogeneousBinary(const TypedDagNode<Tree> *t, bool c, size_t nSites, bool amb, BinaryAscertainmentBias::Coding ty, bool internal, bool gapmatch) :
PhyloCTMCSiteHomogeneousConditional<BinaryState>(  t, 2, c, nSites, amb, AscertainmentBias::Coding(ty), internal, gapmatch )
{
    
}

RevBayesCore::PhyloCTMCSiteHomogeneousBinary* RevBayesCore::PhyloCTMCSiteHomogeneousBinary::clone( void ) const {
    
    return new PhyloCTMCSiteHomogeneousBinary( *this );
}


bool RevBayesCore::PhyloCTMCSiteHomogeneousBinary::isSitePatternCompatible( std::map<size_t, size_t> charCounts )
{
    std::map<size_t, size_t>::iterator zero = charCounts.find(0);
    std::map<size_t, size_t>::iterator one  = charCounts.find(1);
    
    bool compatible = true;
    
    if ( charCounts.size() == 1 )
    {
        if (zero != charCounts.end() && (coding & BinaryAscertainmentBias::NOABSENCESITES) )
        {
            compatible = false;
        }
        else if (one != charCounts.end()  && (coding & BinaryAscertainmentBias::NOPRESENCESITES) )
        {
            compatible = false;
        }
    }
    else
    {
        if (zero != charCounts.end() && zero->second == 1 && (coding & BinaryAscertainmentBias::NOSINGLETONABSENCE) )
        {
            compatible = false;
        }
        else if (one != charCounts.end() && one->second == 1 && (coding & BinaryAscertainmentBias::NOSINGLETONPRESENCE) )
        {
            compatible = false;
        }
    }
    
    return compatible;
}


bool RevBayesCore::PhyloCTMCSiteHomogeneousBinary::isSitePatternCompatible( std::map<RbBitSet, size_t> charCounts )
{
    RbBitSet z(2, false);
    RbBitSet o(2, false);
    
    z.set(0);
    o.set(1);
    
    std::map<RbBitSet, size_t>::iterator zero = charCounts.find(z);
    std::map<RbBitSet, size_t>::iterator one  = charCounts.find(o);
    
    bool compatible = true;
    
    if ( charCounts.size() == 1 )
    {
        if (zero != charCounts.end() && (coding & BinaryAscertainmentBias::NOABSENCESITES) )
        {
            compatible = false;
        }
        else if (one != charCounts.end()  && (coding & BinaryAscertainmentBias::NOPRESENCESITES) )
        {
            compatible = false;
        }
    }
    else
    {
        if (zero != charCounts.end() && zero->second == 1 && (coding & BinaryAscertainmentBias::NOSINGLETONABSENCE) )
        {
            compatible = false;
        }
        else if (one != charCounts.end() && one->second == 1 && (coding & BinaryAscertainmentBias::NOSINGLETONPRESENCE) )
        {
            compatible = false;
        }
    }
    
    return compatible;
}

double RevBayesCore::PhyloCTMCSiteHomogeneousBinary::sumRootLikelihood( void )
{
    double sumPartialProbs = PhyloCTMCSiteHomogeneous<BinaryState>::sumRootLikelihood();
    
    if (coding == BinaryAscertainmentBias::ALL)
        return sumPartialProbs;
    
    // get the root node
    const TopologyNode &root = this->tau->getValue().getRoot();
    
    // get root frequencies
    std::vector<std::vector<double> >   ff;
    getRootFrequencies(ff);
    
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
            const std::vector<double> &f = ff[mixture % ff.size()];
            
            double prob = 0.0;
            
            // iterate over ancestral (non-autapomorphic) states
            for (size_t a = 0; a < this->num_states; a++)
            {
                size_t offset = mixture*correctionmixture_offset + mask*correctionMaskOffset + a*correctionOffset;
                
                std::vector<double>::const_iterator             u = p_node   + offset;
                
                // iterate over combinations of autapomorphic states
                for (size_t c = 0; c < numCorrectionPatterns; c++)
                {
                    // constant site pattern likelihoods
                    std::vector<double>::const_iterator         uc = u  + c*this->num_states;
                    
                    if ( ((coding & BinaryAscertainmentBias::NOABSENCESITES)      && a == 0 && c == 0) ||
                       ((coding & BinaryAscertainmentBias::NOPRESENCESITES)     && a == 1 && c == 0) ||
                       ((coding & BinaryAscertainmentBias::NOSINGLETONPRESENCE) && a == 0 && c == 1
                        && (maskObservationCounts[mask] > 1 || !(coding & BinaryAscertainmentBias::NOPRESENCESITES)) ) ||
                       ((coding & BinaryAscertainmentBias::NOSINGLETONABSENCE)  && a == 1 && c == 1
                        && (maskObservationCounts[mask] > 1 || !(coding & BinaryAscertainmentBias::NOABSENCESITES))
                        && (maskObservationCounts[mask] > 2 || !(coding & BinaryAscertainmentBias::NOSINGLETONPRESENCE)) )
                       )
                    {
                        // iterate over initial states
                        for (size_t ci = 0; ci < this->num_states; ci++)
                        {
                            prob += uc[ci];
                        }
                    }
                }
            }
            
            // impose a per-mixture boundary
            if (prob <= 0.0 || prob >= 1.0)
            {
                prob = RbConstants::Double::nan;
            }
            
            perMaskCorrections[mask] += prob * mixtureProbs[mixture];
            perMaskMixtureCorrections[mask*this->num_site_mixtures + mixture] = (1.0 - prob) * mixtureProbs[mixture];
            
            // add corrections for invariant sites
            double prob_invariant = getPInv();
            if (prob_invariant > 0.0)
            {
                perMaskCorrections[mask] *= (1.0 - prob_invariant);
                perMaskMixtureCorrections[mask*this->num_site_mixtures + mixture] *= (1.0 - prob_invariant);
                
                if (coding & BinaryAscertainmentBias::NOABSENCESITES)
                    perMaskCorrections[mask] += prob_invariant * f[0] * mixtureProbs[mixture];
                
                if (coding & BinaryAscertainmentBias::NOPRESENCESITES)
                    perMaskCorrections[mask] += prob_invariant * f[1] * mixtureProbs[mixture];
            }
            
        }
        
        // add corrections for invariant sites
//        double prob_invariant = getPInv();
//        if (prob_invariant > 0.0)
//        {
//            double mean = 0.0;
//            for (size_t i = 0; i < ff.size(); i++)
//            {
//                if (coding & BinaryAscertainmentBias::NOABSENCESITES)
//                    mean += ff[i][0];
//                
//                if (coding & BinaryAscertainmentBias::NOPRESENCESITES)
//                    mean += ff[i][1];
//            }
//            
//            mean /= ff.size();
//            
//            perMaskCorrections[mask] += mean * prob_invariant;
//        }
        
        // normalize the log-probability
//        perMaskCorrections[mask] /= this->num_site_mixtures;
        
        // impose a per-mask boundary
        if (perMaskCorrections[mask] <= 0.0 || perMaskCorrections[mask] >= 1.0)
            perMaskCorrections[mask] = RbConstants::Double::nan;
        
        perMaskCorrections[mask] = log(1.0 - perMaskCorrections[mask]);
        
        // apply the correction for this correction mask
        sumPartialProbs -= perMaskCorrections[mask]*correctionMaskCounts[mask];
    }
    
    return sumPartialProbs;
}
