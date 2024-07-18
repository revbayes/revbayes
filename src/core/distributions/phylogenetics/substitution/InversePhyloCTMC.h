#ifndef InversePhyloCTMC_H
#define InversePhyloCTMC_H

#include "AbstractHomologousDiscreteCharacterData.h"
#include "DnaState.h"
#include "RateMatrix.h"
#include "RbVector.h"
#include "TopologyNode.h"
#include "TransitionProbabilityMatrix.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
    
    /**
     * Provides the inverse of the probability of a phylogenetic character distribution.
     *
     * @copyright Copyright 2024-
     * @author Martin R. Smith
     * @since 2024-07-14, version 1.2.5
     *
     */

    template<typename valType>
    class InversePhyloCTMC : public TypedDistribution< valType > {

    public:
        InversePhyloCTMC(TypedDistribution< valType > base_dist);
        virtual                                            ~InversePhyloCTMC(void);                                                                   //!< Virtual destructor

        // public member functions
        InversePhyloCTMC*                                   clone(void) const;                                                                          //!< Create an independent clone


    protected:

        virtual void                                        computeRootLikelihood(size_t root, size_t l, size_t r);
        virtual void                                        computeRootLikelihood(size_t root, size_t l, size_t r, size_t m);
        virtual void                                        computeInternalNodeLikelihood(const TopologyNode &n, size_t nIdx, size_t l, size_t r);
        virtual void                                        computeInternalNodeLikelihood(const TopologyNode &n, size_t nIdx, size_t l, size_t r, size_t m);
        virtual void                                        computeTipLikelihood(const TopologyNode &node, size_t nIdx);


    private:
        RevPtr< const RevVariable > base_distribution;
    };

}

#endif