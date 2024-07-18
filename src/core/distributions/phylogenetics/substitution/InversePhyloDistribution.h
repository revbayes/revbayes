#ifndef InversePhyloDistribution_H
#define InversePhyloDistribution_H

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

    class InversePhyloDistribution : public TypedDistribution< AbstractHomologousDiscreteCharacterData > {

    public:
        InversePhyloDistribution(const TypedDistribution<AbstractHomologousDiscreteCharacterData>& d) 
            : TypedDistribution<AbstractHomologousDiscreteCharacterData>(), base_distribution(d) {}
        
        virtual                                            ~InversePhyloDistribution(void);                                                                   //!< Virtual destructor

        // public member functions
        InversePhyloDistribution*                           clone(void) const;                                                                          //!< Create an independent clone
        
        double computeLnProbability(void) {
            return -base_dist.computeLnProbability();
        }

    private:
        RevPtr< const RevVariable > base_distribution;
        // or?
        // const TypedDistribution<AbstractHomologousDiscreteCharacterData>& base_distribution;
    };


}

#endif