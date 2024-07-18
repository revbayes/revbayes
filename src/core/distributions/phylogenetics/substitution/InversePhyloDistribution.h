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
        InversePhyloDistribution(TypedDistribution< AbstractHomologousDiscreteCharacterData >& d)
            : TypedDistribution< AbstractHomologousDiscreteCharacterData >(d), base_distribution(d) {}

        // Virtual destructor
        virtual ~InversePhyloDistribution(void) = default;

        // public member functions
        
        InversePhyloDistribution* clone(void) const override {
            return new InversePhyloDistribution(*this);
        }
        
        double computeLnProbability(void) override {
            return -base_distribution.computeLnProbability();
        }

        void redrawValue(void) override {
            redrawValue();
        }
        
    protected:
        // Parameter management functions
        void swapParameterInternal(const DagNode *oldP, const DagNode *newP) override
        {
            base_distribution.swapParameter( oldP, newP );
        }


    private:
        TypedDistribution< AbstractHomologousDiscreteCharacterData >& base_distribution;
    };


}

#endif