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

    template<typename charType>
    class InversePhyloDistribution : public TypedDistribution< AbstractHomologousDiscreteCharacterData > {

    public:
        // constructor(s)
        // Notice that we don't need a version of
        // InverseDistribution(const TypedDistribution<valType>& d) : TypedDistribution<valType>( new valType() ),
        // as we can't call new AbstractHomologousDiscreteCharacterData() - it's an abstract function
        InversePhyloDistribution(TypedDistribution< AbstractHomologousDiscreteCharacterData >& d)
            : TypedDistribution< AbstractHomologousDiscreteCharacterData >( d ), base_distribution( d.clone() ) {
                
            // add the parameters of the distribution
            for (const auto& parameter : base_distribution->getParameters())
                this->addParameter( parameter );
        }

        // Virtual destructor
        virtual ~InversePhyloDistribution(void) = default;
        
        // Set the current value, e.g. attach an observation (clamp)
        void setValue(valType *v, bool f = false) override {
            base_distribution->setValue(v, f);
        }

        // public member functions
        
        InversePhyloDistribution* clone(void) const override // Create an independent clone
        {
            return new InversePhyloDistribution( *base_distribution );
        }
        
        double computeLnProbability(void) override {
            return -( base_distribution->computeLnProbability() );
        }

        void redrawValue(void) override {
            redrawValue();
        }
        
    protected:
        // Parameter management functions
        void swapParameterInternal(const DagNode *oldP, const DagNode *newP) override
        {
            base_distribution->swapParameter( oldP, newP );
        }


    private:
        std::unique_ptr<TypedDistribution< AbstractHomologousDiscreteCharacterData >> base_distribution;
    };


}

#endif