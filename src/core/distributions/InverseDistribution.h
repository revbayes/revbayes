#ifndef InverseDistribution_h
#define InverseDistribution_h

#include "TypedDistribution.h"
#include "TypedDagNode.h"
#include <memory>

namespace RevBayesCore {
    
    /**
     * Provides the inverse of the probability of a distribution.
     *
     * @copyright Copyright 2024-
     * @author Martin R. Smith
     * @since 2024-07-14, version 1.2.5
     *
     */
    template<typename valType>
    class InverseDistribution :  public TypedDistribution< valType > {
     
    public:
        // constructor(s)
        InverseDistribution(const TypedDistribution<valType>& d)
            : TypedDistribution<valType>( new valType() ),
              dist( d.clone() )
        {
            // add the parameters to our set (in the base class)
            // in that way other class can easily access the set of our parameters
            // this will also ensure that the parameters are not getting deleted before we do
            
            // add the parameters of the distribution
            for (const auto& parameter : dist->getParameters())
                this->addParameter( parameter );
            
            dist->redrawValue();
        }

        InverseDistribution(const InverseDistribution &d)
            : TypedDistribution<valType>( d ),
              dist( d.dist->clone() )
        {
            // add the parameters to our set (in the base class)
            // in that way other class can easily access the set of our parameters
            // this will also ensure that the parameters are not getting deleted before we do
            
            // add the parameters of the distribution
            for (const auto& parameter : dist->getParameters())
                this->addParameter( parameter );
        }

        // Set the current value, e.g. attach an observation (clamp)
        void setValue(valType *v, bool f=false) override {
            dist->setValue(v, f);
        }

        // public member functions
        InverseDistribution* clone(void) const override // Create an independent clone
        {
            return new InverseDistribution( *this );
        }

        double computeLnProbability(void) override
        {
            return -(dist->computeLnProbability());
        }

        void redrawValue(void) override
        {
            dist->redrawValue();
        }
        
    protected:
        // Parameter management functions
        void swapParameterInternal(const DagNode *oldP, const DagNode *newP) override
        {
            dist->swapParameter( oldP, newP );
        }
        
    private:        
        // private members
        std::unique_ptr<TypedDistribution<valType>> dist;
    };
}

#endif /* InverseDistribution_h */