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
    class InverseDistribution : public TypedDistribution< valType > {
     
    public:
        // constructor(s)
        InverseDistribution(const TypedDistribution<valType>& d)
            : TypedDistribution<valType>( nullptr ),
              dist( d.clone() )
        {
            // add the parameters to our set (in the base class)
            // in that way other class can easily access the set of our parameters
            // this will also ensure that the parameters are not getting deleted before we do
            
            // add the parameters of the distribution
            for (const auto& parameter : dist->getParameters())
                this->addParameter( parameter );
            
	        redrawValue();
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

        // functions from 'public methods' section of TypedDistribution.h
        
        template<class variableType>
        const RevBayesCore::RbVector<variableType>& getParameterValues(void) const {
            return dist->getParameterValues();
        }
        
        template<class variableType>
        variableType& getValue( void ) {
            return dist->getValue();
        }

        template<class variableType>
        const variableType& getValue( void ) const {
            return dist->getValue();
        }        
        

        // functions from 'virtual methods' section of TypedDistribution.h

        // Set the current value, e.g. attach an observation (clamp)
        void setValue(valType *v, bool f=false) override {
            dist->setValue(v, f);
        }

        // functions from 'pure virtual methods' section of TypedDistribution.h
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
	        delete this->value;

            if constexpr(std::is_base_of_v<Cloneable, valType>) {
                this->value = dist->getValue().clone();
            } else {
                this->value = new valType(dist->getValue());
            }
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
