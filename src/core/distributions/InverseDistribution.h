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
     * @since 2024-09-17, version 1.2.5
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
            for (const auto& parameter : dist->getParameters()) {
                // Doesn't just add the parameter, but also links it so that updates cause
                // the distribution to update
                this->addParameter( parameter );
            }
            
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
            for (const auto& parameter : dist->getParameters()) {
                this->addParameter( parameter );
            }
        }

        // Functions from Distribution.cpp

        // Overriding `touchSpecialization` was found to be necessary (at 2024-10): 
        // nodes in the associated tree were not being marked as dirty
        // when the value of the distribution was changed indirectly.
        // The override seems not to be necessary for straightforward cases
        // such as:
        //     inv_exp ~ dnInv(dnExp(1))
        //     inv_exp.clamp(42)  # Touched from clamp and from invalidated children
        // or:
        //     k ~ dnUniform(1, 10)
        //     kRates := [ [ 0, k],
        //                 [ 1, 0] ]
        //     qK := fnFreeK(kRates)
        //     invPhy ~ dnInv(dnPhyloCTMC(tree = phylogeny, Q = qK, type = "Standard"))
        //     invPhy.clamp(chars) # Touched likewise
        //     k.clamp(1) # Touched likewise
        // 
        // The override becomes necessary for subsequent calls to e.g.
        //     k.clamp(2)
        //
        // Further investigation into this latter case may yield a better solution.
        void touchSpecialization( const DagNode *toucher, bool touchAll ) override
        {
            dist->touch(toucher, touchAll);
        }

        void keepSpecialization( const DagNode* affecter ) override
        {
            dist->keep(affecter);
        }

        void restoreSpecialization( const DagNode *restorer ) override
        {
            dist->restore(restorer);
        }

        void getAffected(RbOrderedSet<DagNode *> &affected, const DagNode* affecter) override
        {
            dist->getAffected(affected, affecter);
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

            if (this->value != v) {
                delete this->value;
                if constexpr(std::is_base_of_v<Cloneable, valType>) {
                    this->value = dist->getValue().clone();
                } else {
                    this->value = new valType(dist->getValue());
                }
            }
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
            dist->redrawValue();
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