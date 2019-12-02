#ifndef Dist_heterogeneousRateBirthDeath_H
#define Dist_heterogeneousRateBirthDeath_H

#include "HeterogeneousRateBirthDeath.h"
#include "RlTimeTree.h"
#include "RlTypedDistribution.h"

namespace RevLanguage {
    
    /**
     * @file
     * This file contains the declaration of the Rev heterogeneous-rate birth-death tree distribution
     *
     *
     * @author Sebastian Hoehna
     */
    class Dist_heterogeneousRateBirthDeath : public TypedDistribution<TimeTree> {
        
    public:
        Dist_heterogeneousRateBirthDeath( void );        //!< Constructor
        virtual ~Dist_heterogeneousRateBirthDeath();             //!< Virtual destructor
        
        // Basic utility functions
        Dist_heterogeneousRateBirthDeath*               clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                     getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        virtual MethodTable                             getDistributionMethods( void ) const;                                                                       //!< Get the member methods
        const TypeSpec&                                 getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)
        
        
        // Distribution functions you have to override
        RevBayesCore::HeterogeneousRateBirthDeath*      createDistribution(void) const;                                                         //!< Create the core object corresponding to this wrapper
        
    protected:
        
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:
        
        RevPtr<const RevVariable>                       allow_same_category;
        RevPtr<const RevVariable>                       root_age;
        RevPtr<const RevVariable>                       root_state;
        RevPtr<const RevVariable>                       speciation;
        RevPtr<const RevVariable>                       extinction;
        RevPtr<const RevVariable>                       event_rate;
        RevPtr<const RevVariable>                       rho;
        RevPtr<const RevVariable>                       condition;                                                                              //!< The condition of the process (none/survival/#Taxa)
        RevPtr<const RevVariable>                       taxa;
    
    };
    
}

#endif
