#ifndef Dist_conditionedBirthDeathShiftProcessContinuous_H
#define Dist_conditionedBirthDeathShiftProcessContinuous_H

#include "ConditionedBirthDeathShiftProcessContinuous.h"
#include "RlTimeTree.h"
#include "RlTypedDistribution.h"

namespace RevLanguage {
    
    /**
     * @file
     * This file contains the declaration of the Rev conditioned birth-death-shift tree distribution
     *
     *
     * @author Sebastian Hoehna
     */
    class Dist_conditionedBirthDeathShiftProcessContinuous : public TypedDistribution<TimeTree> {
        
    public:
        Dist_conditionedBirthDeathShiftProcessContinuous( void );                                                                                               //!< Constructor
        virtual ~Dist_conditionedBirthDeathShiftProcessContinuous();                                                                                            //!< Virtual destructor
        
        // Basic utility functions
        Dist_conditionedBirthDeathShiftProcessContinuous*               clone(void) const;                                                                      //!< Clone the object
        static const std::string&                                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                                     getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        virtual MethodTable                                             getDistributionMethods( void ) const;                                                   //!< Get the member methods
        const TypeSpec&                                                 getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                                              getParameterRules(void) const;                                                          //!< Get member rules (const)
        
        
        // Distribution functions you have to override
        RevBayesCore::ConditionedBirthDeathShiftProcessContinuous*      createDistribution(void) const;                                                         //!< Create the core object corresponding to this wrapper
        
    protected:
        
        void                                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:
        
        RevPtr<const RevVariable>                                       root_age;
        RevPtr<const RevVariable>                                       root_speciation;
        RevPtr<const RevVariable>                                       root_extinction;
        RevPtr<const RevVariable>                                       speciation;
        RevPtr<const RevVariable>                                       extinction;
        RevPtr<const RevVariable>                                       event_rate;
        RevPtr<const RevVariable>                                       rho;
        RevPtr<const RevVariable>                                       condition;                                                                              //!< The condition of the process (none/survival/#Taxa)
        RevPtr<const RevVariable>                                       taxa;
        
    };
    
}

#endif

