#ifndef Dist_FastBirthDeathShiftProcess_H
#define Dist_FastBirthDeathShiftProcess_H

#include "RlTimeTree.h"
#include "RlTypedDistribution.h"

namespace RevLanguage {
    
    /**
     * @file
     * This file contains the declaration of the Rev character-dependent birth-death tree distribution
     *
     *
     * @author Sebastian Hoehna
     */
    class Dist_FastBirthDeathShiftProcess : public TypedDistribution<TimeTree> {
        
    public:
        Dist_FastBirthDeathShiftProcess( void );        //!< Constructor
        virtual ~Dist_FastBirthDeathShiftProcess();             //!< Virtual destructor
        
        // Basic utility functions
        Dist_FastBirthDeathShiftProcess*                   clone(void) const;                                                                      //!< Clone the object
        static const std::string&                                   getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                                      getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::vector<std::string>                                    getDistributionFunctionAliases(void) const;
        std::string                                                 getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        virtual MethodTable                                         getDistributionMethods( void ) const;                                                                       //!< Get the member methods
        const TypeSpec&                                             getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                                          getParameterRules(void) const;                                                          //!< Get member rules (const)
        
        
        // Distribution functions you have to override
        RevBayesCore::TypedDistribution<RevBayesCore::Tree>*        createDistribution(void) const;                                                         //!< Create the core object corresponding to this wrapper
        
    protected:
        
        void                                                        setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:
        
        RevPtr<const RevVariable>                                   start_age;
        RevPtr<const RevVariable>                                   root_frequencies;
        RevPtr<const RevVariable>                                   speciation_rates;
        RevPtr<const RevVariable>                                   extinction_rates;
        RevPtr<const RevVariable>                                   alpha;
        RevPtr<const RevVariable>                                   beta;
        RevPtr<const RevVariable>                                   rho;
        RevPtr<const RevVariable>                                   condition;                                                                              //!< The condition of the process (none/survival/#Taxa)
        RevPtr<const RevVariable>                                   allow;                                                                              //!< The condition of the process (none/survival/#Taxa)
        RevPtr<const RevVariable>                                   num_time_slices;
        RevPtr<const RevVariable>                                   simulation_condition;
        RevPtr<const RevVariable>                                   min_lineages;
        RevPtr<const RevVariable>                                   max_lineages;
        RevPtr<const RevVariable>                                   exact_lineages;
        RevPtr<const RevVariable>                                   max_time;
        RevPtr<const RevVariable>                                   prune_extinct_lineages;
        std::string                                                 start_condition;                                                                     //!< The start condition of the process (rootAge/originAge)

    };
    
}

#endif
