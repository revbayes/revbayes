#ifndef Dist_FBDSP_H
#define Dist_FBDSP_H

#include "RlBirthDeathProcess.h"

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the constant-rate Fossilized-Birth-Death Process
     *
     * The RevLanguage wrapper of the constant-rate fossilzed-birth-death process connects
     * the variables/parameters of the process and creates the internal ConstantRateFossilizedBirthDeathProcess object.
     * Please read the ConstantRateFossilizedBirthDeathProcess.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-01-26, version 1.0
     *c
     */
    class Dist_FBDSP : public BirthDeathProcess {
        
    public:
        Dist_FBDSP( void );
        
        // Basic utility functions
        Dist_FBDSP*                                             clone(void) const;                                                                       //!< Clone the object
        static const std::string&                               getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                                  getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::vector<std::string>                                getDistributionFunctionAliases(void) const;                                             //!< Get the alternative names used for the constructor function in Rev.
        std::string                                             getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                         getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                                      getParameterRules(void) const;                                                          //!< Get member rules (const)
        
        
        // Distribution functions you have to override
        RevBayesCore::AbstractBirthDeathProcess*                createDistribution(void) const;
        
    protected:
        
        void                                                    setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:

        RevPtr<const RevVariable>                               lambda;                                                                                 //!< The speciation rate(s)
        RevPtr<const RevVariable>                               mu;                                                                                     //!< The extinction rate(s)
        RevPtr<const RevVariable>                               psi;                                                                                    //!< The fossilization rate(s)
        RevPtr<const RevVariable>                               rho;                                                                                    //!< The extant sampling proportion
        RevPtr<const RevVariable>                               timeline;                                                                               //!< The interval times
        RevPtr<const RevVariable>                               fossil_counts;                                                                          //!< The fossil counts
        RevPtr<const RevVariable>                               presence_absence;
        RevPtr<const RevVariable>                               extended;
        std::string                                             start_condition;                                                                        //!< The start condition of the process (rootAge/originAge)
        RevPtr<const RevVariable>                               initial_tree;

    };
    
}

#endif
