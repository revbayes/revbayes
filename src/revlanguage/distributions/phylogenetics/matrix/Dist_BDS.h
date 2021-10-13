#ifndef Dist_BDS_H
#define Dist_BDS_H

#include "BirthDeathRateShiftsProcess.h"
#include "ModelVector.h"
#include "RlMatrixReal.h"
#include "RlTypedDistribution.h"

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the Birth-Death with RateShifts Process
     *
     * The RevLanguage wrapper of the birth-death with rateshifts connects
     * the variables/parameters of the process and creates the internal BirthDeathRateShiftsProcess object.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-01-26, version 1.0
     *c
     */
    class Dist_BDS : public TypedDistribution<MatrixReal> {
        
    public:
        Dist_BDS( void );
        
        // Basic utility functions
        Dist_BDS*                                               clone(void) const;                                                                      //!< Clone the object
        static const std::string&                               getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                                  getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::vector<std::string>                                getDistributionFunctionAliases(void) const;                                             //!< Get the alternative names used for the constructor function in Rev.
        std::string                                             getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                         getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                                      getParameterRules(void) const;                                                          //!< Get member rules (const)
        
        
        // Distribution functions you have to override
        RevBayesCore::BirthDeathRateShiftsProcess*              createDistribution(void) const;


    protected:

        void                                                setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);   //!< Set member variable

        // members
        RevPtr<const RevVariable>                           lambda;                                                                             //!< The speciation rate(s)
        RevPtr<const RevVariable>                           mu;                                                                                 //!< The extinction rate(s)
        RevPtr<const RevVariable>                           psi;                                                                                //!< The fossilization rate(s)
        RevPtr<const RevVariable>                           rho;                                                                                //!< The extant sampling proportion
        RevPtr<const RevVariable>                           timeline;                                                                           //!< The interval times
        RevPtr<const RevVariable>                           taxa;                                                                               //!< The taxa
        RevPtr<const RevVariable>                           condition;                                                                          //!< The condition of the process
        RevPtr<const RevVariable>                           complete;
    };
    
}

#endif
