#ifndef Dist_FBDRange_H
#define Dist_FBDRange_H

#include "FossilizedBirthDeathProcess.h"
#include "RlFossilizedBirthDeathRangeProcess.h"
#include "RlTimeTree.h"

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the Fossilized-Birth-Death Range Process
     *
     * The RevLanguage wrapper of the fossilzed-birth-death process connects
     * the variables/parameters of the process and creates the internal FossilizedBirthDeathRangeProcess object.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-01-26, version 1.0
     *c
     */
    class Dist_FBDRange : public FossilizedBirthDeathRangeProcess<TimeTree> {
        
    public:
        Dist_FBDRange( void );
        
        // Basic utility functions
        Dist_FBDRange*                                          clone(void) const;                                                                      //!< Clone the object
        static const std::string&                               getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                                  getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::vector<std::string>                                getDistributionFunctionAliases(void) const;                                             //!< Get the alternative names used for the constructor function in Rev.
        std::string                                             getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                         getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                                      getParameterRules(void) const;                                                          //!< Get member rules (const)
        
        
        // Distribution functions you have to override
        RevBayesCore::FossilizedBirthDeathProcess*              createDistribution(void) const;
        
    protected:
        
        void                                                    setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:

        RevPtr<const RevVariable>                               lambda_a;                                                                               //!< The anagenetic speciation rate(s)
        RevPtr<const RevVariable>                               beta;                                                                                   //!< The symmetric speciation probability
        RevPtr<const RevVariable>                               extended;                                                                               //!< The flag determining whether to treat tips as sampling or extinction events
        RevPtr<const RevVariable>                               start_age;                                                                              //!< The start time of the process since the origin
        std::string                                             start_condition;                                                                        //!< The start condition of the process (rootAge/originAge)

    };
    
}

#endif
