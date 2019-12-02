//
//  Dist_outgroupBirthDeath.h
//  revbayes-proj
//
//  Created by Michael Landis on 5/5/16.
//  Copyright © 2016 Michael Landis. All rights reserved.
//

#ifndef Dist_outgroupBirthDeath_h
#define Dist_outgroupBirthDeath_h

#ifndef Dist_outgroupBirthDeath_H
#define Dist_outgroupBirthDeath_H

#include "ConstantRateOutgroupBirthDeathProcess.h"
#include "RlBirthDeathProcess.h"
#include "RlTimeTree.h"

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the constant-rate outgroup-sampled birth-death process
     *
     * The RevLanguage wrapper of the constant-rate outgroup-sampled process connects
     * the variables/parameters of the process and creates the internal ConstantRateOutgroupBirthDeathProcess object.
     * Please read the ConstantRateFossilizedBirthDeathProcess.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-01-26, version 1.0
     *c
     */
    class Dist_outgroupBirthDeath : public BirthDeathProcess {
        
    public:
        Dist_outgroupBirthDeath( void );
        
        // Basic utility functions
        Dist_outgroupBirthDeath*                                clone(void) const;                                                                      //!< Clone the object
        static const std::string&                               getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                                  getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::vector<std::string>                                getDistributionFunctionAliases(void) const;                                             //!< Get the alternative names used for the constructor function in Rev.
        std::string                                             getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                         getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                                      getParameterRules(void) const;                                                          //!< Get member rules (const)
        
        
        // Distribution functions you have to override
        RevBayesCore::ConstantRateOutgroupBirthDeathProcess*    createDistribution(void) const;
        
    protected:
        
        void                                                    setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:
        
        RevPtr<const RevVariable>                               lambda;                                                                                 //!< The speciation rate
        RevPtr<const RevVariable>                               mu;                                                                                     //!< The extinction rate
        RevPtr<const RevVariable>                               rhoIngroup;                                                                             //!< The ingroup sampling probability
        RevPtr<const RevVariable>                               rhoOutgroup;                                                                            //!< The outgroup sampling probability
        RevPtr<const RevVariable>                               taxaOutgroup;                                                                           //!< The outgroup sampling probability
        RevPtr<const RevVariable>                               samplingStrategyOutgroup;                                                               //!< The outgroup sampling probability
        RevPtr<const RevVariable>                               samplingStrategyIngroup;                                                                //!< The outgroup sampling probability

        
    };
    
}

#endif


#endif /* Dist_outgroupBirthDeath_h */
