//
//  Dist_sampledSpeciationBirthDeathProcess.h
//  revbayes-proj
//
//  Created by Michael Landis on 2/10/16.
//  Copyright © 2016 Michael Landis. All rights reserved.
//

#ifndef Dist_sampledSpeciationBirthDeathProcess_h
#define Dist_sampledSpeciationBirthDeathProcess_h

#include "SampledSpeciationBirthDeathProcess.h"
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
    class Dist_sampledSpeciationBirthDeathProcess : public TypedDistribution<TimeTree> {
        
    public:
        Dist_sampledSpeciationBirthDeathProcess( void );        //!< Constructor
        virtual ~Dist_sampledSpeciationBirthDeathProcess();             //!< Virtual destructor
        
        // Basic utility functions
        Dist_sampledSpeciationBirthDeathProcess*               clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                     getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        virtual MethodTable                             getDistributionMethods( void ) const;                                                                       //!< Get the member methods
        const TypeSpec&                                 getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)
        
        
        // Distribution functions you have to override
        RevBayesCore::SampledSpeciationBirthDeathProcess*      createDistribution(void) const;                                                         //!< Create the core object corresponding to this wrapper
        
    protected:
        
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:
        
        RevPtr<const RevVariable>                       root_age;
        RevPtr<const RevVariable>                       speciation;
        RevPtr<const RevVariable>                       extinction;
        RevPtr<const RevVariable>                       rho;
        RevPtr<const RevVariable>                       taxa;
        
    };
    
}


#endif /* Dist_sampledSpeciationBirthDeathProcess_h */
