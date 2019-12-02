#ifndef RlBirthDeathProcess_H
#define RlBirthDeathProcess_H

#include "AbstractBirthDeathProcess.h"
#include "RlTimeTree.h"
#include "RlTypedDistribution.h"

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the BirthDeathProcess.
     *
     * The RevLanguage wrapper of the birth-death process is the base-class for all
     * wrappers of homogeneous (i.e. each species has the same rate) birth-death process.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-01-26, version 1.0
     *
     */
    class BirthDeathProcess : public TypedDistribution<TimeTree> {
        
    public:
        BirthDeathProcess( void );
        
        // Basic utility functions you have to overwrite
        virtual BirthDeathProcess*                          clone(void) const = 0;                                                              //!< Clone the object

        // Basic utility functions you may want to overwrite
        const MemberRules&                                  getParameterRules(void) const;                                                      //!< Get member rules (const)
        
        // Basic utility functions
        static const std::string&                           getClassType(void);                                                                 //!< Get Rev type
        static const TypeSpec&                              getClassTypeSpec(void);                                                             //!< Get class type spec
        
        
        // Distribution functions you have to override
        virtual RevBayesCore::AbstractBirthDeathProcess*    createDistribution(void) const = 0;                                                 //!< Create an internal object of the diveristy-dependent pure-birth process.
        
    protected:
        
        void                                                setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);   //!< Set member variable
    
        // members        
        RevPtr<const RevVariable>                           condition;                                                                          //!< The condition of the process (none/survival/#Taxa)
        RevPtr<const RevVariable>                           incomplete_clades;                                                                  //!< The rev taxon objects that will be applied to the initally simulated tree
        RevPtr<const RevVariable>                           sampling_strategy;                                                                   //!< The applied sampling strategy (uniform/diversified)
        RevPtr<const RevVariable>                           start_age;                                                                           //!< The start time of the process since the origin
        RevPtr<const RevVariable>                           taxa;                                                                               //!< The rev taxon objects that will be applied to the initally simulated tree
        
    };
    
}

#endif
