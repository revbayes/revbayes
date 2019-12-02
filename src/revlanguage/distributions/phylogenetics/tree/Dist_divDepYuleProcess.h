#ifndef Dist_divDepYuleProcess_H
#define Dist_divDepYuleProcess_H

#include "DiversityDependentPureBirthProcess.h"
#include "RlTypedDistribution.h"
#include "RlTimeTree.h"

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the Dist_divDepYuleProcess
     *
     * The RevLanguage wrapper of the diversity-dependent pure-birth process connects
     * the variables/parameters of the process and creates the internal Dist_divDepYuleProcess object.
     * Please read the Dist_divDepYuleProcess.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-01-26, version 1.0
     *
     */
    class Dist_divDepYuleProcess : public TypedDistribution<TimeTree> {
        
    public:
        Dist_divDepYuleProcess( void );
        
        // Basic utility functions
        Dist_divDepYuleProcess*                             clone(void) const;                                                                      //!< Clone the object
        static const std::string&                           getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                              getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                         getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                     getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                                  getParameterRules(void) const;                                                          //!< Get member rules (const)
        
        
        // Distribution functions you have to override
        RevBayesCore::DiversityDependentPureBirthProcess*   createDistribution(void) const;                                                         //!< Create an internal object of the diveristy-dependent pure-birth process.
        
    protected:
    
        void                                                setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:
        
        // members        
        RevPtr<const RevVariable>                           rootAge;                                                                                //!< The time of the process since the origin
        RevPtr<const RevVariable>                           condition;                                                                              //!< The condition of the process (none/survival/#Taxa)
        RevPtr<const RevVariable>                           taxa;                                                                             //!< The taxon names that will be applied to the initally simulated tree
        RevPtr<const RevVariable>                           initialLambda;
        RevPtr<const RevVariable>                           capacity;
        
    };
    
}

#endif
