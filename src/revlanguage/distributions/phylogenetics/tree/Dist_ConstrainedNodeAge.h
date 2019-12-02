#ifndef Dist_ConstrainedNodeAge_H
#define Dist_ConstrainedNodeAge_H

#include "NodeAgeConstrainedTreeDistribution.h"
#include "RlTypedDistribution.h"
#include "RlTimeTree.h"

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the Dist_ConstrainedNodeAge
     *
     * The RevLanguage wrapper of the ultrametric tree distribution with constrained node age connects
     * the variables/parameters of the process and creates the internal Dist_ConstrainedNodeAge object.
     * Please read the Dist_ConstrainedNodeAge.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-01-26, version 1.0
     *
     */
    class Dist_ConstrainedNodeAge : public TypedDistribution<TimeTree> {
        
    public:
        Dist_ConstrainedNodeAge(void);
        
        // Basic utility functions
        Dist_ConstrainedNodeAge*                            clone(void) const;                                                                      //!< Clone the object
        static const std::string&                           getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                              getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::vector<std::string>                            getDistributionFunctionAliases(void) const;                                             //!< Get the alternative names used for the constructor function in Rev.
        std::string                                         getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                     getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                                  getParameterRules(void) const;                                                          //!< Get member rules (const)
        
        
        // Distribution functions you have to override
        RevBayesCore::NodeAgeConstrainedTreeDistribution*   createDistribution(void) const;
        
    protected:
        void                                                setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
    private:
        
        RevPtr<const RevVariable>                           baseDistribution;
        RevPtr<const RevVariable>                           constraints;                                                                            //!< The base tree distribution
    };
}

#endif
