#ifndef Dist_UltrametricTree_H
#define Dist_UltrametricTree_H

#include "RlTypedDistribution.h"
#include "RlTimeTree.h"
#include "UltrametricTreeDistribution.h"

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the Dist_UltrametricTree
     *
     * The RevLanguage wrapper of the constrained-topology distribution connects
     * the variables/parameters of the process and creates the internal Dist_UltrametricTree object.
     * Please read the Dist_UltrametricTree.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-01-26, version 1.0
     *
     */
    class Dist_UltrametricTree : public TypedDistribution<TimeTree> {
        
    public:
        Dist_UltrametricTree(void);
        
        // Basic utility functions
        Dist_UltrametricTree*                               clone(void) const;                                                                      //!< Clone the object
        static const std::string&                           getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                              getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::vector<std::string>                            getDistributionFunctionAliases(void) const;                                             //!< Get the alternative names used for the constructor function in Rev.
        std::string                                         getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        virtual MethodTable                                 getDistributionMethods( void ) const;                                                                       //!< Get the member methods
        const TypeSpec&                                     getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                                  getParameterRules(void) const;                                                          //!< Get member rules (const)
        
        // Distribution functions you have to override
        RevBayesCore::UltrametricTreeDistribution*          createDistribution(void) const;
        
    protected:
        void                                                setConstParameter(const std::string& name, const RevPtr<const RevVariable>& var);       //!< Set member variable
        
    private:
        RevPtr<const RevVariable>                           tree_prior;
        RevPtr<const RevVariable>                           rate_prior;
        RevPtr<const RevVariable>                           trees;
        RevPtr<const RevVariable>                           root_age;
        RevPtr<const RevVariable>                           root_branch_fraction;
        RevPtr<const RevVariable>                           density;
    };
}

#endif

