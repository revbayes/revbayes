#ifndef Dist_BranchRateTree_H
#define Dist_BranchRateTree_H

#include "RlTypedDistribution.h"
#include "RlTimeTree.h"
#include "BranchRateTreeDistribution.h"

namespace RevLanguage {

    /**
     * The RevLanguage wrapper of the Dist_BranchRateTree
     *
     * The RevLanguage wrapper of the constrained-topology distribution connects
     * the variables/parameters of the process and creates the internal Dist_BranchRateTree object.
     * Please read the Dist_BranchRateTree.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-01-26, version 1.0
     *
     */
    class Dist_BranchRateTree : public TypedDistribution<Tree> {

    public:
        Dist_BranchRateTree(void);

        // Basic utility functions
        Dist_BranchRateTree*                                clone(void) const;                                                                      //!< Clone the object
        static const std::string&                           getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                              getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::vector<std::string>                            getDistributionFunctionAliases(void) const;                                             //!< Get the alternative names used for the constructor function in Rev.
        std::string                                         getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        virtual MethodTable                                 getDistributionMethods( void ) const;                                                   //!< Get the member methods
        const TypeSpec&                                     getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                                  getParameterRules(void) const;                                                          //!< Get member rules (const)

        // Distribution functions you have to override
        RevBayesCore::BranchRateTreeDistribution*           createDistribution(void) const;

    protected:
        void                                                setConstParameter(const std::string& name, const RevPtr<const RevVariable>& var);       //!< Set member variable

    private:
        RevPtr<const RevVariable>                           branch_rate_prior;
        RevPtr<const RevVariable>                           time_tree;
    };
}

#endif
