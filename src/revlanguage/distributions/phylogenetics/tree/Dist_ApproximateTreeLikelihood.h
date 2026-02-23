#ifndef Dist_ApproximateTreeLikelihood_H
#define Dist_ApproximateTreeLikelihood_H

#include "RlTypedDistribution.h"
#include "RlBranchLengthTree.h"
#include "ApproximateTreeLikelihood.h"

namespace RevLanguage {

    /**
     * The RevLanguage wrapper of the ApproximateTreeLikelihood "distribution"
     *
     * The RevLanguage wrapper of the approximate tree likelihood function connects
     * its variables/parameters between the Rev interface and the Core backend, and creates
     * an internal ApproximateTreeLikelihood object.
     *
     *
     * @copyright Copyright 2009-
     * @author David Cerny, Sebastian Hoehna
     * @since 2026-02-23, version 1.4.0-preview
     *
     */
    class Dist_ApproximateTreeLikelihood : public TypedDistribution<BranchLengthTree> {

    public:
        Dist_ApproximateTreeLikelihood(void);

        // Basic utility functions
        Dist_ApproximateTreeLikelihood*                     clone(void) const;                                                                      //!< Clone the object
        static const std::string&                           getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                              getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                         getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                     getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                                  getParameterRules(void) const;                                                          //!< Get member rules (const)

        // Distribution functions you have to override
        RevBayesCore::ApproximateTreeLikelihood*            createDistribution(void) const;

    protected:
        void                                                setConstParameter(const std::string& name, const RevPtr<const RevVariable>& var);       //!< Set member variable

    private:
        RevPtr<const RevVariable>                           time_tree;
        RevPtr<const RevVariable>                           branch_rates;
        RevPtr<const RevVariable>                           gradients;
        RevPtr<const RevVariable>                           hessian;
        RevPtr<const RevVariable>                           transform;
    };
}

#endif
