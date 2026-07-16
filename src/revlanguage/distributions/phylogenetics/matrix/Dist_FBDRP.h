#ifndef Dist_FBDRP_H
#define Dist_FBDRP_H

#include "FossilizedBirthDeathRangeProcess.h"
#include "RlFossilizedBirthDeathRangeProcess.h"
#include "ModelVector.h"
#include "RlMatrixReal.h"

namespace RevLanguage {

    /**
     * The RevLanguage wrapper of the Fossilized-Birth-Death Range matrix SKELETON (dnFBDRP)
     *
     * The bare birth-death-range process over a matrix of (birth, death) times: the diversification
     * half of the fossilized-birth-death range model. It carries the fossil occurrences (they bound
     * the ranges) and the augmented extreme ages tau_1/tau_last, but NOT the fossil-record term
     * Pr(occurrences | skeleton) -- that is a separate downstream dnFossilRecord node, which also
     * supplies the reporting model. Hence this distribution takes no reporting arguments.
     *
     * Registered as dnFossilizedBirthDeathRange (canonical) and dnFBDRP. dnFBDRMatrix is the
     * deprecated fused form (skeleton + record in one node), which the canonical name used to
     * refer to.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-01-26, version 1.0
     *
     */
    class Dist_FBDRP : public FossilizedBirthDeathRangeProcess<MatrixReal> {

    public:
        Dist_FBDRP( void );

        // Basic utility functions
        Dist_FBDRP*                                             clone(void) const;                                                                      //!< Clone the object
        static const std::string&                               getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                                  getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::vector<std::string>                                getDistributionFunctionAliases(void) const;                                             //!< Get the alternative names used for the constructor function in Rev.
        std::string                                             getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                         getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                                      getParameterRules(void) const;                                                          //!< Get member rules (const)


        // Distribution functions you have to override
        RevBayesCore::FossilizedBirthDeathRangeProcess*         createDistribution(void) const;

    protected:

        void                                                    setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable

        RevPtr<const RevVariable>                               bds;
        RevPtr<const RevVariable>                               origin;
    };

}

#endif
