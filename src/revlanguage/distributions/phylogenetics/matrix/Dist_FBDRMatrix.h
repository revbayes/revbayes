#ifndef Dist_FBDRMatrix_H
#define Dist_FBDRMatrix_H

#include "FossilizedBirthDeathRangeProcess.h"
#include "RlFossilizedBirthDeathRangeProcess.h"
#include "ModelVector.h"
#include "RlMatrixReal.h"

namespace RevLanguage {

    /**
     * The RevLanguage wrapper of the FUSED Fossilized-Birth-Death Range Matrix Process (DEPRECATED)
     *
     * The historical fused form: one node carrying both the birth-death-range skeleton and the
     * fossil-record term Pr(occurrences | skeleton), with the occurrences smuggled in through the
     * `taxa` constructor argument rather than clamped. It is retained as a deprecated facade,
     * warns, and is byte-identical to the pre-split process. The factored replacement is
     * dnFBDRP (skeleton) + dnFossilRecord (record).
     *
     * Reachable as dnFBDRMatrix only: the canonical dnFossilizedBirthDeathRange name has moved to
     * the skeleton, so a script that called the fused process by its long name and passed no
     * reporting arguments now builds a bare skeleton and needs a dnFossilRecord node added.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-01-26, version 1.0
     *
     */
    class Dist_FBDRMatrix : public FossilizedBirthDeathRangeProcess<MatrixReal> {

    public:
        Dist_FBDRMatrix( void );

        // Basic utility functions
        Dist_FBDRMatrix*                                        clone(void) const;                                                                      //!< Clone the object
        static const std::string&                               getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                                  getClassTypeSpec(void);                                                                 //!< Get class type spec
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
