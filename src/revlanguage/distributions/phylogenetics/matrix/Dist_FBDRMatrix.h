#ifndef Dist_FBDRMatrix_H
#define Dist_FBDRMatrix_H

#include "FossilizedBirthDeathRangeProcess.h"
#include "RlFossilizedBirthDeathRangeProcess.h"
#include "ModelVector.h"
#include "RlMatrixReal.h"

namespace RevLanguage {

    /**
     * The RevLanguage wrapper of the fused Fossilized-Birth-Death Range Matrix Process (deprecated)
     *
     * The fused form of the fossilized-birth-death range process: a single node carrying both the
     * birth-death range process and the fossil-record term Pr(occurrences | ranges), taking the
     * occurrences as a constructor argument rather than as clamped data. Deprecated in favour of
     * dnFBDRP (the range process) plus dnFossilRecord (the record), and kept only so that existing
     * scripts continue to run.
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

        RevPtr<const RevVariable>                               origin;
    };

}

#endif
