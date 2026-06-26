#ifndef Dist_BDS_H
#define Dist_BDS_H

#include "FossilizedBirthDeathRangeProcess.h"
#include "RlFossilizedBirthDeathRangeProcess.h"
#include "ModelVector.h"
#include "RlMatrixReal.h"

namespace RevLanguage {

    /**
     * The RevLanguage wrapper of the BDS (birth-death-sampling) matrix process of
     * Silvestro et al. (2019).
     *
     * This is the same matrix process as the fossilized-birth-death range process
     * (Dist_FBDRP), but with the BDS likelihood always enabled: it exposes the same
     * parameters minus the BDS flag and simply passes use_bds=true into the core
     * FossilizedBirthDeathRangeProcess constructor.
     */
    class Dist_BDS : public FossilizedBirthDeathRangeProcess<MatrixReal> {

    public:
        Dist_BDS( void );

        // Basic utility functions
        Dist_BDS*                                       clone(void) const;                                  //!< Clone the object
        static const std::string&                       getClassType(void);                                 //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                             //!< Get class type spec
        std::string                                     getDistributionFunctionName(void) const;            //!< Get the Rev-name for this distribution.
        const TypeSpec&                                 getTypeSpec(void) const;                            //!< Get the type spec of the instance

        // Distribution functions you have to override
        RevBayesCore::FossilizedBirthDeathRangeProcess* createDistribution(void) const;
    };

}

#endif
