#ifndef Func_boundaryCrosserRates_H
#define Func_boundaryCrosserRates_H

#include <string>

#include "Func_FossilRateEstimator.h"

namespace RevLanguage {

    /**
     * @brief Foote's boundary-crosser rates, which drop the taxa confined to one interval.
     */
    class Func_boundaryCrosserRates : public Func_FossilRateEstimator {

    public:
        Func_boundaryCrosserRates( void );

        Func_boundaryCrosserRates*                                  clone(void) const override;
        static const std::string&               getClassType(void);
        static const TypeSpec&                  getClassTypeSpec(void);
        std::string                             getFunctionName(void) const override;
        const TypeSpec&                         getTypeSpec(void) const override;

    protected:
        RevBayesCore::FossilRateEstimatorFunction::Method method(void) const override;
    };

}

#endif
