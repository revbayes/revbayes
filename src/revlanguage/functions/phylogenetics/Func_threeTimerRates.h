#ifndef Func_threeTimerRates_H
#define Func_threeTimerRates_H

#include <string>

#include "Func_FossilRateEstimator.h"

namespace RevLanguage {

    /**
     * @brief Alroy's three-timer rates, which correct the counts by a pooled sampling probability.
     */
    class Func_threeTimerRates : public Func_FossilRateEstimator {

    public:
        Func_threeTimerRates( void );

        Func_threeTimerRates*                                  clone(void) const override;
        static const std::string&               getClassType(void);
        static const TypeSpec&                  getClassTypeSpec(void);
        std::string                             getFunctionName(void) const override;
        const TypeSpec&                         getTypeSpec(void) const override;

    protected:
        RevBayesCore::FossilRateEstimatorFunction::Method method(void) const override;
    };

}

#endif
