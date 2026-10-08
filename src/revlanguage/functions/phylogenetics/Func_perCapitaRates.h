#ifndef Func_perCapitaRates_H
#define Func_perCapitaRates_H

#include <string>

#include "Func_FossilRateEstimator.h"

namespace RevLanguage {

    /**
     * @brief Foote's per-taxon rates: the proportion of first or last appearances in each interval.
     */
    class Func_perCapitaRates : public Func_FossilRateEstimator {

    public:
        Func_perCapitaRates( void );

        Func_perCapitaRates*                                  clone(void) const override;
        static const std::string&               getClassType(void);
        static const TypeSpec&                  getClassTypeSpec(void);
        std::string                             getFunctionName(void) const override;
        const TypeSpec&                         getTypeSpec(void) const override;

    protected:
        RevBayesCore::FossilRateEstimatorFunction::Method method(void) const override;
    };

}

#endif
