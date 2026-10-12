#ifndef Func_perTaxonRates_H
#define Func_perTaxonRates_H

#include <string>

#include "Func_FossilRateEstimator.h"

namespace RevLanguage {

    /**
     * @brief Foote's per-taxon rates: the proportion of first or last appearances in each interval.
     */
    class Func_perTaxonRates : public Func_FossilRateEstimator {

    public:
        Func_perTaxonRates( void );

        Func_perTaxonRates*                     clone(void) const override;
        static const std::string&               getClassType(void);
        static const TypeSpec&                  getClassTypeSpec(void);
        std::string                             getFunctionName(void) const override;
        const TypeSpec&                         getTypeSpec(void) const override;

    protected:
        RevBayesCore::FossilRateEstimatorFunction::Method method(void) const override;
    };

}

#endif
