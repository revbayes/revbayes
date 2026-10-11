#ifndef Func_FossilRateEstimator_H
#define Func_FossilRateEstimator_H

#include <string>

#include "FossilRateEstimatorFunction.h"
#include "RlMatrixReal.h"
#include "RlTypedFunction.h"

namespace RevLanguage {

    /**
     * @brief Shared base for the classical per-interval rate estimators.
     */
    class Func_FossilRateEstimator : public TypedFunction< MatrixReal > {

    public:
        Func_FossilRateEstimator( void );

        const ArgumentRules&                    getArgumentRules(void) const override;                  //!< Get argument rules
        RevBayesCore::TypedFunction< RevBayesCore::MatrixReal >* createFunction(void) const override;   //!< Create the internal function object

    protected:
        virtual RevBayesCore::FossilRateEstimatorFunction::Method method(void) const = 0;
    };

}

#endif
