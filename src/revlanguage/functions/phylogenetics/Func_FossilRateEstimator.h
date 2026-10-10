#ifndef Func_FossilRateEstimator_H
#define Func_FossilRateEstimator_H

#include <string>

#include "FossilRateEstimatorFunction.h"
#include "RlMatrixReal.h"
#include "RlTypedFunction.h"

namespace RevLanguage {

    /**
     * @brief Shared base for the classical per-interval rate estimators.
     *
     * The three methods take the same arguments and differ only in the formula, so they share
     * their argument rules and their createFunction, and each derived class supplies its Rev
     * name and its method. The value is a two by l matrix, speciation then extinction.
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
