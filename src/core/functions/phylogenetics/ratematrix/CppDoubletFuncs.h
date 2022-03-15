#ifndef CppCodonFuncs_H
#define CppCodonFuncs_H

#include "Simplex.h"
#include "MatrixReal.h"
#include "ConcreteTimeReversibleRateMatrix.h"

namespace RevBayesCore
{
    typedef ConcreteTimeReversibleRateMatrix CGTR;

    Simplex F2x4(const Simplex& nuc_pi1, const Simplex& nuc_pi2);

    CGTR X2X2(const TimeReversibleRateMatrix& q1, const TimeReversibleRateMatrix& q2);

    CGTR X2(const TimeReversibleRateMatrix& q);
}


#endif
