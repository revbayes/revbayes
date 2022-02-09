#include <stddef.h>
#include <cmath>
#include <complex>
#include <iosfwd>
#include <vector>

#include "CodonState.h"
#include "MatrixComplex.h"
#include "MatrixReal.h"
#include "ConcreteTimeReversibleRateMatrix.h"
#include "RbException.h"
#include "TransitionProbabilityMatrix.h"
#include "AminoAcidState.h"
#include "Assignable.h"
#include "Cloneable.h"
#include "DiscreteCharacterState.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "TimeReversibleRateMatrix.h"

using namespace RevBayesCore;
using std::vector;

/** Construct rate matrix with n states */
ConcreteTimeReversibleRateMatrix::ConcreteTimeReversibleRateMatrix( const vector<double>& er, const vector<double>& pi)
    : TimeReversibleRateMatrix( pi.size() )
{
    int n = pi.size();
    assert( er.size() == n*(n-1)/2 );

    setExchangeabilityRates(er);
    setStationaryFrequencies(pi);
    update();
}

ConcreteTimeReversibleRateMatrix& ConcreteTimeReversibleRateMatrix::assign(const Assignable &m)
{
    const ConcreteTimeReversibleRateMatrix *rm = dynamic_cast<const ConcreteTimeReversibleRateMatrix*>(&m);
    if ( rm != NULL )
    {
        return operator=(*rm);
    }
    else
    {
        throw RbException("Could not assign rate matrix.");
    }
}

ConcreteTimeReversibleRateMatrix* ConcreteTimeReversibleRateMatrix::clone( void ) const
{
    return new ConcreteTimeReversibleRateMatrix( *this );
}

/** Calculate the transition probabilities */
void ConcreteTimeReversibleRateMatrix::calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const
{
    double t = rate * (startAge - endAge);
    exponentiateMatrixByScalingAndSquaring(t, P);
}


void ConcreteTimeReversibleRateMatrix::update( void )
{
    if ( needs_update )
    {
        // compute the off-diagonal values
        computeOffDiagonal();

        // compute the diagonal values
        setDiagonal();

        // rescale
        rescaleToAverageRate( 1.0 );

        // clean flags
        needs_update = false;
    }
}


