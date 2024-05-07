#include <cstddef>
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
#include "Cloneable.h"
#include "DiscreteCharacterState.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "TimeReversibleRateMatrix.h"

using namespace RevBayesCore;
using std::vector;

namespace RevBayesCore
{
vector<double> compute_flattened_exchange_rates( const MatrixReal& Q, const vector<double>& pi )
{
    int n = pi.size();

    assert(Q.getNumberOfRows() == n);
    assert(Q.getNumberOfColumns() == n);

    vector<double> er(n*(n-1)/2);

    int k=0;
    for(int i=0; i<n; i++)
    {
        assert( pi[i] >= 0.0 );

        for(int j=i+1; j<n; j++)
        {
            assert( Q[i][j] >= 0.0 );

#ifndef NDEBUG
            double relative_to = std::abs( Q[i][j]*pi[j] ) + std::abs( Q[j][i]*pi[i] );
            double reversibility_error = std::abs( pi[i]*Q[i][j] - pi[j]*Q[j][i] );
            assert( reversibility_error <= relative_to * 1.0e-9 );
#endif
            if (pi[j] > 0)
                er[k++] = Q[i][j] / pi[j];
            else
                er[k++] = 0;
        }
    }

    return er;
}

vector<double> flatten_exchange_rates( const MatrixReal& ER )
{
    int n = ER.getNumberOfRows();

    assert(ER.getNumberOfColumns() == n);

    vector<double> er(n*(n-1)/2);

    int k=0;
    for(int i=0; i<n; i++)
    {
        for(int j=i+1; j<n; j++)
        {
            er[k++] = ER[i][j];
        }
    }

    assert( k == n*(n-1) / 2 );
    
    return er;
}
}

/** Construct rate matrix with n states */
ConcreteTimeReversibleRateMatrix::ConcreteTimeReversibleRateMatrix( const vector<double>& er, const vector<double>& pi, boost::optional<double> r)
    : TimeReversibleRateMatrix( pi.size() ), _rate(r)
{
    int n = pi.size();
    assert( er.size() == n*(n-1)/2 );

    setExchangeabilityRates(er);
    setStationaryFrequencies(pi);
    update();
}

ConcreteTimeReversibleRateMatrix::ConcreteTimeReversibleRateMatrix( const MatrixReal& ER, const vector<double>& pi, boost::optional<double> r)
    : ConcreteTimeReversibleRateMatrix( flatten_exchange_rates(ER), pi, r)
{
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
        if (_rate)
            rescaleToAverageRate( *_rate );

        // clean flags
        needs_update = false;
    }
}


