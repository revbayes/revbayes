#include "RateMatrix_Rational.h"

#include <cmath>
#include <fstream>
#include <algorithm>
#include <cstddef>

#include "MatrixReal.h"
#include "RandomNumberGenerator.h"
#include "RandomNumberFactory.h"
#include "RbConstants.h"
#include "RbException.h"
#include "RbMathMatrix.h"
#include "DistributionPoisson.h"
#include "TransitionProbabilityMatrix.h"
#include "Cloneable.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "RlUserInterface.h"

using namespace RevBayesCore;

/** Construct rate matrix with n states */
RateMatrix_Rational::RateMatrix_Rational(size_t n) : RateMatrix(n),
//    the_rate_matrix( NULL ),
    needs_update( true )
{

    // I cannot call a pure virtual function from the constructor (Sebastian)
    update();
}



/** Copy constructor */
RateMatrix_Rational::RateMatrix_Rational(const RateMatrix_Rational& m) : RateMatrix(m),
//    the_rate_matrix( NULL ),
    needs_update( true )
{

}


/** Destructor */
RateMatrix_Rational::~RateMatrix_Rational(void)
{

//    delete the_rate_matrix;
}


RateMatrix_Rational& RateMatrix_Rational::operator=(const RateMatrix_Rational &r)
{

    if (this != &r)
    {
        // delegate to parent class
        RateMatrix::operator=( r );

//        delete the_rate_matrix;
//
//        the_rate_matrix      = new MatrixRational( *r.the_rate_matrix );
        needs_update         = true;

    }

    return *this;
}


/** Calculate the average rate for the rate matrix */
double RateMatrix_Rational::averageRate(void) const
{
    
    mpq_class ave = 0;
    for (size_t i=0; i<num_states; ++i)
    {
//        ave += -stationary_frequencies[i] * (*the_rate_matrix)[i][i];
    }
    return ave.get_d();
}


/** This function calculates the stationary frequencies of the rate matrix. The
 rate matrix, Q, is the infinitesimal generator of the Markov chain. It is an
 n X n matrix whose off-diagonal elements are q_ij >= 0 and whose diagonal elements
 are specified such that each row sums to zero. The rate matrix is finite (has
 a fixed number of states) and we assume that the input matrix is irreducible, as
 is the usual case for substitution models. Because Q is irreducible and finite,
 it has a stationary distribution, pi, which is a row vector of n probabilities.
 The stationary probabilities can be calculated by solving the homogeneous system
 of equations, pi*Q = 0, where 0 is a vector of zeros.

 We do the following to calculate the stationary frequencies.

 1. We perform an LU decomposition of the transpose of the matrix Q.

 Q' = LU

 2. Now we set Ux = z (x will eventually hold the stationary probabilities).
 Because L is nonsingular, we have z = 0. We proceed to back substitute on
 Ux = z = 0. When u_nn = 0, we can put in any solution for x. Here, we put
 in x_n = 1. We then solve the other values of x through back substitution.

 3. The solution obtained in 2 is not a probability vector. We normalize the
 vector such that the sum of the elements is 1.

 Note that the only time we need to use this function is when we don't
 know the stationary frequencies of the rate matrix beforehand. For most
 substitution models used in molecular evolution, the stationary frequencies
 are built into the rate matrix itself. These models are time-reversible.
 This function is useful for the non-reversible models.

 For more information on the fascinating topic of calculating the stationary
 probabilities of a rate matrix, see:

 Stewart, W. J. 1999. Numerical methods for computing stationary distributions of
 finite irreducible Markov chains. In "Advances in Computational
 Probability", W. Grassmann, ed. Kluwer Academic Publishers. */
std::vector<double> RateMatrix_Rational::calculateStationaryFrequencies(void) const
{

    // transpose the rate matrix and put into QT
    MatrixReal QT(num_states, num_states);
    for (size_t i=0; i<num_states; i++)
    {
        for (size_t j=0; j<num_states; j++)
        {
//            QT[i][j] = (*the_rate_matrix)[j][i].get_d();
        }
    }

    // compute the LU decomposition of the transposed rate matrix
    MatrixReal L(num_states, num_states);
    MatrixReal U(num_states, num_states);
    RbMath::computeLandU(QT, L, U);

    // back substitute into z = 0 to find un-normalized stationary frequencies, starting with x_n = 1.0
    std::vector<double> pi(num_states, 0.0);
    pi[num_states-1] = 1.0;
    size_t i=num_states-1;
    while ( i > 0 )
    {
        i--;
        double dotProduct = 0.0;
        for (size_t j=i+1; j<num_states; j++)
        {
            dotProduct += U[i][j] * pi[j];
        }
        pi[i] = (0.0 - dotProduct) / U[i][i];
    }

    // normalize the solution vector
    double sum = 0.0;
    for (size_t i=0; i<num_states; i++)
    {
        sum += pi[i];
    }

    for (size_t i=0; i<num_states; i++)
    {
        pi[i] /= sum;
    }

    // return the stationary frequencies
    return pi;
}


/** This function checks that the rate matrix is time reversible. It takes as
 input the rate matrix, a, and the stationary frequencies of the process, f.
 It checks that f[i] * q[i][j] = f[j] * q[j][i] for all i != j. It does this
 by accumulating the difference | f[i] * q[i][j] - f[j] * q[j][i] | for all
 off-diagonal comparisons. If this difference is less than tolerance,
 it reports that the rate matrix is time-reversible. If the flag isRev
 is set to true, then we do not need to check because then we have determined
 previously that the rate matrix is reversible. */
bool RateMatrix_Rational::checkTimeReversibity( void )
{

    const std::vector<mpq_class>& theStationaryFreqs = stationary_frequencies;

    for (size_t i=0; i<num_states; i++)
    {

        for (size_t j=i+1; j<num_states; j++)
        {
//            if ( theStationaryFreqs[i] * (*the_rate_matrix)[i][j] != theStationaryFreqs[j] * (*the_rate_matrix)[j][i] )
            {
                return false;
            }
        }

    }

    return true;
}


RateMatrix_Rational* RateMatrix_Rational::clone( void ) const
{
    return new RateMatrix_Rational( *this );
}


std::vector<double> RateMatrix_Rational::getStationaryFrequencies(void) const
{
    std::vector<double> tmp = std::vector<double>( num_states, 0.0 );
    for (size_t i=0; i<num_states; ++i)
    {
        tmp[i] = stationary_frequencies[i].get_d();
    }
    
    return tmp;
}


double RateMatrix_Rational::getRate(size_t from, size_t to, double rate) const
{
    if ( from >= num_states || to > num_states )
    {
        throw RbException( "Index to RateMatrix.getRate() out of bounds" );
    }

//    return (*the_rate_matrix)[from][to].get_d() * rate;
    return -1;
}



double RateMatrix_Rational::getRate(size_t from, size_t to, double age, double rate) const
{
    if ( from >= num_states || to > num_states )
    {
        throw RbException( "Index to RateMatrix.getRate() out of bounds" );
    }

//    return (*the_rate_matrix)[from][to].get_d() * rate;
    return -1;
}

//MatrixRational RateMatrix_Rational::getRateMatrix() const
//{
//    return *the_rate_matrix;
//}



/** Rescale the rates such that the average rate is r */
void RateMatrix_Rational::rescaleToAverageRate(double r)
{

    double cur_ave = averageRate();
    double scale_factor = r / cur_ave;
    for (size_t i=0; i<num_states; i++)
    {
        for (size_t j=0; j<num_states; j++)
        {
//            (*the_rate_matrix)[i][j] *= scale_factor;
        }
    }

    // set flags
    needs_update = true;

}


/** Set the diagonal of the rate matrix such that each row sums to zero */
void RateMatrix_Rational::setDiagonal(void)
{

    for (size_t i=0; i<num_states; ++i)
    {
        double sum = 0.0;
        for (size_t j=0; j<num_states; ++j)
        {

            if (i != j)
            {
//                sum += (*the_rate_matrix)[i][j].get_d();
            }

        }
//        (*the_rate_matrix)[i][i] = -sum;
    }

    // set flags
    needs_update = true;
}

std::vector<int> RateMatrix_Rational::get_emitted_letters() const
{
    std::vector<int> emit(num_states);
    for(int i=0;i<num_states;i++)
        emit[i] = i;

    return emit;
}
