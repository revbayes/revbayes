#include "AbstractRateMatrix.h"

#include <cmath>
#include <fstream>
#include <algorithm>
#include <cstddef>

#include "EigenSystem.h"
#include "MatrixReal.h"
#include "RandomNumberGenerator.h"
#include "RandomNumberFactory.h"
#include "RbConstants.h"
#include "RbException.h"
#include "RbMathMatrix.h"
#include "RbMathLogic.h"
#include "RbSettings.h"
#include "DistributionPoisson.h"
#include "TransitionProbabilityMatrix.h"
#include "Cloneable.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "RlUserInterface.h"

using namespace RevBayesCore;

/** Construct rate matrix with n states */
AbstractRateMatrix::AbstractRateMatrix(size_t n, bool r, METHOD m) : RateMatrix(n),
    the_rate_matrix( new MatrixReal(num_states, num_states, 1.0) ),
    the_eigen_system( NULL ),
    needs_update( true ),
    rescale( r ),
    my_method( m )
{
    
    if ( my_method == EIGEN )
    {
        // create the eigen system so the destructor has something to delete
        the_eigen_system       = new EigenSystem(the_rate_matrix);
        c_ijk.resize(num_states * num_states * num_states);
        cc_ijk.resize(num_states * num_states * num_states);
    }
    
    matrixProducts = new std::vector<MatrixReal>();

    // I cannot call a pure virtual function from the constructor (Sebastian)
    //    update();
}



/** Copy constructor */
AbstractRateMatrix::AbstractRateMatrix(const AbstractRateMatrix& m) : RateMatrix(m),
    the_rate_matrix( new MatrixReal(*m.the_rate_matrix) ),
    the_eigen_system( NULL ),
    needs_update( true ),
    rescale( m.rescale )
{
    my_method               = m.my_method;
    
    if ( my_method == EIGEN )
    {
        the_eigen_system        = new EigenSystem( *m.the_eigen_system );
        c_ijk                   = m.c_ijk;
        cc_ijk                  = m.cc_ijk;
        
        the_eigen_system->setRateMatrixPtr(the_rate_matrix);
    }
    

}


/** Destructor */
AbstractRateMatrix::~AbstractRateMatrix(void)
{

    delete the_rate_matrix;
    delete the_eigen_system;
}


AbstractRateMatrix& AbstractRateMatrix::operator=(const AbstractRateMatrix &r)
{

    if (this != &r)
    {
        // delegate to parent class
        RateMatrix::operator=( r );

        delete the_rate_matrix;

        the_rate_matrix         = new MatrixReal( *r.the_rate_matrix );
        needs_update            = true;
        rescale                 = r.rescale;
        
        my_method               = r.my_method;
        
        if ( my_method == EIGEN )
        {
            the_eigen_system        = new EigenSystem( *r.the_eigen_system );
            c_ijk                   = r.c_ijk;
            cc_ijk                  = r.cc_ijk;
            
            the_eigen_system->setRateMatrixPtr(the_rate_matrix);
        }
        

    }

    return *this;
}


/** Do precalculations on eigenvectors */
void AbstractRateMatrix::calculateCijk(void)
{
    
    if ( the_eigen_system->isComplex() == false )
    {
        // real case
        const MatrixReal& ev  = the_eigen_system->getEigenvectors();
        const MatrixReal& iev = the_eigen_system->getInverseEigenvectors();
        double* pc = &c_ijk[0];
        for (size_t i=0; i<num_states; i++)
        {
            
            for (size_t j=0; j<num_states; j++)
            {
                
                for (size_t k=0; k<num_states; k++)
                {
                    *(pc++) = ev[i][k] * iev[k][j];
                }
                
            }
            
        }
        
    }
    else
    {
        // complex case
        const MatrixComplex& cev  = the_eigen_system->getComplexEigenvectors();
        const MatrixComplex& ciev = the_eigen_system->getComplexInverseEigenvectors();
        std::complex<double>* pc = &cc_ijk[0];
        for (size_t i=0; i<num_states; i++)
        {
            
            for (size_t j=0; j<num_states; j++)
            {
                
                for (size_t k=0; k<num_states; k++)
                {
                    *(pc++) = cev[i][k] * ciev[k][j];
                }
                
            }
            
        }
        
    }
    
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
std::vector<double> AbstractRateMatrix::calculateStationaryFrequencies(void) const
{

    // transpose the rate matrix and put into QT
    MatrixReal QT(num_states, num_states);
    for (size_t i=0; i<num_states; i++)
    {
        for (size_t j=0; j<num_states; j++)
        {
            QT[i][j] = (*the_rate_matrix)[j][i];
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


/** Calculate the transition probabilities */
void AbstractRateMatrix::calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const
{
    assert(num_states == P.num_states);

    // The eigensystem code was returning NaN likelihood values when transition rates
    // were close to 0.0, so now we use the scaling and squaring method.
    double t = rate * (startAge - endAge);
    if (my_method == SCALING_AND_SQUARING)
    {
        exponentiateMatrixByScalingAndSquaring(t, P);
    }
    else if (my_method == SCALING_AND_SQUARING_PADE || my_method == SCALING_AND_SQUARING_TAYLOR)
    {
        tiProbsScalingAndSquaring(t, P);
    }
    else if (my_method == UNIFORMIZATION)
    {
        tiProbsUniformization(t, P);
    }
    else if (my_method == EIGEN)
    {
        if ( the_eigen_system->isComplex() == false )
        {
            tiProbsEigens(t, P);
        }
        else
        {
            tiProbsComplexEigens(t, P);
        }
    }

    for(int i=0; i < num_states; i++)
    {
        double total = 0;
        for(int j=0; j < num_states; j++)
        {
#ifdef DEBUG_MATRIX_EXPONENTIAL
            if (P[i][j] < -1.0e-6)
                std::cerr<<my_method<<" exp(M)["<<i<<","<<j<<"] = "<<P[i][j]<<"\n";
#endif
            P[i][j] = std::max(0.0, P[i][j]);
            total += P[i][j];
        }
#ifdef DEBUG_MATRIX_EXPONENTIAL
        if (std::abs(total - 1.0) > num_states * 1.0e-6)
            std::cerr<<my_method<<" exp(M)["<<i<<"]: row sum = "<<total<<"\n";
#endif
        double scale = 1.0/total;
        for(int j=0; j < num_states; j++)
            P[i][j] *= scale;
    }
}


/** This function computes the transition probability matrix that is used for
 stochastic mapping. In general, we use the matrix uniformization method
 for stochastic mapping, which requires P(t) = exp(Qt). This function may be
 overridden if an alternative method is needed to compute P(t) for the model.
 See RateMatrix_DECRateMatrix for an example.
 */
void AbstractRateMatrix::calculateTransitionProbabilitiesForStochasticMapping(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const
{
    calculateTransitionProbabilities(startAge, endAge, rate, P);
}


void AbstractRateMatrix::computeDominatingRate(void) const
{
    dominating_rate = 0.0;
    for (size_t i = 0; i < num_states; i++)
    {
        double r = -(*the_rate_matrix)[i][i];
        if ( r > dominating_rate)
        {
            dominating_rate = r;
        }
    }
}


void AbstractRateMatrix::checkMatrixDiff(MatrixReal x, double tolerance, bool& diff) const
{
    for (size_t i = 0; i < num_states; ++i)
    {
        for (size_t j = 0; j < num_states; ++j)
        {
            diff = (std::abs(x[i][j]) < tolerance) ? true : false;
            if (diff == false)
            {
                return;
            }
        }
    }
}


void AbstractRateMatrix::checkMatrixIrreducible(double tolerance, TransitionProbabilityMatrix &P) const
{
    // check if the Q matrix is irreducible by checking if there is any element in the P matrix
    // that is smaller than some specified tolerance
    // and if that's the case, fill in the P matrix with all zeros
    // so that the current proposal will certainly get rejected
    // here we assume that all the states in the Q matrix exist in the observed data
    bool irreducible = true;
    
    for (size_t i = 0; i < num_states - 1; ++i)
    {
        for (size_t j = i + 1; j < num_states; ++j)
        {
            if (P[i][j] < tolerance && P[j][i] < tolerance)
            {
                irreducible = false;
                break;
            }
        }
        if (irreducible == false)
        {
            break;
        }
    }
    
    if (irreducible == false)
    {
        for (size_t i = 0; i < num_states; ++i)
        {
            P[i][i] = 1.0;
            for (size_t j = i + 1; j < num_states; ++j)
            {
                P[i][j] = 0.0;
                P[j][i] = 0.0;
            }
        }
    }
}


/** This function checks that the rate matrix is time reversible. It takes as
 input the rate matrix, a, and the stationary frequencies of the process, f.
 It checks that f[i] * q[i][j] = f[j] * q[j][i] for all i != j. It does this
 by accumulating the difference | f[i] * q[i][j] - f[j] * q[j][i] | for all
 off-diagonal comparisons. If this difference is less than tolerance,
 it reports that the rate matrix is time-reversible. If the flag isRev
 is set to true, then we do not need to check because then we have determined
 previously that the rate matrix is reversible. */
bool AbstractRateMatrix::checkTimeReversibity(double tolerance)
{

    std::vector<double> theStationaryFreqs = getStationaryFrequencies();
	double diff = 0.0;
	for (size_t i=0; i<num_states; i++)
    {

        for (size_t j=i+1; j<num_states; j++)
        {
			diff += fabs( theStationaryFreqs[i] * (*the_rate_matrix)[i][j] - theStationaryFreqs[j] * (*the_rate_matrix)[j][i] );
        }

    }
    //    reversibilityChecked = true;
	if (diff < tolerance)
    {
        return true;
    }

	return false;
}


void AbstractRateMatrix::expandUniformization(int truncation, double tolerance) const
{

    int n = (int)matrixProducts->size();
    int d = truncation - n;
    int i = 0;
    
    for (; i < d; ++i)
    {
        // add terms of the power series to matrix products until the difference between the last two terms is smaller than the tolerance
        MatrixReal diffMatrix = matrixProducts->at(n - 1 + i) - matrixProducts->at(n - 2 + i);
        
        bool diff = true;
        checkMatrixDiff(diffMatrix, tolerance, diff);
        if (diff == true)
        {
            break;
        }

        MatrixReal m = singleStepMatrix * matrixProducts->at(n - 1 + i);
        matrixProducts->push_back(m);
    }
    
    // if the current size of the matrix products is still smaller than the truncation, fill all the remaining terms with the same converged matrix
    MatrixReal m = matrixProducts->at(n - 1 + i);
    for (int j = i; j < d; ++j)
    {
        matrixProducts->push_back(m);
    }
    
}

void AbstractRateMatrix::exponentiateMatrixByScalingAndSquaring(double t, TransitionProbabilityMatrix& p) const
{
    
    MatrixReal m = *the_rate_matrix;
    m *= t;
    
    MatrixReal result(num_states);
    double tol = RbSettings::userSettings().getTolerance();
    
    
    exponentiateMatrixTaylor(m, result, tol);
    
//    assert(std::isnan(t) or t >= 0);
//    assert(p.num_states == p.getNumberOfStates());
//
//    // Here we use the scaling and squaring method with a 4th order Taylor approximant as described in:
//    //
//    // Moler, C., & Van Loan, C. 2003. Nineteen dubious ways to compute the exponential of a
//    // matrix, twenty-five years later. SIAM review, 45(1), 3-49.
//    //
//    // I tested this implementation against the Eigen C++ package and a scaling parameter s = 6 had similar time
//    // efficiency and returned the same results with about 10^-9 accuracy. The scaling parameter could be
//    // increased for better accuracy.
//    // -- Will Freyman 11/27/16
//
//    // Note that Will had previously chosen a fixed value of 8 scalings/squarings, but this
//    // didn't work for large branch lengths.
//
//    // We need the norm to be small enough for the 4-th order Taylor series to be a good approximation.
//    // This code gives the number of scalings needed to bring the norm down below 2^-10 = 0.000977
//    // This is a heuristic choice, but note that that ((2^-10)^4)/24 = 3.79e-14
//    // The fifth order term would be ((2^-10)^5)/120 = 7.4e-18
//    // Compare 10 with the number 12 in RateMatrix_FreeK::expMatrixTaylor( ) in RateMatrix_FreeK.cpp
//    // -- Ben Redelings 10/23/22
//    double norm = L1norm(*the_rate_matrix) * std::abs(t);
//    int s = 0;
//    std::frexp(norm, &s);
//    s = std::max(10 + s, 0);
//
//    // first scale the matrix
//    double scale = t / pow(2, s);
//    for ( size_t i = 0; i < num_states; i++ )
//    {
//        for ( size_t j = 0; j < num_states; j++ )
//        {
//            p[i][j] = (*the_rate_matrix)[i][j] * scale;
//        }
//    }
//
//    // compute the 4th order Taylor approximant
//
//    // BDR: Note more recent work uses Pade approximants, which I think are better.
//    // Eigen implements that approach, which goes up to a (13,13) Pade approximant, depending
//    //   on whether the extra work is needed.
//    // Its kind of complicated, so we should probably use Eigen if we switch to that.
//
//    // calculate the scaled matrix raised to powers 2, 3 and 4
//    TransitionProbabilityMatrix p_2 = p * p;
//
//    TransitionProbabilityMatrix p_3 = p * p_2;
//
//    TransitionProbabilityMatrix p_4 = p_2 * p_2;
//
//    // add k=0 (the identity matrix) and k=1 terms
//    for ( size_t i = 0; i < num_states; i++ )
//    {
//        p[i][i] += 1;
//    }
//
//    // add the k=2, k=3, k=4 terms of the Taylor series
//    for ( size_t i = 0; i < num_states; i++ )
//    {
//        for ( size_t j = 0; j < num_states; j++ )
//        {
//            p[i][j] += ( ( p_2[i][j] / 2 ) + ( p_3[i][j] / 6 ) + ( p_4[i][j] / 24 ) );
//        }
//    }
//
//    // Make sure that our Taylor approximation is now a stochastic matrix
//    // BEFORE we start squaring it.
//    ensure_nonnegative(p);
//    normalize_rows(p);
//
//    // now perform the repeated squaring
//    TransitionProbabilityMatrix r(num_states);
//    for (size_t i = 0; i < s; i++)
//    {
//        // We could do p = p * p, but that allocates memory.
//        p.multiplyTo(p, r);
//        p = std::move(r);
//
//        // Handle roundoff-error.
//        normalize_rows(p);
//    }
}


void AbstractRateMatrix::exponentiateMatrixTaylor(MatrixReal &A, MatrixReal &F, double tolerance) const
{
    // here use the global tolerance to determine the truncation order for Taylor series
    // alternatively, it can be determined according to Table 1 in Moler and Van Loan, 2003
    
    // here we dynamically determine the scaling parameter according to the norm of Q*t
    // first compute the norm of Q*t
    double normA = 0.0;
    for (size_t i = 0; i < num_states; ++i)
    {
        double x = std::abs(A[i][i]);
        if (x > normA)
            normA = x;
    }
    normA *= 2.0;
    
    // check if the Q matrix is irreducible
    // if that is not the case, directly exit the current function
    // as otherwise the following taylor series or repeated squaring loops may not be finite
    // here we assume that all the states in the Q matrix exist in the observed data
    if (RbMath::isNan(normA) || normA == 0.0)
    {
        return;
    }
    
    // then use the norm to determine the appropriate scaling parameter s
    // here plus 4 is a decision made under the guidance of Table 1 in Moler and Van Loan, 2003 and my experiments
    // the original implementation in R expm package for Pade approximation was s = max(0, [log2(norm)]+1), which seems not enough for Taylor series
    // and it turns out the intermediate norm values are those that tends to be the most problematic ones under the formula
    // but they actually would be the most common ones for the usage in phylogenetics (branch length around 0.1)
    // so I decide to jack it up by 11 here, which results in a transition probability matrix
    // whose sum differ from the ones computed under the eigen method by about 1e-9 under 0.1 branch length
    // if that turns out to be still insufficient, a larger number should be considered
    // Jiansi Gao 09/07/2017
    int e = ceil(log2(normA)) + 12;
    int s = (e < 0) ? 0 : e;
    
    // scale the matrix by 2^s
    double scale = pow(2, s);
    A *= 1.0/scale;
    
    // the first term of the taylor series is the identity matrix
    for (size_t i = 0; i < num_states; ++i)
    {
        F[i][i] = 1.0;
    }
    
    // add remaining terms of taylor series until a given term is smaller than the tolerance
    int it = 1;
    int fact = 1;
    while (true)
    {
        MatrixReal m = A * (1.0/fact);
        F += m;
        
        bool diff = true;
        checkMatrixDiff(m, tolerance, diff);
        if (diff == true)
        {
            break;
        }
        
        ++it;
        fact *= it;
        A *= A;
    }
    
    // now repeated squaring result s times
    for (size_t i = 0; i < s; i++)
    {
        F *= F;
    }
    
}


void AbstractRateMatrix::fillRateMatrix( void )
{
    
    MatrixReal& m = *the_rate_matrix;
    
//    // fill the rate matrix
//    for (size_t i=0, k=0; i<num_states; i++)
//    {
//        double sum = 0.0;
//        
//        // off-diagonal
//        for (size_t j=0; j<num_states; j++)
//        {
//            if (i==j)
//            {
//                continue;
//            }
//            double r = transition_rates[k];
//            sum += r;
//            m[i][j] = r;
//            k++;
//        }
//        
//        // diagonal
//        m[i][i] = -sum;
//    }
    
    // set flags
    needs_update = true;
}



double AbstractRateMatrix::getDominatingRate(void) const
{
    return dominating_rate;
}


std::vector<int> AbstractRateMatrix::get_emitted_letters( void ) const
{
    std::vector<int> emit(num_states);
    for (int i=0;i<num_states;i++)
    {
        emit[i] = i;
    }
    
    return emit;
}


double AbstractRateMatrix::getRate(size_t from, size_t to, double rate) const
{
    if ( from >= num_states || to > num_states )
    {
        throw RbException( "Index to RateMatrix.getRate() out of bounds" );
    }

    return (*the_rate_matrix)[from][to] * rate;
}



double AbstractRateMatrix::getRate(size_t from, size_t to, double age, double rate) const
{
    if ( from >= num_states || to > num_states )
    {
        throw RbException( "Index to RateMatrix.getRate() out of bounds" );
    }

    return (*the_rate_matrix)[from][to] * rate;
}

MatrixReal AbstractRateMatrix::getRateMatrix() const
{
    return *the_rate_matrix;
}


std::vector<double> AbstractRateMatrix::getStationaryFrequencies( void ) const
{
    return calculateStationaryFrequencies();
}

MatrixReal AbstractRateMatrix::getStochasticMatrix(size_t n) const
{
    if (n >= stochastic_matrix.size())
    {
        std::stringstream ss;
        ss << "Cannot access stochastic matrix R^" << n;
        throw RbException(ss.str());
    }
    return stochastic_matrix[n];
}

void AbstractRateMatrix::computeStochasticMatrix(size_t n) const
{
    MatrixReal r;
    if (n == 0) {
        // identity matrix, R^0
        r = MatrixReal(num_states, num_states);
        for (size_t i = 0; i < num_states; i++) {
            r[i][i] = 1.0;
        }
    }
    else if (n == 1) {
        // stochastic matrix, R^1
        r = (*the_rate_matrix) * (1.0/dominating_rate) + stochastic_matrix[0];
    }
    else {
        // stochastic matrix, R^n = R^(n-1) * R^1
        MatrixReal r_n_minus_1 = stochastic_matrix[n-1];

        // helps manage machine precision error/underflow for large R^n
        double smallest_non_zero = 1.0; //RbConstants::Double::inf;
        for (size_t i = 0; i < this->num_states; i++)
        {
            if (r_n_minus_1[i][i] < smallest_non_zero && r_n_minus_1[i][i] > 0.0)
            {
                smallest_non_zero = r_n_minus_1[i][i];
            }
        }
        r_n_minus_1 *= (1.0 / smallest_non_zero);

//        r = stochastic_matrix[n-1] * stochastic_matrix[1];
        r = r_n_minus_1 * stochastic_matrix[1];
        r *= smallest_non_zero;

    }

    if (stochastic_matrix.size() <= n) {
        stochastic_matrix.push_back(r);
    }
}

/** Rescale the rates such that the average rate is r */
void AbstractRateMatrix::rescaleToAverageRate(double r)
{

    double cur_ave = averageRate();
    double scale_factor = r / cur_ave;
    for (size_t i=0; i<num_states; i++)
    {
        for (size_t j=0; j<num_states; j++)
        {
            (*the_rate_matrix)[i][j] *= scale_factor;
        }
    }

    // set flags
    needs_update = true;

}


/** Set the diagonal of the rate matrix such that each row sums to zero */
void AbstractRateMatrix::setDiagonal(void)
{

    for (size_t i=0; i<num_states; ++i)
    {
        double sum = 0.0;
        for (size_t j=0; j<num_states; ++j)
        {

            if (i != j)
            {
                sum += (*the_rate_matrix)[i][j];
            }

        }
        (*the_rate_matrix)[i][i] = -sum;
    }

    // set flags
    needs_update = true;
}


bool AbstractRateMatrix::simulateStochasticMapping(double startAge, double endAge, double rate,std::vector<size_t>& transition_states, std::vector<double>& transition_times) const
{
    // start and end states
    size_t start_state = transition_states[0];
    size_t end_state = transition_states[1];
    double branch_length = (startAge - endAge);

    if (branch_length == 0.0) {
        transition_states = std::vector<size_t>(1, start_state);
        transition_times = std::vector<double>(1, 0.0);
        return true;
    }

    // transition probabilities
    TransitionProbabilityMatrix P(num_states);
    calculateTransitionProbabilitiesForStochasticMapping(startAge, endAge, rate, P);
//    exponentiateMatrixByScalingAndSquaring(branch_length * rate, P);
    stochastic_matrix = std::vector<MatrixReal>();

    // dominating rate
    computeDominatingRate();

    // sample number of events
    size_t num_events = 0;
    double lambda = branch_length * rate * dominating_rate;
    double prob_transition_ctmc = P[start_state][end_state];

//    if (prob_transition_ctmc < 1e-20) {
//        std::cout << P << "\n";
//    }
    double u = GLOBAL_RNG->uniform01() * (1.0 - 1e-3);
    double g = u * prob_transition_ctmc;
    double prob_num_events_sum = 0.0;
    double prob_total_sum = 0.0;

    while (g > 0.0) {

        // probability for num_events
        double prob_num_events = RbStatistics::Poisson::pdf(lambda, (int)num_events);
        prob_num_events_sum += prob_num_events;

        // add the R^n stochastic matrix
        computeStochasticMatrix(num_events);

        // probability of start_state -> end_state after num_events
        const MatrixReal& R_n = getStochasticMatrix(num_events);
        double prob_transition_dtmc = R_n[start_state][end_state];

        // update sampling prob
        double prob_total = prob_num_events * prob_transition_dtmc;

        prob_total_sum += prob_total;

        /*
        // Debugging statements to be removed
        std::cout << "num_events " << num_events << "\n";
        std::cout << start_state << " -> " << end_state << "\n";
        std::cout << "P(n) = " << prob_num_events << "\n";
        std::cout << "P(b|a,n) = " << prob_transition_dtmc << "\n";
        std::cout << "P(b,n|a) = " << prob_total << "\n";
        std::cout << "sum_n P(b|a) = " << prob_total_sum << "\n";
        std::cout << "P(b|a) = " << prob_transition_ctmc << "\n";
        std::cout << "brlen " << branch_length << "\n";
        std::cout << "rate " << rate << "\n";
        std::cout << "dom_rate " << dominating_rate << "\n";
        std::cout << "lambda " << lambda << "\n";
        std::cout << "u " << u << "\n";
        std::cout << "log P(b|a) = " << log(prob_transition_ctmc) << "\n";
        std::cout << "P\n";
        std::cout << P << "\n";

        if (num_events > 0) {
            std::cout << "R^1\n";
            std::cout << getStochasticMatrix(1) << "\n";
        }
        std::cout << "R^" << num_events << "\n";
        std::cout << R_n << "\n";

        std::cout << "Q\n";
        std::cout << *the_rate_matrix << "\n";
        std::cout << "\n";
        */

        g -= prob_total;
        if (g <= 0.0)
        {
            break;
        }

        // make sure we simulate plenty of events
        size_t max_events = 20 + 5 * lambda;
        if (num_events > max_events) {
            return false;
            // throw RbException("AbstractRateMatrix::simulateStochasticMapping fails to converge when sampling num_events!");
        }

        num_events += 1;

    }

    // sample event types per interval
    transition_states = std::vector<size_t>(1, start_state);
    for (size_t n = 0; n < num_events; n++)
    {
        size_t prev_state = transition_states[n];
        size_t num_events_left = num_events - n - 1;

        const MatrixReal& R_1 = getStochasticMatrix(1);
        const MatrixReal& R_n = getStochasticMatrix(num_events_left);

        // get the normalization constant for sampling
        double p_sum = 0.0;
        for (size_t j = 0; j < num_states; j++)
        {
            p_sum += R_1[prev_state][j] * R_n[j][end_state];
        }

        // sample transition to next state
        size_t next_state = 0;
        double u = GLOBAL_RNG->uniform01() * p_sum;
        while (u > 0) {

            double p = R_1[prev_state][next_state] * R_n[next_state][end_state];
            u -= p;
            if (u <= 0) {
                break;
            }
            next_state += 1;
        }

        // sample the next state for the n+1 event
        transition_states.push_back(next_state);

    }

    // sample event times
    for (size_t i = 0; i < num_events; i++)
    {
        transition_times.push_back(GLOBAL_RNG->uniform01() * branch_length);
    }
    transition_times.push_back(0.0);
    std::sort( transition_times.begin(), transition_times.end() );


    // filter out the virtual events
    std::vector<size_t> save_states(1, transition_states[0]);
    std::vector<double> save_times(1, 0.0);
    size_t prev_state = transition_states[0];
    for (size_t i = 1; i < transition_states.size(); i++) {
        if (prev_state != transition_states[i]) {
            save_states.push_back(transition_states[i]);
            save_times.push_back(transition_times[i]);
            prev_state = transition_states[i];
        }
    }
    transition_states = save_states;
    transition_times = save_times;

    // convert relative times along branches into incremental times
    double sum_transition_times = 0.0;
    for (size_t i = 0; i < transition_times.size()-1; i++)
    {
        transition_times[i] = transition_times[i+1] - sum_transition_times;
        sum_transition_times += transition_times[i];
    }
    transition_times[transition_times.size()-1] = branch_length - sum_transition_times;

    // done!
    return true;
}


/** Calculate the transition probabilities for the real case */
void AbstractRateMatrix::tiProbsEigens(double t, TransitionProbabilityMatrix& P) const
{
    
    // get a reference to the eigenvalues
    const std::vector<double>& eigenValue = the_eigen_system->getRealEigenvalues();
    
    // precalculate the product of the eigenvalue and the branch length
    std::vector<double> eigValExp(num_states);
    for (size_t s=0; s<num_states; s++)
    {
        eigValExp[s] = exp(eigenValue[s] * t);
    }
    
    // calculate the transition probabilities
    const double* ptr = &c_ijk[0];
    double*         p = P.theMatrix;
    for (int i=0; i<num_states; i++)
    {
        for (int j=0; j<num_states; j++)
        {
            double sum = 0.0;
            for (int s=0; s<num_states; s++)
                sum += (*ptr++) * eigValExp[s];
            
            p[j] = (sum < 0.0) ? 0.0 : sum;
        }

        p += num_states;
    }
    
//    double tol = RbSettings::userSettings().getTolerance();
//    checkMatrixIrreducible(tol, P);

}


/** Calculate the transition probabilities for the complex case */
void AbstractRateMatrix::tiProbsComplexEigens(double t, TransitionProbabilityMatrix& P) const
{
    
    // get a reference to the eigenvalues
    const std::vector<double>& eigenValueReal = the_eigen_system->getRealEigenvalues();
    const std::vector<double>& eigenValueComp = the_eigen_system->getImagEigenvalues();
    
    // precalculate the product of the eigenvalue and the branch length
    std::vector<std::complex<double> > ceigValExp(num_states);
    for (size_t s=0; s<num_states; s++)
    {
        std::complex<double> ev = std::complex<double>(eigenValueReal[s], eigenValueComp[s]);
        ceigValExp[s] = exp(ev * t);
    }
    
    // calculate the transition probabilities
    const std::complex<double>* ptr = &cc_ijk[0];
    for (size_t i=0; i<num_states; i++)
    {
        for (size_t j=0; j<num_states; j++)
        {
            std::complex<double> sum = std::complex<double>(0.0, 0.0);
            for (size_t s=0; s<num_states; s++)
            {
                sum += (*ptr++) * ceigValExp[s];
            }

            P[i][j] = (sum.real() < 0.0) ? 0.0 : sum.real();
        }
    }
    
//    double tol = RbSettings::userSettings().getTolerance();
//    checkMatrixIrreducible(tol, P);
    
}


/** Calculate the transition probabilities with scaling and squaring */
void AbstractRateMatrix::tiProbsScalingAndSquaring(double t, TransitionProbabilityMatrix& P) const
{
    // Moler, C., & Van Loan, C. 2003. Nineteen dubious ways to compute the exponential of a
    // matrix, twenty-five years later. SIAM review, 45(1), 3-49.
    // Higham, N. 2009. The Scaling and Squaring Method for the Matrix Exponential Revisited.
    // SIAM review, 51(4), 747-764.

    MatrixReal m = *the_rate_matrix;
    m *= t;
    MatrixReal result(num_states);
    double tol = RbSettings::userSettings().getTolerance();
    
    if (my_method == SCALING_AND_SQUARING_PADE)
    {
        // the value of truncation computed by findPadeQValue is 5 under RevBayes default tolerance (1e-9)
        // which seems a bit too generous comparing with the value given in Table 1 of Moler and Van Loan, 2003
        // here we use max(4, findPadeQValue(tol)) as the truncation
        // since findPadeQValue will return 4 if the tolerance is 1e-8 and 4 is also recommended in Table 1
        // if that turns out to be insufficient or if a higher accuracy is desired, a larger number should be considered
        // Jiansi Gao 09/07/2017
        int truncation = RbMath::findPadeQValue(tol);
        if (truncation > 4)
        {
            truncation = 4;
        }
        RbMath::expMatrixPade(m, result, truncation);
    }
    else if (my_method == SCALING_AND_SQUARING_TAYLOR)
    {
        exponentiateMatrixTaylor(m, result, tol);
    }
    
    // fill in P from result
    for (size_t i = 0; i < num_states; ++i)
    {
        for (size_t j = 0; j < num_states; ++j)
        {
            P[i][j] = (result[i][j] < 0.0) ? 0.0 : result[i][j];
        }
    }

// Probably we should not make the assumption that the Q matrix should be irreducible as the probability could still be positive
// under reducible Q matrix if not all the states of the Q matrix appear at the tip, so it probably makes more sense just let
// the likelihood take care of the reducibility, despite of the somewhat tiny possibility of numerical instability.
// Jiansi Gao 06/22/2018
//    checkMatrixIrreducible(tol, P);

}


/** Calculate the transition probabilities with uniformization */
void AbstractRateMatrix::tiProbsUniformization(double t, TransitionProbabilityMatrix& P) const
{
    
    // compute the appropriate truncation given t
    // the formula is taken from Tataru and Hobolth, 2011, BMC Bioinformatics
    // which seems to be pretty generous in most cases, so it should be sufficient for now
    // if not, a larger number should be considered
    // Jiansi Gao 09/07/2017
    MatrixReal result(num_states);
    double tol = RbSettings::userSettings().getTolerance();
    
    // check if the Q matrix is irreducible
    // if that is not the case, directly fill in the P matrix with all zeros
    // as otherwise the following loop which expands the power series of the uniformized matrix may not be finite
    // here we assume that all the states in the Q matrix exist in the observed data
    if ((RbMath::isNan(maxRate) || maxRate >= 0.0) == false)
    {
        double lambda = -maxRate * t;
        int truncation = std::ceil(4 + 6 * sqrt(lambda) + lambda);
        
        // update the matrix poisson products sum vector if necessary
        expandUniformization(truncation, tol);
        
        // compute the transition probability by weighted average
        for (size_t i = 0; i < truncation; ++i)
        {
            
            // compute the poisson probability
            double p = RbStatistics::Poisson::pdf(lambda, (int)i);
            
            // add the weighted terms of taylor series until a given term is smaller than the tolerance
            result += matrixProducts->at(i) * p;
            
        }
    }
    
    // fill in P from result
    for (size_t i = 0; i < num_states; ++i)
    {
        for (size_t j = 0; j < num_states; ++j)
        {
            P[i][j] = (result[i][j] < 0.0) ? 0.0 : result[i][j];
        }
    }
    
}


void AbstractRateMatrix::update( void )
{
    
    if ( needs_update )
    {
        // assign all rate matrix elements
//        fillRateMatrix();
        updateInternalRateMatrix();
        
        // rescale
        if ( rescale == true )
        {
            rescaleToAverageRate( 1.0 );
        }

        // update the uniformization system if necessary
        if (my_method == UNIFORMIZATION)
        {
            updateUniformization();
        }
        // update the eigensystem if necessary
        if (my_method == EIGEN)
        {
            updateEigenSystem();
        }
        
        // clean flags
        needs_update = false;
    }
    
}

void AbstractRateMatrix::updateInternalRateMatrix( void )
{
    // empty
    // for derived classes to overwrite if wanted.
}


/** Update the eigen system */
void AbstractRateMatrix::updateEigenSystem(void)
{
    
    the_eigen_system->update();
    calculateCijk();
    
}


void AbstractRateMatrix::updateUniformization(void)
{
    
    // we need to fill in the single-step transition probability matrix
    
    // find the diagonial element of the matrix with the maximal value
    MatrixReal m = *the_rate_matrix;
    maxRate = m[0][0];
    for (size_t i = 1; i < num_states; ++i)
    {
        if (m[i][i] < maxRate )
        {
            maxRate = m[i][i];
        }
    }
    
    // for the given max rate, fill in the single-step transition probability matrix
    singleStepMatrix = MatrixReal(num_states);
    for (size_t i = 0; i < num_states; ++i)
    {
        singleStepMatrix[i][i] = 1 - m[i][i] / maxRate;
        for (size_t j = i + 1; j < num_states; ++j)
        {
            singleStepMatrix[i][j] = -m[i][j] / maxRate;
            singleStepMatrix[j][i] = -m[j][i] / maxRate;
        }
    }
    
    // finally, clear the vector of matrix products
    matrixProducts->clear();
    
    // add the identity matrix (the first one) and the singleStepMatrix (the second one)
    MatrixReal identity_matrix(num_states);
    for (size_t i = 0; i < num_states; ++i)
    {
        identity_matrix[i][i] = 1.0;
    }
    matrixProducts->push_back(identity_matrix);
    matrixProducts->push_back(singleStepMatrix);

}






double L1norm(const MatrixReal& M)
{
    // ||A||_1 = sup x ||A x||_1 / ||x||_1
    //         = max 1 <= j<= n \sum_i |A_{i,j}|
    //         = max over columns of (the sum of the absolute values in each column).
    double norm = 0;
    for(int j=0; j < M.getNumberOfColumns(); j++)
    {
        double sum = 0;
        for(int i=0; i < M.getNumberOfRows(); i++)
            sum += std::abs(M[i][j]);
        norm = std::max(norm, sum);
    }
    assert(norm >= 0);
    return norm;
}

double L1norm(const TransitionProbabilityMatrix& M)
{
    // ||A||_1 = sup x ||A x||_1 / ||x||_1
    //         = max 1 <= j<= n \sum_i |A_{i,j}|
    //         = max over columns of (the sum of the absolute values in each column).
    double norm = 0;
    for(int j=0; j < M.getNumberOfStates(); j++)
    {
        double sum = 0;
        for(int i=0; i < M.getNumberOfStates(); i++)
            sum += std::abs(M[i][j]);
        norm = std::max(norm, sum);
    }
    assert(norm >= 0);
    return norm;
}
