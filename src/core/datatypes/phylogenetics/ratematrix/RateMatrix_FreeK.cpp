#include <stddef.h>
#include <cmath>
#include <string>
#include <complex>
#include <iosfwd>
#include <vector>

#include "DistributionPoisson.h"
#include "EigenSystem.h"
#include "MatrixComplex.h"
#include "MatrixReal.h"
#include "RateMatrix_FreeK.h"
#include "RbException.h"
#include "RbMathLogic.h"
#include "RbMathMatrix.h"
#include "RbSettings.h"
#include "TransitionProbabilityMatrix.h"
#include "Assignable.h"
#include "GeneralRateMatrix.h"
#include "RbVector.h"
#include "RbVectorImpl.h"


using namespace RevBayesCore;

RateMatrix_FreeK::METHOD method_from_string(const std::string& method)
{
    // determine the type of matrix exponentiation
    if (method == "scalingAndSquaring")
    {
        return RateMatrix_FreeK::SCALING_AND_SQUARING;
    }
    else if (method == "scalingAndSquaringPade")
    {
        return RateMatrix_FreeK::SCALING_AND_SQUARING_PADE;
    }
    else if (method == "scalingAndSquaringTaylor")
    {
        return RateMatrix_FreeK::SCALING_AND_SQUARING_TAYLOR;
    }
    else if (method == "uniformization")
    {
        return RateMatrix_FreeK::UNIFORMIZATION;
    }
    else if (method == "eigen")
    {
        return RateMatrix_FreeK::EIGEN;
    }
    else
        return RateMatrix_FreeK::EIGEN; /// We haven't been complaining here, but maybe we should.
}

/** Construct rate matrix with n states */
RateMatrix_FreeK::RateMatrix_FreeK(size_t n) : RateMatrix_FreeK(n, true, EIGEN)
{ }


RateMatrix_FreeK::RateMatrix_FreeK(size_t n, bool r) : RateMatrix_FreeK(n, r, EIGEN)
{ }


RateMatrix_FreeK::RateMatrix_FreeK(size_t n, bool r, std::string method) : RateMatrix_FreeK( n, r, method_from_string(method))
{ }

RateMatrix_FreeK::RateMatrix_FreeK(size_t n, bool r, METHOD method) : GeneralRateMatrix( n ),
    rescale(r),
    my_method( method )
{
    // create the eigen system so the destructor has something to delete
    theEigenSystem       = new EigenSystem(the_rate_matrix);
    c_ijk.resize(num_states * num_states * num_states);
    cc_ijk.resize(num_states * num_states * num_states);
    
    matrixProducts = new std::vector<MatrixReal>();
    
    // Initialize emit_letters to [0...N-1]
    emit_letters.resize(num_states);
    for(int i=0;i<num_states;i++)
        emit_letters[i] = i;

    update();
}


/** Copy constructor */
RateMatrix_FreeK::RateMatrix_FreeK(const RateMatrix_FreeK& m) : GeneralRateMatrix( m )
{
    
    rescale               = m.rescale;
    my_method             = m.my_method;
    
    matrixProducts        = new std::vector<MatrixReal>( *m.matrixProducts );
    singleStepMatrix      = m.singleStepMatrix;
    maxRate               = m.maxRate;
    
    theEigenSystem        = new EigenSystem( *m.theEigenSystem );
    c_ijk                 = m.c_ijk;
    cc_ijk                = m.cc_ijk;

    theEigenSystem->setRateMatrixPtr(the_rate_matrix);
}


/** Destructor */
RateMatrix_FreeK::~RateMatrix_FreeK(void)
{
    
    delete theEigenSystem;
    delete matrixProducts;
}


RateMatrix_FreeK& RateMatrix_FreeK::operator=(const RateMatrix_FreeK &r)
{
    
    if (this != &r) 
    {
        GeneralRateMatrix::operator=( r );
       
        delete theEigenSystem;
        delete matrixProducts;
        
        rescale               = r.rescale;
        my_method             = r.my_method;
        
        matrixProducts        = new std::vector<MatrixReal>( *r.matrixProducts );
        singleStepMatrix      = r.singleStepMatrix;
        maxRate               = r.maxRate;
        
        theEigenSystem       = new EigenSystem( *r.theEigenSystem );
        c_ijk                = r.c_ijk;
        cc_ijk               = r.cc_ijk;
        
        theEigenSystem->setRateMatrixPtr(the_rate_matrix);

    }
    
    return *this;
}


RateMatrix_FreeK& RateMatrix_FreeK::assign(const Assignable &m)
{
    const RateMatrix_FreeK *rm = dynamic_cast<const RateMatrix_FreeK*>(&m);
    if ( rm != NULL )
    {
        return operator=(*rm);
    }
    else
    {
        throw RbException("Could not assign rate matrix.");
    }
    
}


/** Do precalculations on eigenvectors */
void RateMatrix_FreeK::calculateCijk(void)
{
    
    if ( theEigenSystem->isComplex() == false )
    {
        // real case
        const MatrixReal& ev  = theEigenSystem->getEigenvectors();
        const MatrixReal& iev = theEigenSystem->getInverseEigenvectors();
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
        const MatrixComplex& cev  = theEigenSystem->getComplexEigenvectors();
        const MatrixComplex& ciev = theEigenSystem->getComplexInverseEigenvectors();
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


/** Calculate the transition probabilities */
void RateMatrix_FreeK::calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const
{
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
        if ( theEigenSystem->isComplex() == false )
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


void RateMatrix_FreeK::checkMatrixIrreducible(double tolerance, TransitionProbabilityMatrix &P) const
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


void RateMatrix_FreeK::checkMatrixDiff(MatrixReal x, double tolerance, bool& diff) const
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


RateMatrix_FreeK* RateMatrix_FreeK::clone( void ) const
{
    return new RateMatrix_FreeK( *this );
}


void RateMatrix_FreeK::expandUniformization(int truncation, double tolerance) const
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


void RateMatrix_FreeK::expMatrixTaylor(MatrixReal &A, MatrixReal &F, double tolerance) const
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
    
    //now repeated squaring result s times
    for (size_t i = 0; i < s; i++)
    {
        F *= F;
    }
    
}


void RateMatrix_FreeK::exponentiateMatrixByScalingAndSquaring(double t,  TransitionProbabilityMatrix& p) const {
    
    // Here we use the scaling and squaring method with a 4th order Taylor approximant as described in:
    //
    // Moler, C., & Van Loan, C. 2003. Nineteen dubious ways to compute the exponential of a
    // matrix, twenty-five years later. SIAM review, 45(1), 3-49.
    //
    // I tested this implementation against the Eigen C++ package and a scaling parameter s = 6 had similar time
    // efficiency and returned the same results with about 10^-9 accuracy. The scaling parameter could be
    // increased for better accuracy.
    // -- Will Freyman 11/27/16
    size_t s = 6;
    
    // first scale the matrix
    double scale = t / pow(2, s);
    for ( size_t i = 0; i < num_states; i++ )
    {
        for ( size_t j = 0; j < num_states; j++ )
        {
            p[i][j] = (*the_rate_matrix)[i][j] * scale;
        }
    }
    
    // compute the 4th order Taylor approximant
    
    // calculate the scaled matrix raised to powers 2, 3 and 4
    TransitionProbabilityMatrix p_2(num_states);
    multiplyMatrices(p, p, p_2);
    
    TransitionProbabilityMatrix p_3(num_states);
    multiplyMatrices(p, p_2, p_3);
    
    TransitionProbabilityMatrix p_4(num_states);
    multiplyMatrices(p, p_3, p_4);
    
    // add k=0 (the identity matrix) and k=1 terms
    for ( size_t i = 0; i < num_states; i++ )
    {
        p[i][i] += 1;
    }
    
    // add the k=2, k=3, k=4 terms of the Taylor series
    for ( size_t i = 0; i < num_states; i++ )
    {
        for ( size_t j = 0; j < num_states; j++ )
        {
            p[i][j] += ( ( p_2[i][j] / 2 ) + ( p_3[i][j] / 6 ) + ( p_4[i][j] / 24 ) );
        }
    }
    
    // now perform the repeated squaring
    for (size_t i = 0; i < s; i++)
    {
        TransitionProbabilityMatrix r(num_states);
        multiplyMatrices(p, p, r);
        p = r;
    }
}


void RateMatrix_FreeK::fillRateMatrix( void )
{
    
    MatrixReal& m = *the_rate_matrix;
    
    // fill the rate matrix
    for (size_t i=0, k=0; i<num_states; i++)
    {
        double sum = 0.0;
        
        // off-diagonal
        for (size_t j=0; j<num_states; j++)
        {
            if (i==j)
            {
                continue;
            }
            double r = transition_rates[k];
            sum += r;
            m[i][j] = r;
            k++;
        }
        
        // diagonal
        m[i][i] = -sum;
    }
    
    // set flags
    needs_update = true;
}
    

inline void RateMatrix_FreeK::multiplyMatrices(TransitionProbabilityMatrix& p,  TransitionProbabilityMatrix& q,  TransitionProbabilityMatrix& r) const {
    
    // could probably use boost::ublas here, for the moment we do it ourselves.
    for ( size_t i = 0; i < num_states; i++ )
    {
        for ( size_t j = 0; j < num_states; j++ )
        {
            r[i][j] = 0;
            for ( size_t k = 0; k < num_states; k++ )
            {
                r[i][j] += p[i][k] * q[k][j];
            }
        }
    }
}


/** Calculate the transition probabilities for the real case */
void RateMatrix_FreeK::tiProbsEigens(double t, TransitionProbabilityMatrix& P) const
{
    
    // get a reference to the eigenvalues
    const std::vector<double>& eigenValue = theEigenSystem->getRealEigenvalues();
    
    // precalculate the product of the eigenvalue and the branch length
    std::vector<double> eigValExp(num_states);
    for (size_t s=0; s<num_states; s++)
    {
        eigValExp[s] = exp(eigenValue[s] * t);
    }
    
    // calculate the transition probabilities
    const double* ptr = &c_ijk[0];
    double*         p = P.theMatrix;
    for (size_t i=0; i<num_states; i++)
    {
        for (size_t j=0; j<num_states; j++, ++p)
        {
            double sum = 0.0;
            for (size_t s=0; s<num_states; s++)
            {
                sum += (*ptr++) * eigValExp[s];
            }
            
            //                  P[i][j] = (sum < 0.0) ? 0.0 : sum;
            (*p) = (sum < 0.0) ? 0.0 : sum;
        }
        
    }
    
//    double tol = RbSettings::userSettings().getTolerance();
//    checkMatrixIrreducible(tol, P);

}


/** Calculate the transition probabilities for the complex case */
void RateMatrix_FreeK::tiProbsComplexEigens(double t, TransitionProbabilityMatrix& P) const
{
    
    // get a reference to the eigenvalues
    const std::vector<double>& eigenValueReal = theEigenSystem->getRealEigenvalues();
    const std::vector<double>& eigenValueComp = theEigenSystem->getImagEigenvalues();
    
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
void RateMatrix_FreeK::tiProbsScalingAndSquaring(double t, TransitionProbabilityMatrix& P) const
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
        expMatrixTaylor(m, result, tol);
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
void RateMatrix_FreeK::tiProbsUniformization(double t, TransitionProbabilityMatrix& P) const
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
    
//    checkMatrixIrreducible(tol, P);
    
}


/** Update the eigen system */
void RateMatrix_FreeK::updateEigenSystem(void)
{
    
    theEigenSystem->update();
    calculateCijk();
    
}


void RateMatrix_FreeK::updateUniformization(void)
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


void RateMatrix_FreeK::update( void )
{
    
    if ( needs_update )
    {
        // assign all rate matrix elements
        fillRateMatrix();
        
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

void RateMatrix_FreeK::set_emitted_letters(const std::vector<int>& emit)
{
    assert(emit.size() == num_states);

    emit_letters = emit;
}

std::vector<int> RateMatrix_FreeK::get_emitted_letters() const
{
    return emit_letters;
}
