#include <cstdlib>
#include <cmath>
#include <string>
#include <complex>
#include <iosfwd>
#include <vector>

#include "EigenSystem.h"
#include "MatrixComplex.h"
#include "MatrixReal.h"
#include "RateMatrix_GTR.h"
#include "RbException.h"
#include "TransitionProbabilityMatrix.h"
#include "Assignable.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "StringUtilities.h"
#include "TimeReversibleRateMatrix.h"

using namespace RevBayesCore;

/** Construct rate matrix with n states */
RateMatrix_GTR::RateMatrix_GTR(size_t n) : TimeReversibleRateMatrix( n )
{
    
    theEigenSystem       = new EigenSystem(the_rate_matrix);
    c_ijk.resize(num_states * num_states * num_states);
    cc_ijk.resize(num_states * num_states * num_states);
    
    update();
}


/** Copy constructor */
RateMatrix_GTR::RateMatrix_GTR(const RateMatrix_GTR& m) : TimeReversibleRateMatrix( m )
{
    
    theEigenSystem       = new EigenSystem( *m.theEigenSystem );
    c_ijk                = m.c_ijk;
    cc_ijk               = m.cc_ijk;
    
    theEigenSystem->setRateMatrixPtr(the_rate_matrix);
}


/** Destructor */
RateMatrix_GTR::~RateMatrix_GTR(void)
{
    
    delete theEigenSystem;
}


RateMatrix_GTR& RateMatrix_GTR::operator=(const RateMatrix_GTR &r)
{
    
    if (this != &r)
    {
        TimeReversibleRateMatrix::operator=( r );
        
        delete theEigenSystem;
        
        theEigenSystem       = new EigenSystem( *r.theEigenSystem );
        c_ijk                = r.c_ijk;
        cc_ijk               = r.cc_ijk;
        
        theEigenSystem->setRateMatrixPtr(the_rate_matrix);
    }
    
    return *this;
}


RateMatrix_GTR& RateMatrix_GTR::assign(const Assignable &m)
{
    
    const RateMatrix_GTR *rm = dynamic_cast<const RateMatrix_GTR*>(&m);
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
void RateMatrix_GTR::calculateCijk(void)
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
void RateMatrix_GTR::calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const
{
    double t = rate * (startAge - endAge);
    if ( theEigenSystem->isComplex() == false )
    {
        tiProbsEigens(t, P);
    }
    else
    {
        tiProbsComplexEigens(t, P);
    }
}


RateMatrix_GTR* RateMatrix_GTR::clone( void ) const
{
    return new RateMatrix_GTR( *this );
}


/** Calculate the transition probabilities for the real case */
void RateMatrix_GTR::tiProbsEigens(double t, TransitionProbabilityMatrix& P) const
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
        double rowsum = 0.0;
        for (size_t j=0; j<num_states; j++, ++p) 
        {
            double sum = 0.0;
            for (size_t s=0; s<num_states; s++)
            {
                sum += (*ptr++) * eigValExp[s];
            }
            
            sum = (sum < 0.0) ? 0.0 : sum;
            rowsum += sum;
            (*p) = sum;
        }

        // Normalize transition probabilities for row to sum to 1.0
        double* p2 = p - num_states;
        for (size_t j=0; j<num_states; j++, ++p2)
            *p2 /= rowsum;
    }
}


void RateMatrix_GTR::initFromString(const std::string &s)
{

    std::string tmp = s;
    StringUtilities::replaceSubstring(tmp, "[", "");
    StringUtilities::replaceSubstring(tmp, " ", "");
    
    std::vector<std::string> elements;
    StringUtilities::stringSplit(tmp, ",", elements);
    
    size_t n = size_t( sqrt(elements.size()) );
    
    delete the_rate_matrix;
    the_rate_matrix = new MatrixReal(n);
    for (size_t i=0; i<n; ++i)
    {
        for (size_t j=0; j<n; ++j)
        {
            (*the_rate_matrix)[i][j] = atof( elements[i*n+j].c_str() );
        }
    }
    
    stationary_freqs = calculateStationaryFrequencies();
    
    
    MatrixReal& m = *the_rate_matrix;
    // set the off-diagonal portions of the rate matrix
    for (size_t i=0, k=0; i<num_states; ++i)
    {
        for (size_t j=i+1; j<num_states; ++j)
        {
            double a = m[i][j] / stationary_freqs[j];
//            double b = m[j][i] / stationary_freqs[i];
//            
//            if ( a != b )
//            {
//                throw RbException("Unequal rates.");
//            }
            exchangeability_rates[k] = a;
            k++;
        }
    }
    
    needs_update = true;
    
}


/** Calculate the transition probabilities for the complex case */
void RateMatrix_GTR::tiProbsComplexEigens(double t, TransitionProbabilityMatrix& P) const
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
        double rowsum = 0.0;
        for (size_t j=0; j<num_states; j++) 
        {
            std::complex<double> sum = std::complex<double>(0.0, 0.0);
            for (size_t s=0; s<num_states; s++)
                sum += (*ptr++) * ceigValExp[s];

            double real_sum = (sum.real() < 0.0) ? 0.0 : sum.real();
            P[i][j] = real_sum;
            rowsum += real_sum;
        }
        // Normalize transition probabilities for row to sum to 1.0
        for (size_t j=0; j<num_states; j++)
            P[i][j] /= rowsum;
    }
}


/** Update the eigen system */
void RateMatrix_GTR::updateEigenSystem(void)
{
    
    theEigenSystem->update();
    calculateCijk();
    
}


void RateMatrix_GTR::update( void )
{
    
    if ( needs_update ) 
    {
        // compute the off-diagonal values
        computeOffDiagonal();
        
        // set the diagonal values
        setDiagonal();
        
        // rescale 
        rescaleToAverageRate( 1.0 );
        
        // now update the eigensystem
        updateEigenSystem();
        
        // clean flags
        needs_update = false;
    }
    
}


