#include <cstddef>
#include <cmath>
#include <complex>
#include <vector>

#include "EigenSystem.h"
#include "MatrixComplex.h"
#include "MatrixReal.h"
#include "RateMatrix_Empirical.h"
#include "RbException.h"
#include "TransitionProbabilityMatrix.h"
#include "Assignable.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "TimeReversibleRateMatrix.h"

using namespace RevBayesCore;

/** Construct rate matrix with n states */
RateMatrix_Empirical::RateMatrix_Empirical(size_t n) : TimeReversibleRateMatrix( n ){
    
    theEigenSystem       = new EigenSystem(the_rate_matrix);
    c_ijk.resize(num_states * num_states * num_states);
    cc_ijk.resize(num_states * num_states * num_states);
    
    update();
}


/** Copy constructor */
RateMatrix_Empirical::RateMatrix_Empirical(const RateMatrix_Empirical& m) : TimeReversibleRateMatrix( m )
{
    
    theEigenSystem       = new EigenSystem( *m.theEigenSystem );
    c_ijk                = m.c_ijk;
    cc_ijk               = m.cc_ijk;
    
    theEigenSystem->setRateMatrixPtr(the_rate_matrix);
}


/** Destructor */
RateMatrix_Empirical::~RateMatrix_Empirical(void)
{
    
    delete theEigenSystem;
}


RateMatrix_Empirical& RateMatrix_Empirical::operator=(const RateMatrix_Empirical &r)
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


/**
 * Assign the value of m to this instance. This function is our mechanism to call the assignment operator.
 *
 *
 */
RateMatrix_Empirical& RateMatrix_Empirical::assign(const Assignable &m)
{
    
    const RateMatrix_Empirical *rm = dynamic_cast<const RateMatrix_Empirical*>(&m);
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
void RateMatrix_Empirical::calculateCijk(void) {
    
    if ( theEigenSystem->isComplex() == false )
    {
        // real case
        const MatrixReal& ev  = theEigenSystem->getEigenvectors();
        const MatrixReal& iev = theEigenSystem->getInverseEigenvectors();
        double* pc = &c_ijk[0];
        for (size_t i=0; i<num_states; i++)
            for (size_t j=0; j<num_states; j++)
                for (size_t k=0; k<num_states; k++)
                    *(pc++) = ev[i][k] * iev[k][j];   
    }
    else
    {
        // complex case
        const MatrixComplex& cev  = theEigenSystem->getComplexEigenvectors();
        const MatrixComplex& ciev = theEigenSystem->getComplexInverseEigenvectors();
        std::complex<double>* pc = &cc_ijk[0];
        for (size_t i=0; i<num_states; i++)
            for (size_t j=0; j<num_states; j++)
                for (size_t k=0; k<num_states; k++)
                    *(pc++) = cev[i][k] * ciev[k][j];
    }
}


/** Calculate the transition probabilities */
void RateMatrix_Empirical::calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const {
    
    double t = rate * (startAge - endAge);
	if ( theEigenSystem->isComplex() == false )
		tiProbsEigens(t, P);
	else
		tiProbsComplexEigens(t, P);
}


RateMatrix_Empirical* RateMatrix_Empirical::clone( void ) const {
    return new RateMatrix_Empirical( *this );
}



/** Calculate the transition probabilities for the real case */
void RateMatrix_Empirical::tiProbsEigens(double t, TransitionProbabilityMatrix& P) const {
    
    // get a reference to the eigenvalues
    const std::vector<double>& eigenValue = theEigenSystem->getRealEigenvalues();
    
    // precalculate the product of the eigenvalue and the branch length
    std::vector<double> eigValExp(num_states);
	for (size_t s=0; s<num_states; s++)
		eigValExp[s] = exp(eigenValue[s] * t);
    
    // calculate the transition probabilities
	const double* ptr = &c_ijk[0];
	for (size_t i=0; i<num_states; i++) 
    {
		for (size_t j=0; j<num_states; j++) 
        {
			double sum = 0.0;
			for (size_t s=0; s<num_states; s++)
				sum += (*ptr++) * eigValExp[s];
			P[i][j] = (sum < 0.0) ? 0.0 : sum;
        }
    }
}


/** Calculate the transition probabilities for the complex case */
void RateMatrix_Empirical::tiProbsComplexEigens(double t, TransitionProbabilityMatrix& P) const {
    
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
				sum += (*ptr++) * ceigValExp[s];
			P[i][j] = (sum.real() < 0.0) ? 0.0 : sum.real();
        }
    }
}


/** Update the eigen system */
void RateMatrix_Empirical::updateEigenSystem(void) {
    
    theEigenSystem->update();
    calculateCijk();
    
}


void RateMatrix_Empirical::update( void ) {
    
}


