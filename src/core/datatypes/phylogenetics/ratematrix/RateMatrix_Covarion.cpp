//
//  RateMatrix_Covarion.cpp
//  revbayes-proj
//
//  Created by Michael Landis on 3/4/17.
//  Copyright © 2017 Michael Landis. All rights reserved.
//

#include <cstddef>
#include <cmath>
#include <complex>
#include <vector>

#include "EigenSystem.h"
#include "MatrixComplex.h"
#include "MatrixReal.h"
#include "RateMatrix_Covarion.h"
#include "RbException.h"
#include "TransitionProbabilityMatrix.h"
#include "GeneralRateMatrix.h"
#include "RbVector.h"
#include "RbVectorImpl.h"


using namespace RevBayesCore;

/** Construct rate matrix with n states */
RateMatrix_Covarion::RateMatrix_Covarion(size_t n, size_t k) : GeneralRateMatrix( n*k ),
rescale(true),
num_classes(k),
num_states_per_class(n),
rate_matrices(RbVector<MatrixReal>(k,MatrixReal(n,n,1))),
switch_rates(MatrixReal(k,k,1)),
clock_rates(RbVector<double>(k,1))
{
    
    theEigenSystem       = new EigenSystem(the_rate_matrix);
    c_ijk.resize(num_states * num_states * num_states);
    cc_ijk.resize(num_states * num_states * num_states);
    
    *the_rate_matrix = MatrixReal(this->num_states, this->num_states, 0.0);

    update();
}

RateMatrix_Covarion::RateMatrix_Covarion(size_t n, size_t k, bool r) : GeneralRateMatrix( n*k ),
rescale(r),
num_classes(k),
num_states_per_class(n),
rate_matrices(RbVector<MatrixReal>(k,MatrixReal(n,n,1))),
switch_rates(MatrixReal(k,k,1)),
clock_rates(RbVector<double>(k,1))
{
    
    theEigenSystem       = new EigenSystem(the_rate_matrix);
    c_ijk.resize(num_states * num_states * num_states);
    cc_ijk.resize(num_states * num_states * num_states);
    
    *the_rate_matrix = MatrixReal(this->num_states, this->num_states, 0.0);
    
    update();
}


/** Copy constructor */
RateMatrix_Covarion::RateMatrix_Covarion(const RateMatrix_Covarion& m) : GeneralRateMatrix( m )
{
    
    rescale              = m.rescale;
    theEigenSystem       = new EigenSystem( *m.theEigenSystem );
    c_ijk                = m.c_ijk;
    cc_ijk               = m.cc_ijk;
    num_states_per_class = m.num_states_per_class;
    num_classes          = m.num_classes;
    rate_matrices        = m.rate_matrices;
    switch_rates         = m.switch_rates;
    clock_rates          = m.clock_rates;
    
    theEigenSystem->setRateMatrixPtr(the_rate_matrix);
}


/** Destructor */
RateMatrix_Covarion::~RateMatrix_Covarion(void)
{
    
    delete theEigenSystem;
}


RateMatrix_Covarion& RateMatrix_Covarion::operator=(const RateMatrix_Covarion &r)
{
    
    if (this != &r)
    {
        GeneralRateMatrix::operator=( r );
        
        delete theEigenSystem;
        
        theEigenSystem       = new EigenSystem( *r.theEigenSystem );
        c_ijk                = r.c_ijk;
        cc_ijk               = r.cc_ijk;
        num_states_per_class = r.num_states_per_class;
        num_classes          = r.num_classes;
        rate_matrices        = r.rate_matrices;
        switch_rates         = r.switch_rates;
        clock_rates          = r.clock_rates;

        
        theEigenSystem->setRateMatrixPtr(the_rate_matrix);
        
    }
    
    return *this;
}


std::vector<int> RateMatrix_Covarion::get_emitted_letters() const
{
    std::vector<int> emit(num_states);
    for(int i=0;i<num_states;i++)
        emit[i] = i % num_states_per_class;
    return emit;
}

/** Do precalculations on eigenvectors */
void RateMatrix_Covarion::calculateCijk(void)
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
void RateMatrix_Covarion::calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const
{
    // The eigensystem code was returning NaN likelihood values when transition rates
    // were close to 0.0, so now we use the scaling and squaring method.
    double t = rate * (startAge - endAge);
    exponentiateMatrixByScalingAndSquaring(t, P);
    
    //    double t = rate * (startAge - endAge);
    //	if ( theEigenSystem->isComplex() == false )
    //    {
    //		tiProbsEigens(t, P);
    //    }
    //	else
    //    {
    //		tiProbsComplexEigens(t, P);
    //    }
    
}

RateMatrix_Covarion* RateMatrix_Covarion::clone( void ) const
{
    return new RateMatrix_Covarion( *this );
}


void RateMatrix_Covarion::fillRateMatrix( void )
{

    *the_rate_matrix = MatrixReal(this->num_states, this->num_states, 0.0);
    
    MatrixReal& m = *the_rate_matrix;

    
    std::vector<double> row_sums(this->num_states, 0.0);
    
    // fill in the diagonal blocks
    for (size_t i = 0; i < num_classes; i++)
    {
        // fill in diagonal block matrices except the diagonal term, which is the negative sum of rates
        size_t block_i = i * num_states_per_class;
        for (size_t j = 0; j < num_states_per_class; j++)
        {
            for (size_t k = 0; k < num_states_per_class; k++)
            {
                if (j == k)
                    continue;
                
                double v = rate_matrices[i][j][k] * clock_rates[i];
                m[block_i + j][block_i + k] = v;
                row_sums[block_i + j] += v;
            }
        }
        
        // switch rates on the diagonals of the off-diagonal blocks
        for (size_t j = 0; j < num_classes; j++)
        {
            if (i == j)
                continue;
            
            size_t block_j = j * num_states_per_class;
            for (size_t k = 0; k < num_states_per_class; k++)
            {
                double v = switch_rates[i][j];
                m[block_i + k][block_j + k] = v;
                row_sums[block_i + k] += v;
            }
        }
    }
    
    // set diagonal negative sums of rates
    for (size_t i = 0; i < row_sums.size(); i++)
    {
        m[i][i] = -row_sums[i];
    }
    
    // set flags
    needs_update = true;
}

void RateMatrix_Covarion::setRateMatrices(const RbVector<MatrixReal> &rm)
{
    rate_matrices = rm;
    num_states_per_class = rm[0].getDim();
    // set flags
    needs_update = true;
}

void RateMatrix_Covarion::setSwitchRates(const MatrixReal &sr)
{
    switch_rates = sr;
    num_classes = sr.getDim();
    // set flags
    needs_update = true;
}

void RateMatrix_Covarion::setClockRates(const RbVector<double>& cr)
{
    clock_rates = cr;
    num_classes = cr.size();
    // set flags
    needs_update = true;
}


/** Calculate the transition probabilities for the real case */
void RateMatrix_Covarion::tiProbsEigens(double t, TransitionProbabilityMatrix& P) const
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
            
            //			P[i][j] = (sum < 0.0) ? 0.0 : sum;
            (*p) = (sum < 0.0) ? 0.0 : sum;
        }
        
    }
    
}


/** Calculate the transition probabilities for the complex case */
void RateMatrix_Covarion::tiProbsComplexEigens(double t, TransitionProbabilityMatrix& P) const
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
    
}


/** Update the eigen system */
void RateMatrix_Covarion::updateEigenSystem(void)
{
    
    theEigenSystem->update();
    calculateCijk();
    
}


void RateMatrix_Covarion::update( void )
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
        
        // now update the eigensystem
        //        updateEigenSystem();
        
        // clean flags
        needs_update = false;
    }
    
}

