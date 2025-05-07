#include <cstddef>
#include <cmath>
#include <complex>
#include <iosfwd>
#include <vector>

#include "CodonState.h"
#include "EigenSystem.h"
#include "MatrixComplex.h"
#include "MatrixReal.h"
#include "RateMatrix_CodonSynonymousNonsynonymous.h"
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

/** Construct rate matrix with n states */
RateMatrix_CodonSynonymousNonsynonymous::RateMatrix_CodonSynonymousNonsynonymous( void ) : TimeReversibleRateMatrix( 61 ),
    codon_freqs(61,1.0/61)
{
    
    eigen_system       = new EigenSystem(the_rate_matrix);
    c_ijk.resize(num_states * num_states * num_states);
    cc_ijk.resize(num_states * num_states * num_states);
    
    update();
}


/** Copy constructor */
RateMatrix_CodonSynonymousNonsynonymous::RateMatrix_CodonSynonymousNonsynonymous(const RateMatrix_CodonSynonymousNonsynonymous& m) : TimeReversibleRateMatrix( m ),
    codon_freqs( m.codon_freqs )
{
    
    eigen_system        = new EigenSystem( *m.eigen_system );
    c_ijk               = m.c_ijk;
    cc_ijk              = m.cc_ijk;
    
    eigen_system->setRateMatrixPtr(the_rate_matrix);
}


/** Destructor */
RateMatrix_CodonSynonymousNonsynonymous::~RateMatrix_CodonSynonymousNonsynonymous(void)
{
    
    delete eigen_system;
}


RateMatrix_CodonSynonymousNonsynonymous& RateMatrix_CodonSynonymousNonsynonymous::operator=(const RateMatrix_CodonSynonymousNonsynonymous &r)
{
    
    if (this != &r)
    {
        TimeReversibleRateMatrix::operator=( r );
        
        delete eigen_system;
        
        eigen_system        = new EigenSystem( *r.eigen_system );
        c_ijk               = r.c_ijk;
        cc_ijk              = r.cc_ijk;
        codon_freqs         = r.codon_freqs;
        
        eigen_system->setRateMatrixPtr(the_rate_matrix);
    }
    
    return *this;
}


RateMatrix_CodonSynonymousNonsynonymous& RateMatrix_CodonSynonymousNonsynonymous::assign(const Assignable &m)
{
    
    const RateMatrix_CodonSynonymousNonsynonymous *rm = dynamic_cast<const RateMatrix_CodonSynonymousNonsynonymous*>(&m);
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
void RateMatrix_CodonSynonymousNonsynonymous::calculateCijk(void)
{
    
    if ( eigen_system->isComplex() == false )
    {
        // real case
        const MatrixReal& ev  = eigen_system->getEigenvectors();
        const MatrixReal& iev = eigen_system->getInverseEigenvectors();
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
        const MatrixComplex& cev  = eigen_system->getComplexEigenvectors();
        const MatrixComplex& ciev = eigen_system->getComplexInverseEigenvectors();
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
void RateMatrix_CodonSynonymousNonsynonymous::calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const
{
    double t = rate * (startAge - endAge);
    if ( eigen_system->isComplex() == false )
    {
        tiProbsEigens(t, P);
    }
    else
    {
        tiProbsComplexEigens(t, P);
    }
}


RateMatrix_CodonSynonymousNonsynonymous* RateMatrix_CodonSynonymousNonsynonymous::clone( void ) const
{
    return new RateMatrix_CodonSynonymousNonsynonymous( *this );
}

/**
 *
 *		0: codon changes in more than one codon position (or stop codons)
 *		1: synonymous transition
 *		2: synonymous transversion
 *		3: non-synonymous transition
 *		4: non-synonymous transversion
 */
void RateMatrix_CodonSynonymousNonsynonymous::computeOffDiagonal( void )
{
    
    MatrixReal& m = *the_rate_matrix;
    
    std::vector<double> rate = std::vector<double>(5,0);
    rate[0] = 0.0;
    rate[1] = 1.0;
    rate[2] = omega;

    int rate_class = 0;
    
    // set the off-diagonal portions of the rate matrix
    for (size_t i=0; i<num_states; ++i)
    {
        
        CodonState c1 = CodonState( CodonState::CODONS[i] );
        std::vector<unsigned int> codon_from = c1.getTripletStates();
        unsigned int codon_from_pos_1 = codon_from[0];
        unsigned int codon_from_pos_2 = codon_from[1];
        unsigned int codon_from_pos_3 = codon_from[2];
        
        AminoAcidState aa_from = c1.getAminoAcidState();
        
        for (size_t j=i+1; j<num_states; ++j)
        {

            CodonState c2 = CodonState( CodonState::CODONS[j] );
            
            std::vector<unsigned int> codon_to = c2.getTripletStates();
            unsigned int codon_to_pos_1 = codon_to[0];
            unsigned int codon_to_pos_2 = codon_to[1];
            unsigned int codon_to_pos_3 = codon_to[2];
            
            AminoAcidState aa_to = c2.getAminoAcidState();
            
            rate_class = -1;
            if (codon_from_pos_1 != codon_to_pos_1)
            {
                rate_class = 1;
            }
            
            if (codon_from_pos_2 != codon_to_pos_2)
            {
                
                if (rate_class == -1)
                {
                    rate_class = 1;
                }
                else
                {
                    rate_class = 0; // Codon changes at more than one position
                }
                
            }
            
            if (codon_from_pos_3 != codon_to_pos_3)
            {
                
                if (rate_class == -1)
                {
                    rate_class = 1;
                    
                }
                else
                {
                    rate_class = 0; // Codon changes at more than one position
                }
                
            }
            
            if (rate_class > 0 && aa_from != aa_to)
            {
                rate_class = 2; // Is a non-synonymous change
            }
            
            if ( rate_class > 0 )
            {
                m[i][j] = rate[rate_class] * codon_freqs[j];
                m[j][i] = rate[rate_class] * codon_freqs[i];
            }
            else
            {
                m[i][j] = 0.0;
                m[j][i] = 0.0;
            }
            
        }
    }
    
    // set flags
    needs_update = true;
}



/** Calculate the transition probabilities for the real case */
void RateMatrix_CodonSynonymousNonsynonymous::tiProbsEigens(double t, TransitionProbabilityMatrix& P) const
{
    
    // get a reference to the eigenvalues
    const std::vector<double>& eigenValue = eigen_system->getRealEigenvalues();
    
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
void RateMatrix_CodonSynonymousNonsynonymous::tiProbsComplexEigens(double t, TransitionProbabilityMatrix& P) const
{
    
    // get a reference to the eigenvalues
    const std::vector<double>& eigenValueReal = eigen_system->getRealEigenvalues();
    const std::vector<double>& eigenValueComp = eigen_system->getImagEigenvalues();
    
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


void RateMatrix_CodonSynonymousNonsynonymous::setOmega(double o)
{
    
    omega = o;
    
    // set flags
    needs_update = true;
    
}


void RateMatrix_CodonSynonymousNonsynonymous::setCodonFrequencies( const std::vector<double> &f )
{
    codon_freqs = f;
    
    // set flags
    needs_update = true;
}


/** Update the eigen system */
void RateMatrix_CodonSynonymousNonsynonymous::updateEigenSystem(void)
{
    
    eigen_system->update();
    calculateCijk();
    
}


void RateMatrix_CodonSynonymousNonsynonymous::update( void )
{
    
    if ( needs_update )
    {
        // compute the off-diagonal values
        computeOffDiagonal();
        
        // we also need to update the stationary frequencies
        this->stationary_freqs = codon_freqs;
        
        // set the diagonal values
        setDiagonal();
        
        // rescale
        rescaleToAverageRate( 3.0 );
        
        // now update the eigensystem
        updateEigenSystem();
        
        // clean flags
        needs_update = false;
    }
}


