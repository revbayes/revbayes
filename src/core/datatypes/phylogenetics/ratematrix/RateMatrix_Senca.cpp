#include <cstddef>
#include <cmath>
#include <complex>
#include <vector>

#include "CodonState.h"
#include "EigenSystem.h"
#include "MatrixComplex.h"
#include "MatrixReal.h"
#include "RateMatrix_Senca.h"
#include "RbException.h"
#include "TransitionProbabilityMatrix.h"
#include "AminoAcidState.h"
#include "Cloneable.h"
#include "DiscreteCharacterState.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "TimeReversibleRateMatrix.h"

using namespace RevBayesCore;

/** Construct rate matrix with n states */
RateMatrix_Senca::RateMatrix_Senca(size_t n) : TimeReversibleRateMatrix( n )
{
    
    theEigenSystem       = new EigenSystem(the_rate_matrix);
    c_ijk.resize(num_states * num_states * num_states);
    cc_ijk.resize(num_states * num_states * num_states);
    
    update();
}


/** Copy constructor */
RateMatrix_Senca::RateMatrix_Senca(const RateMatrix_Senca& m) : TimeReversibleRateMatrix( m )
{
    
    theEigenSystem       = new EigenSystem( *m.theEigenSystem );
    c_ijk                = m.c_ijk;
    cc_ijk               = m.cc_ijk;
    
    theEigenSystem->setRateMatrixPtr(the_rate_matrix);
}


/** Destructor */
RateMatrix_Senca::~RateMatrix_Senca(void)
{
    
    delete theEigenSystem;
}


RateMatrix_Senca& RateMatrix_Senca::operator=(const RateMatrix_Senca &r)
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


/** Do precalculations on eigenvectors */
void RateMatrix_Senca::calculateCijk(void)
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
void RateMatrix_Senca::calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const
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


RateMatrix_Senca* RateMatrix_Senca::clone( void ) const
{
    return new RateMatrix_Senca( *this );
}


double RateMatrix_Senca::computePreferenceRatio(const CodonState &from, const CodonState &to) const
{
    
    AminoAcidState aa_from = from.getAminoAcidState();
    AminoAcidState aa_to   = to.getAminoAcidState();
    
    // get the frequencies of the amino acids
//    double aa_freq_from = aa_freqs[aa_from.getState()];
//    double aa_freq_to   = aa_freqs[aa_to.getState()];
    
    // we need to get the frequency of the codon within that amino acid
//    double codon_freq_from = codon_freqs[aa_from.getState()][  ];
//    double codon_freq_to   = codon_freqs[aa_to.getState()][];
    
    // or
    double codon_freq_from = codon_freqs[from.getStateIndex()];
    double codon_freq_to   = codon_freqs[to.getStateIndex()];
    
    double x = codon_freq_from;
    double y = codon_freq_to;
    
    double mu = 0.0;
    if ( x == y )
    {
        mu = 1.0;
    }
    else if ( x == 0.0 )
    {
        mu = 100.0;
    }
    else if ( y == 0 )
    {
        mu = 0.0;
    }
    else
    {
        mu = -( log(x/y) / (1.0-(x/y)) );
    }
    
    return mu;
    
}


/**
 *
 *        0: codon changes in more than one codon position (or stop codons)
 *        1: synonymous transition
 *        2: synonymous transversion
 *        3: non-synonymous transition
 *        4: non-synonymous transversion
 */
void RateMatrix_Senca::computeOffDiagonal( void )
{
    
    MatrixReal& m = *the_rate_matrix;
    
    std::vector<double> rate = std::vector<double>(5,0);
    rate[0] = 0;
    rate[1] = kappa;
    rate[2] = 1;
    rate[3] = kappa*omega;
    rate[4] = omega;
    
    size_t rateClass = 0;
    
    // set the off-diagonal portions of the rate matrix
    for (size_t i=0; i<num_states; ++i)
    {
        CodonState c1 = CodonState(i);
        std::vector<unsigned int> codon_from = c1.getTripletStates();
        unsigned int codon_from_pos_1 = codon_from[0];
        unsigned int codon_from_pos_2 = codon_from[1];
        unsigned int codon_from_pos_3 = codon_from[2];
        
        AminoAcidState aa_from = c1.getAminoAcidState();
        
        for (size_t j=i+1; j<num_states; j++)
        {
            CodonState c2 = CodonState(j);
            
            std::vector<unsigned int> codon_to = c2.getTripletStates();
            unsigned int codon_to_pos_1 = codon_to[0];
            unsigned int codon_to_pos_2 = codon_to[1];
            unsigned int codon_to_pos_3 = codon_to[2];
            
            AminoAcidState aa_to = c2.getAminoAcidState();
            
            rateClass = 0;
            if (codon_from_pos_1 != codon_to_pos_1)
            {
                if ( (codon_from_pos_1 == 0 && codon_to_pos_1 == 2) || (codon_from_pos_1 == 2 && codon_to_pos_1 == 0) || // A <-> G
                    (codon_from_pos_1 == 1 && codon_to_pos_1 == 3) || (codon_from_pos_1 == 3 && codon_to_pos_1 == 1) )  // C <-> T
                {
                    rateClass = 1; // Transition at position 1
                }
                else
                {
                    rateClass = 2; // Transversion at position 1
                }
            }
            
            if (codon_from_pos_2 != codon_to_pos_2)
            {
                
                if (rateClass == -1)
                {
                    if ( (codon_from_pos_2 == 0 && codon_to_pos_2 == 2) || (codon_from_pos_2 == 2 && codon_to_pos_2 == 0) || // A <-> G
                        (codon_from_pos_2 == 1 && codon_to_pos_2 == 3) || (codon_from_pos_2 == 3 && codon_to_pos_2 == 1) )  // C <-> T
                    {
                        rateClass = 1; // Transition
                    }
                    else
                    {
                        rateClass = 2; // Transversion
                    }
                }
                else
                {
                    rateClass = 0; // Codon changes at more than one position
                }
                
            }
            
            if (codon_from_pos_3 != codon_to_pos_3)
            {
                
                if (rateClass == -1)
                {
                    if ( (codon_from_pos_3 == 0 && codon_to_pos_3 == 2) || (codon_from_pos_3 == 2 && codon_to_pos_3 == 0) || // A <-> G
                        (codon_from_pos_3 == 1 && codon_to_pos_3 == 3) || (codon_from_pos_3 == 3 && codon_to_pos_3 == 1) )  // C <-> T
                    {
                        rateClass = 1; // Transition
                    }
                    else
                    {
                        rateClass = 2; // Transversion
                    }
                    
                }
                else
                {
                    rateClass = 0; // Codon changes at more than one position
                }
                
            }
            
            if (rateClass != 0)
            {
                
                if (aa_from != aa_to)
                {
                    rateClass += 2; // Is a non-synonymous change
                }
                
            }
            
            double f = computePreferenceRatio( c1, c2 );
//            m[i][j] = rate[rateClass] * stationary_freqs[j] * f;
//            m[j][i] = rate[rateClass] * stationary_freqs[i] * f;
            m[i][j] = rate[rateClass] * f;
            m[j][i] = rate[rateClass] * f;
            
        }
    }
    
    // set flags
    needs_update = true;
}


void RateMatrix_Senca::computeStateFrequencies( void )
{
    
    throw RbException("Missing implementation in RateMatrix_senca::computeStateFrequencies.");
    
}


/** Calculate the transition probabilities for the real case */
void RateMatrix_Senca::tiProbsEigens(double t, TransitionProbabilityMatrix& P) const
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
            
            //            P[i][j] = (sum < 0.0) ? 0.0 : sum;
            (*p) = (sum < 0.0) ? 0.0 : sum;
        }
    }
}


/** Calculate the transition probabilities for the complex case */
void RateMatrix_Senca::tiProbsComplexEigens(double t, TransitionProbabilityMatrix& P) const {
    
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
void RateMatrix_Senca::updateEigenSystem(void)
{
    
    theEigenSystem->update();
    calculateCijk();
    
}


void RateMatrix_Senca::update( void )
{
    
    if ( needs_update )
    {
        // recompute the state frequencies
        computeStateFrequencies();
        
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



