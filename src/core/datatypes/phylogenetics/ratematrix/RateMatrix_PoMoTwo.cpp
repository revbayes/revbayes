#include "CodonState.h"
#include "EigenSystem.h"
#include "MatrixComplex.h"
#include "MatrixReal.h"
#include "RateMatrix_PoMoTwo.h"
#include "RbException.h"
#include "RbMathMatrix.h"
#include "TransitionProbabilityMatrix.h"

#include <cmath>
#include <string>
#include <iomanip>

using namespace RevBayesCore;

/** Construct rate matrix with n states */
RateMatrix_PoMoTwo::RateMatrix_PoMoTwo( void ) : TimeReversibleRateMatrix( 10 ),
    rho( 6, 1.0 ), 
	pi(4, 0.25), 
	N( 10 )
{
    
    eigen_system       = new EigenSystem(the_rate_matrix);
    c_ijk.resize(num_states * num_states * num_states);
    cc_ijk.resize(num_states * num_states * num_states);
    
    update();
}


/** Copy constructor */
RateMatrix_PoMoTwo::RateMatrix_PoMoTwo(const RateMatrix_PoMoTwo& m) : TimeReversibleRateMatrix( m ),
    rho( m.rho ),
    pi( m.pi ),
    N( m.N )
{
    
    eigen_system        = new EigenSystem( *m.eigen_system );
    c_ijk               = m.c_ijk;
    cc_ijk              = m.cc_ijk;
    
    eigen_system->setRateMatrixPtr(the_rate_matrix);
}


/** Destructor */
RateMatrix_PoMoTwo::~RateMatrix_PoMoTwo(void)
{
    
    delete eigen_system;
}


RateMatrix_PoMoTwo& RateMatrix_PoMoTwo::operator=(const RateMatrix_PoMoTwo &r)
{
    
    if (this != &r)
    {
        TimeReversibleRateMatrix::operator=( r );
        
        delete eigen_system;
        
        eigen_system        = new EigenSystem( *r.eigen_system );
        c_ijk               = r.c_ijk;
        cc_ijk              = r.cc_ijk;
        rho               	= r.rho;
        pi               	= r.pi;
        N    				= r.N;
        
        eigen_system->setRateMatrixPtr(the_rate_matrix);
    }
    
    return *this;
}


RateMatrix_PoMoTwo& RateMatrix_PoMoTwo::assign(const Assignable &m)
{
    
    const RateMatrix_PoMoTwo *rm = dynamic_cast<const RateMatrix_PoMoTwo*>(&m);
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
void RateMatrix_PoMoTwo::calculateCijk(void)
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
void RateMatrix_PoMoTwo::calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const
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


RateMatrix_PoMoTwo* RateMatrix_PoMoTwo::clone( void ) const
{
    return new RateMatrix_PoMoTwo( *this );
}


void RateMatrix_PoMoTwo::computeOffDiagonal( void )
{
    
    MatrixReal& m = *the_rate_matrix;
    
	// parameter transformation
    // rho2 are the scaled parameters rho
    std::vector<double> rho2 = transformRho();

    double rate = getExpectedNumberEvents(rho2);

    for (int i=0; i<10 ; i++){
        for (int j=0; j<10; j++ ){
            m[i][j] = 0.0;
        }
    }
	
	
	// mutations
    m[0][4] = rho2[0]*pi[1]/rate;    //mutation AC
    m[0][5] = rho2[1]*pi[2]/rate;    //mutation AG
    m[0][6] = rho2[2]*pi[3]/rate;    //mutation AT

    m[1][4] = rho2[0]*pi[0]/rate;    //mutation CA
    m[1][7] = rho2[3]*pi[2]/rate;    //mutation CG
    m[1][8] = rho2[4]*pi[3]/rate;    //mutation CT

    m[2][5] = rho2[1]*pi[0]/rate;    //mutation GA
    m[2][7] = rho2[3]*pi[1]/rate;    //mutation GC
    m[2][9] = rho2[5]*pi[3]/rate;    //mutation GT

    m[3][6] = rho2[2]*pi[0]/rate;    //mutation TA
    m[3][8] = rho2[4]*pi[1]/rate;    //mutation TC
    m[3][9] = rho2[5]*pi[2]/rate;    //mutation TG



    //fixations: genetic drift and selection
    m[4][0]   = 1.0/(2.0 * rate);             //A fixed
    m[4][1]   = 1.0/(2.0 * rate);             //C fixed

    m[5][0]   = 1.0/(2.0 * rate);             //A fixed
    m[5][2]   = 1.0/(2.0 * rate);             //G fixed

    m[6][0]   = 1.0/(2.0 * rate);             //A fixed
    m[6][3]   = 1.0/(2.0 * rate);             //T fixed

    m[7][1]   = 1.0/(2.0 * rate);             //C fixed
    m[7][2]   = 1.0/(2.0 * rate);             //G fixed

    m[8][1]   = 1.0/(2.0 * rate);             //C fixed
    m[8][3]   = 1.0/(2.0 * rate);             //T fixed
  
    m[9][2]   = 1.0/(2.0 * rate);             //G fixed
    m[9][3]   = 1.0/(2.0 * rate);             //T fixed

    // set flags
    needs_update = true;

}

std::vector<double> RateMatrix_PoMoTwo::transformRho( void ) const
{
	//parameter scaling for rho
    std::vector<double> rho2(6, 0.0);
    
    double harmonic_number = 0.0;

    for (int n = 1; n < N; n++){

        harmonic_number = harmonic_number + 1.0/n;

    }

    rho2[0] = rho[0]*harmonic_number;
    rho2[1] = rho[1]*harmonic_number;
    rho2[2] = rho[2]*harmonic_number;
    rho2[3] = rho[3]*harmonic_number;
    rho2[4] = rho[4]*harmonic_number;
    rho2[5] = rho[5]*harmonic_number;

    return rho2;
}


std::vector<double> RateMatrix_PoMoTwo::getStationaryFrequencies( void ) const
{
    // parameter transformation
    // rho2 are the parameters rho recalculated
    std::vector<double> rho2   = transformRho();


    //normalization constant
    double nConstant = 0.0;

    nConstant =pi[0]+
               pi[1]+
               pi[2]+
               pi[3]+
               pi[0]*pi[1]*rho2[0]*2.0+
               pi[0]*pi[2]*rho2[1]*2.0+
               pi[0]*pi[3]*rho2[2]*2.0+
               pi[1]*pi[2]*rho2[3]*2.0+
               pi[1]*pi[3]*rho2[4]*2.0+
               pi[2]*pi[3]*rho2[5]*2.0;

    //stationary vector
    std::vector<double> stationaryVector(10,0.0);

    stationaryVector[0]  = pi[0]/nConstant;
    stationaryVector[1]  = pi[1]/nConstant;
    stationaryVector[2]  = pi[2]/nConstant;
    stationaryVector[3]  = pi[3]/nConstant;
    stationaryVector[4]  = pi[0]*pi[1]*rho2[0]*2.0/nConstant;
    stationaryVector[5]  = pi[0]*pi[2]*rho2[1]*2.0/nConstant;
    stationaryVector[6]  = pi[0]*pi[3]*rho2[2]*2.0/nConstant;
    stationaryVector[7]  = pi[1]*pi[2]*rho2[3]*2.0/nConstant;
    stationaryVector[8]  = pi[1]*pi[3]*rho2[4]*2.0/nConstant;
    stationaryVector[9]  = pi[2]*pi[3]*rho2[5]*2.0/nConstant;

  return stationaryVector;

}


double RateMatrix_PoMoTwo::getExpectedNumberEvents( std::vector<double> rho2 ) const
{


    double rate = 0.0;

    rate =     pi[0]*pi[1]*rho2[0]*2.0+
               pi[0]*pi[2]*rho2[1]*2.0+
               pi[0]*pi[3]*rho2[2]*2.0+
               pi[1]*pi[2]*rho2[3]*2.0+
               pi[1]*pi[3]*rho2[4]*2.0+
               pi[2]*pi[3]*rho2[5]*2.0;

    //normalization constant
    double nConstant = 0.0;

    nConstant =pi[0]+
               pi[1]+
               pi[2]+
               pi[3]+
               pi[0]*pi[1]*rho2[0]*2.0+
               pi[0]*pi[2]*rho2[1]*2.0+
               pi[0]*pi[3]*rho2[2]*2.0+
               pi[1]*pi[2]*rho2[3]*2.0+
               pi[1]*pi[3]*rho2[4]*2.0+
               pi[2]*pi[3]*rho2[5]*2.0;

  rate = 2.0*rate/nConstant;
  //std::cout << "rate:" << rate << "\t";
  //std::cout << nConstant << "\t";
  return rate;
}


/** Calculate the transition probabilities for the real case */
void RateMatrix_PoMoTwo::tiProbsEigens(double t, TransitionProbabilityMatrix& P) const
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
            
            //            P[i][j] = (sum < 0.0) ? 0.0 : sum;
            (*p) = (sum < 0.0) ? 0.0 : sum;
        }
    }
}


/** Calculate the transition probabilities for the complex case */
void RateMatrix_PoMoTwo::tiProbsComplexEigens(double t, TransitionProbabilityMatrix& P) const
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


void RateMatrix_PoMoTwo::setPopulationSize(double ps)
{
    N = ps;
    // set flags
    needs_update = true;
}


void RateMatrix_PoMoTwo::setEquilibriumFrequencies( const std::vector<double> &f )
{
    pi = f;
    // set flags
    needs_update = true;
}

void RateMatrix_PoMoTwo::setExchangeabilities( const std::vector<double> &r )
{
    rho = r;
    // set flags
    needs_update = true;
}

/** Update the eigen system */
void RateMatrix_PoMoTwo::updateEigenSystem(void)
{
    
    eigen_system->update();
    calculateCijk();
    
}


void RateMatrix_PoMoTwo::update( void )
{
    
    if ( needs_update )
    {
        // compute the off-diagonal values
        computeOffDiagonal();
        
        // set the diagonal values
        setDiagonal();
        
        // rescale
        //double e_rate = getExpectedNumberEvents();
        //rescaleToAverageRate(e_rate);
        
        // now update the eigensystem
        updateEigenSystem();
        
        // clean flags
        needs_update = false;
    }
}



