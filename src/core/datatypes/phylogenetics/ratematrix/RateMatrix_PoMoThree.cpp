#include "CodonState.h"
#include "EigenSystem.h"
#include "MatrixComplex.h"
#include "MatrixReal.h"
#include "RateMatrix_PoMoThree.h"
#include "RbException.h"
#include "RbMathMatrix.h"
#include "TransitionProbabilityMatrix.h"

#include <cmath>
#include <string>
#include <iomanip>

using namespace RevBayesCore;

/** Construct rate matrix with n states */
RateMatrix_PoMoThree::RateMatrix_PoMoThree( void ) : TimeReversibleRateMatrix( 16 ),
    rho( 6, 1.0 ), 
	pi(4, 0.25), 
	sigma(4, 0.0),
	N( 10 )
{
    
    eigen_system       = new EigenSystem(the_rate_matrix);
    c_ijk.resize(num_states * num_states * num_states);
    cc_ijk.resize(num_states * num_states * num_states);
    
    update();
}


/** Copy constructor */
RateMatrix_PoMoThree::RateMatrix_PoMoThree(const RateMatrix_PoMoThree& m) : TimeReversibleRateMatrix( m ),
    rho( m.rho ),
    pi( m.pi ),
	sigma( m.sigma ),
    N( m.N )
{
    
    eigen_system        = new EigenSystem( *m.eigen_system );
    c_ijk               = m.c_ijk;
    cc_ijk              = m.cc_ijk;
    
    eigen_system->setRateMatrixPtr(the_rate_matrix);
}


/** Destructor */
RateMatrix_PoMoThree::~RateMatrix_PoMoThree(void)
{
    
    delete eigen_system;
}


RateMatrix_PoMoThree& RateMatrix_PoMoThree::operator=(const RateMatrix_PoMoThree &r)
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
		sigma               = r.sigma;
        N    				= r.N;
        
        eigen_system->setRateMatrixPtr(the_rate_matrix);
    }
    
    return *this;
}


RateMatrix_PoMoThree& RateMatrix_PoMoThree::assign(const Assignable &m)
{
    
    const RateMatrix_PoMoThree *rm = dynamic_cast<const RateMatrix_PoMoThree*>(&m);
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
void RateMatrix_PoMoThree::calculateCijk(void)
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
void RateMatrix_PoMoThree::calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const
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


RateMatrix_PoMoThree* RateMatrix_PoMoThree::clone( void ) const
{
    return new RateMatrix_PoMoThree( *this );
}


void RateMatrix_PoMoThree::computeOffDiagonal( void )
{
    
    MatrixReal& m = *the_rate_matrix;
    
	// parameter transformation
    // sigma2 and rho2 are the scaled parameters sigma and rho
    std::vector<double> sigma2 = transformSigma();
    std::vector<double> rho2   = transformRho();

    double rate = getExpectedNumberEvents(rho2,sigma2);

    for (int i=0; i<16 ; i++){
        for (int j=0; j<16; j++ ){
            m[i][j] = 0.0;
        }
    }
	
	
	// mutations
    m[0][5]  = rho2[0]*pi[1]/rate;    //mutation AC
    m[0][7]  = rho2[1]*pi[2]/rate;    //mutation AG
    m[0][9]  = rho2[2]*pi[3]/rate;    //mutation AT

    m[1][4]  = rho2[0]*pi[0]/rate;    //mutation CA
    m[1][11] = rho2[3]*pi[2]/rate;    //mutation CG
    m[1][13] = rho2[4]*pi[3]/rate;    //mutation CT

    m[2][6]  = rho2[1]*pi[0]/rate;    //mutation GA
    m[2][10] = rho2[3]*pi[1]/rate;    //mutation GC
    m[2][15] = rho2[5]*pi[3]/rate;    //mutation GT

    m[3][8]  = rho2[2]*pi[0]/rate;    //mutation TA
    m[3][12] = rho2[4]*pi[1]/rate;    //mutation TC
    m[3][14] = rho2[5]*pi[2]/rate;    //mutation TG



    //fixations: genetic drift and selection
    m[4][1]   = 2.0*(1+sigma2[1])/(3.0 * rate);             //C fixed
    m[4][5]   = 2.0*(1+sigma2[0])/(3.0 * rate);             //A increases in frequency

    m[5][0]   = 2.0*(1+sigma2[0])/(3.0 * rate);             //A fixed
    m[5][4]   = 2.0*(1+sigma2[1])/(3.0 * rate);             //C increases in frequency

    m[6][2]   = 2.0*(1+sigma2[2])/(3.0 * rate);             //G fixed
    m[6][7]   = 2.0*(1+sigma2[0])/(3.0 * rate);             //A increases in frequency

    m[7][0]   = 2.0*(1+sigma2[0])/(3.0 * rate);             //A fixed
    m[7][6]   = 2.0*(1+sigma2[2])/(3.0 * rate);             //G increases in frequency

    m[8][3]   = 2.0*(1+sigma2[3])/(3.0 * rate);             //T fixed
    m[8][9]   = 2.0*(1+sigma2[0])/(3.0 * rate);             //A increases in frequency
  
    m[9][0]   = 2.0*(1+sigma2[0])/(3.0 * rate);             //A fixed
    m[9][8]   = 2.0*(1+sigma2[3])/(3.0 * rate);             //T increases in frequency

    m[10][2]  = 2.0*(1+sigma2[2])/(3.0 * rate);             //G fixed
    m[10][11] = 2.0*(1+sigma2[1])/(3.0 * rate);             //C increases in frequency

    m[11][1]  = 2.0*(1+sigma2[1])/(3.0 * rate);             //C fixed
    m[11][10] = 2.0*(1+sigma2[2])/(3.0 * rate);             //G increases in frequency

    m[12][3]  = 2.0*(1+sigma2[3])/(3.0 * rate);             //T fixed
    m[12][13] = 2.0*(1+sigma2[1])/(3.0 * rate);             //C increases in frequency

    m[13][1]  = 2.0*(1+sigma2[1])/(3.0 * rate);             //C fixed
    m[13][12] = 2.0*(1+sigma2[3])/(3.0 * rate);             //T increases in frequency

    m[14][3]  = 2.0*(1+sigma2[3])/(3.0 * rate);             //T fixed
    m[14][15] = 2.0*(1+sigma2[2])/(3.0 * rate);             //G increases in frequency

    m[15][2]  = 2.0*(1+sigma2[2])/(3.0 * rate);             //G fixed
    m[15][14] = 2.0*(1+sigma2[3])/(3.0 * rate);             //T increases in frequency


    // set flags
    needs_update = true;
}


std::vector<double> RateMatrix_PoMoThree::transformSigma( void ) const
{
	// parameter scaling for sigma
    std::vector<double> sigma2(4, 0.0);

    sigma2[0] = exp(log(1+sigma[0])*(N-1.0)/2.0)-1.0;
    sigma2[1] = exp(log(1+sigma[1])*(N-1.0)/2.0)-1.0;
    sigma2[2] = exp(log(1+sigma[2])*(N-1.0)/2.0)-1.0;
    sigma2[3] = exp(log(1+sigma[3])*(N-1.0)/2.0)-1.0;

    return sigma2;
}


std::vector<double> RateMatrix_PoMoThree::transformRho( void ) const
{
	//parameter scaling for rho
    std::vector<double> rho2(6, 0.0);

    for (int n = 1; n < N; n++){
        rho2[0] = rho2[0] + exp(log(1+sigma[0])*(n-1.0) + log(1+sigma[1])*(N-n-1.0))*N/(n*(N-n));
        rho2[1] = rho2[1] + exp(log(1+sigma[0])*(n-1.0) + log(1+sigma[2])*(N-n-1.0))*N/(n*(N-n));
        rho2[2] = rho2[2] + exp(log(1+sigma[0])*(n-1.0) + log(1+sigma[3])*(N-n-1.0))*N/(n*(N-n));
        rho2[3] = rho2[3] + exp(log(1+sigma[1])*(n-1.0) + log(1+sigma[2])*(N-n-1.0))*N/(n*(N-n));
        rho2[4] = rho2[4] + exp(log(1+sigma[1])*(n-1.0) + log(1+sigma[3])*(N-n-1.0))*N/(n*(N-n));
        rho2[5] = rho2[5] + exp(log(1+sigma[2])*(n-1.0) + log(1+sigma[3])*(N-n-1.0))*N/(n*(N-n));
    }

    rho2[0] = rho[0]*rho2[0]/( 3.0*exp((N-1.0)*log(1+sigma[0])/2.0)/2.0 + 3.0*exp((N-1.0)*log(1+sigma[1])/2.0)/2.0 );
    rho2[1] = rho[1]*rho2[1]/( 3.0*exp((N-1.0)*log(1+sigma[0])/2.0)/2.0 + 3.0*exp((N-1.0)*log(1+sigma[2])/2.0)/2.0 );
    rho2[2] = rho[2]*rho2[2]/( 3.0*exp((N-1.0)*log(1+sigma[0])/2.0)/2.0 + 3.0*exp((N-1.0)*log(1+sigma[3])/2.0)/2.0 );
    rho2[3] = rho[3]*rho2[3]/( 3.0*exp((N-1.0)*log(1+sigma[1])/2.0)/2.0 + 3.0*exp((N-1.0)*log(1+sigma[2])/2.0)/2.0 );
    rho2[4] = rho[4]*rho2[4]/( 3.0*exp((N-1.0)*log(1+sigma[1])/2.0)/2.0 + 3.0*exp((N-1.0)*log(1+sigma[3])/2.0)/2.0 );
    rho2[5] = rho[5]*rho2[5]/( 3.0*exp((N-1.0)*log(1+sigma[2])/2.0)/2.0 + 3.0*exp((N-1.0)*log(1+sigma[3])/2.0)/2.0 );

    return rho2;
}


std::vector<double> RateMatrix_PoMoThree::getStationaryFrequencies( void ) const
{
  // TRIED RHO2 AND SIGMA2 AS INPUTS; REVBAYES RETURNS ERRORS
    // parameter transformation
    // sigma2 and rho2 are the parameters sigma and rho recalculated
    std::vector<double> sigma2 = transformSigma();
    std::vector<double> rho2   = transformRho();


    //normalization constant
    double nConstant;

    nConstant =pi[0]*pow(1+sigma2[0],2)+
               pi[1]*pow(1+sigma2[1],2)+
               pi[2]*pow(1+sigma2[2],2)+
               pi[3]*pow(1+sigma2[3],2)+
               pi[0]*pi[1]*rho2[0]*(1+sigma2[1])*3.0/2.0+
               pi[0]*pi[1]*rho2[0]*(1+sigma2[0])*3.0/2.0+
               pi[0]*pi[2]*rho2[1]*(1+sigma2[2])*3.0/2.0+
               pi[0]*pi[2]*rho2[1]*(1+sigma2[0])*3.0/2.0+
               pi[0]*pi[3]*rho2[2]*(1+sigma2[3])*3.0/2.0+
               pi[0]*pi[3]*rho2[2]*(1+sigma2[0])*3.0/2.0+
               pi[1]*pi[2]*rho2[3]*(1+sigma2[2])*3.0/2.0+
               pi[1]*pi[2]*rho2[3]*(1+sigma2[1])*3.0/2.0+
               pi[1]*pi[3]*rho2[4]*(1+sigma2[3])*3.0/2.0+
               pi[1]*pi[3]*rho2[4]*(1+sigma2[1])*3.0/2.0+
               pi[2]*pi[3]*rho2[5]*(1+sigma2[3])*3.0/2.0+
               pi[2]*pi[3]*rho2[5]*(1+sigma2[2])*3.0/2.0;

    //stationary vector
    std::vector<double> stationaryVector(16,0.0);

    stationaryVector[0]  = pi[0]*pow(1+sigma2[0],2)/nConstant;
    stationaryVector[1]  = pi[1]*pow(1+sigma2[1],2)/nConstant;
    stationaryVector[2]  = pi[2]*pow(1+sigma2[2],2)/nConstant;
    stationaryVector[3]  = pi[3]*pow(1+sigma2[3],2)/nConstant;

    stationaryVector[4]  = pi[0]*pi[1]*rho2[0]*(1+sigma2[1])*3.0/(2.0*nConstant);
    stationaryVector[5]  = pi[0]*pi[1]*rho2[0]*(1+sigma2[0])*3.0/(2.0*nConstant);
    stationaryVector[6]  = pi[0]*pi[2]*rho2[1]*(1+sigma2[2])*3.0/(2.0*nConstant);
    stationaryVector[7]  = pi[0]*pi[2]*rho2[1]*(1+sigma2[0])*3.0/(2.0*nConstant);
    stationaryVector[8]  = pi[0]*pi[3]*rho2[2]*(1+sigma2[3])*3.0/(2.0*nConstant);
    stationaryVector[9]  = pi[0]*pi[3]*rho2[2]*(1+sigma2[0])*3.0/(2.0*nConstant);
    stationaryVector[10] = pi[1]*pi[2]*rho2[3]*(1+sigma2[2])*3.0/(2.0*nConstant);
    stationaryVector[11] = pi[1]*pi[2]*rho2[3]*(1+sigma2[1])*3.0/(2.0*nConstant);
    stationaryVector[12] = pi[1]*pi[3]*rho2[4]*(1+sigma2[3])*3.0/(2.0*nConstant);
    stationaryVector[13] = pi[1]*pi[3]*rho2[4]*(1+sigma2[1])*3.0/(2.0*nConstant);
    stationaryVector[14] = pi[2]*pi[3]*rho2[5]*(1+sigma2[3])*3.0/(2.0*nConstant);
    stationaryVector[15] = pi[2]*pi[3]*rho2[5]*(1+sigma2[2])*3.0/(2.0*nConstant);

  return stationaryVector;

}


double RateMatrix_PoMoThree::getExpectedNumberEvents( std::vector<double> rho2, std::vector<double> sigma2 ) const
{

    //std::vector<double> sigma2 = transformSigma();
    //std::vector<double> rho2   = transformRho();

    double rate = 0.0;

    for (int n = 1; n < 4; n++){
        rate = rate + pi[0]*pi[1]*rho2[0]*exp(log(1+sigma2[0])*(n-1.0) + log(1+sigma2[1])*(3.0-n));
        rate = rate + pi[0]*pi[2]*rho2[1]*exp(log(1+sigma2[0])*(n-1.0) + log(1+sigma2[2])*(3.0-n));
        rate = rate + pi[0]*pi[3]*rho2[2]*exp(log(1+sigma2[0])*(n-1.0) + log(1+sigma2[3])*(3.0-n));
        rate = rate + pi[1]*pi[2]*rho2[3]*exp(log(1+sigma2[1])*(n-1.0) + log(1+sigma2[2])*(3.0-n));
        rate = rate + pi[1]*pi[3]*rho2[4]*exp(log(1+sigma2[1])*(n-1.0) + log(1+sigma2[3])*(3.0-n));
        rate = rate + pi[2]*pi[3]*rho2[5]*exp(log(1+sigma2[2])*(n-1.0) + log(1+sigma2[3])*(3.0-n));
    }

    //normalization constant
    double nConstant = 0.0;
    nConstant =pi[0]*pow(1+sigma2[0],2)+
               pi[1]*pow(1+sigma2[1],2)+
               pi[2]*pow(1+sigma2[2],2)+
               pi[3]*pow(1+sigma2[3],2)+
               pi[0]*pi[1]*rho2[0]*(1+sigma2[1])*3.0/2.0+
               pi[0]*pi[1]*rho2[0]*(1+sigma2[0])*3.0/2.0+
               pi[0]*pi[2]*rho2[1]*(1+sigma2[2])*3.0/2.0+
               pi[0]*pi[2]*rho2[1]*(1+sigma2[0])*3.0/2.0+
               pi[0]*pi[3]*rho2[2]*(1+sigma2[3])*3.0/2.0+
               pi[0]*pi[3]*rho2[2]*(1+sigma2[0])*3.0/2.0+
               pi[1]*pi[2]*rho2[3]*(1+sigma2[2])*3.0/2.0+
               pi[1]*pi[2]*rho2[3]*(1+sigma2[1])*3.0/2.0+
               pi[1]*pi[3]*rho2[4]*(1+sigma2[3])*3.0/2.0+
               pi[1]*pi[3]*rho2[4]*(1+sigma2[1])*3.0/2.0+
               pi[2]*pi[3]*rho2[5]*(1+sigma2[3])*3.0/2.0+
               pi[2]*pi[3]*rho2[5]*(1+sigma2[2])*3.0/2.0;

  rate = 2.0*rate/nConstant;
  //std::cout << "rate:" << rate << "\t";
  //std::cout << nConstant << "\t";
  return rate;
}


/** Calculate the transition probabilities for the real case */
void RateMatrix_PoMoThree::tiProbsEigens(double t, TransitionProbabilityMatrix& P) const
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
void RateMatrix_PoMoThree::tiProbsComplexEigens(double t, TransitionProbabilityMatrix& P) const
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


void RateMatrix_PoMoThree::setPopulationSize(double ps)
{
    N = ps;
    // set flags
    needs_update = true;
}


void RateMatrix_PoMoThree::setEquilibriumFrequencies( const std::vector<double> &f )
{
    pi = f;
    // set flags
    needs_update = true;
}

void RateMatrix_PoMoThree::setExchangeabilities( const std::vector<double> &r )
{
    rho = r;
    // set flags
    needs_update = true;
}

void RateMatrix_PoMoThree::setSelectionCoefficients( const std::vector<double> &s )
{
    sigma = s;
    // set flags
    needs_update = true;
}


/** Update the eigen system */
void RateMatrix_PoMoThree::updateEigenSystem(void)
{
    
    eigen_system->update();
    calculateCijk();
    
}


void RateMatrix_PoMoThree::update( void )
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



