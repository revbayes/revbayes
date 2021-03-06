#include "CodonState.h"
#include "EigenSystem.h"
#include "MatrixComplex.h"
#include "MatrixReal.h"
#include "RateMatrix_PoMo4N.h"
#include "RbException.h"
#include "RbMathMatrix.h"
#include "TransitionProbabilityMatrix.h"

#include <cmath>
#include <string>
#include <iomanip>

using namespace RevBayesCore;




/** Construct rate matrix with n states */
RateMatrix_PoMo4N::RateMatrix_PoMo4N( long num_states, long in_n ) : TimeReversibleRateMatrix( num_states ),
    N( in_n ),
    mu( 12 , 0.01 ),
    phi( 4 , 1.0 )
{
    
    eigen_system       = new EigenSystem(the_rate_matrix);
    c_ijk.resize(num_states * num_states * num_states);
    cc_ijk.resize(num_states * num_states * num_states);
    
    update();
}


/** Copy constructor */
RateMatrix_PoMo4N::RateMatrix_PoMo4N(const RateMatrix_PoMo4N& m) : TimeReversibleRateMatrix( m ),
    N( m.N ),
    mu( m.mu ),
    phi( m.phi )
{
    
    eigen_system        = new EigenSystem( *m.eigen_system );
    c_ijk               = m.c_ijk;
    cc_ijk              = m.cc_ijk;
    
    eigen_system->setRateMatrixPtr(the_rate_matrix);
}


/** Destructor */
RateMatrix_PoMo4N::~RateMatrix_PoMo4N(void)
{
    
    delete eigen_system;
}


RateMatrix_PoMo4N& RateMatrix_PoMo4N::operator=(const RateMatrix_PoMo4N &r)
{
    
    if (this != &r)
    {
        TimeReversibleRateMatrix::operator=( r );
        
        delete eigen_system;
        
        eigen_system        = new EigenSystem( *r.eigen_system );
        c_ijk               = r.c_ijk;
        cc_ijk              = r.cc_ijk;
        N                   = r.N;
        mu                  = r.mu;
        phi                 = r.phi;

        
        eigen_system->setRateMatrixPtr(the_rate_matrix);
    }
    
    return *this;
}


RateMatrix_PoMo4N& RateMatrix_PoMo4N::assign(const Assignable &m)
{
    
    const RateMatrix_PoMo4N *rm = dynamic_cast<const RateMatrix_PoMo4N*>(&m);
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
void RateMatrix_PoMo4N::calculateCijk(void)
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
void RateMatrix_PoMo4N::calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const
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


RateMatrix_PoMo4N* RateMatrix_PoMo4N::clone( void ) const
{
    return new RateMatrix_PoMo4N( *this );
}


/*populating the rate matrix*/
void RateMatrix_PoMo4N::computeOffDiagonal( void )
{
    
  MatrixReal& m = *the_rate_matrix;

  /*  
  INFORMATION ABOUT THE PoMo4N RATE MATRIX

  It includes both fixed and polymorphic sites, but only biallelic states are considered.

  The pomo rate matrices defined here first list the fixed states {Na0}, {Na1} ...,
  these occupying positions 0:3, and then polymorphic states.

  4 alleles comprise 6 pairwise combinations of alleles.
  This is the number of edges in the PoMo4 state-space. Each edge comprises N-1 polymorphic states, 
  each of which represents a state of allelic frequencies in the population (summing to N).
  Example for a random allele pair aiaj: {(N-1)ai,1aj}, {(N-2)ai,2aj}, ..., {1ai,(N-1)aj}

  The polymorphic edges are listed in the following order a0a1, a0a2, a0a3, ..., a(K-2)aK-1
  We say the a0a1 polymorphic states sit at edge 0. Thus, {(N-1)a0,1a1} sits at position 4 and 
  {1a0,(N-1)a1} sits N-2 positions further (i.e., 4+N-2).
  More generally, the polymorphic states of the allele pair sitting at edge E occupy the positions
  [4+E*N-E]:[4+(E+1)*(N-1)-1].
  */


  //populate rate matrix with 0.0
  // **total waste of time with sparse matrices like pomos**
  for (int i=0; i< num_states; i++){
    for (int j=0; j< num_states; j++){
      m[i][j] = 0.0;        
    }
  }


  //reciprocal of the population size
  double rN = 1.0/N;


  //mutations
  //AC
  m[0][4]      = mu[0];    //{NA} -> {(N-1)A,1C}
  m[1][N+2]    = mu[1];    //{NC} -> {1A,(N-1)C}

  //AG
  m[0][N+3]    = mu[2];    //{NA} -> {(N-1)A,1G}
  m[2][2*N+1]  = mu[3];    //{NG} -> {1A,(N-1)G}

  //AT
  m[0][2*N+2]  = mu[4];    //{NA} -> {(N-1)A,1T}
  m[3][3*N]    = mu[5];    //{NT} -> {1A,(N-1)T}

  //CG
  m[1][3*N+1]  = mu[6];    //{NC} -> {(N-1)C,1aG}
  m[2][4*N-1]  = mu[7];    //{NG} -> {1C,(N-1)aG}

  //CT
  m[1][4*N]    = mu[8];    //{NC} -> {(N-1)C,1T}
  m[3][5*N-2]  = mu[9];    //{NT} -> {1C,(N-1)T}

  //GT
  m[2][5*N-1]  = mu[10];    //{NG} -> {(N-1)G,1T}
  m[3][6*N-3]  = mu[11];    //{NT} -> {1G,(N-1)T}


  //fixations
  //AC
  m[4]    [0]   = (N-1.0)*phi[0]*rN;  //{(N-1)A,1C} -> {NA} 
  m[N+2]  [1]   = (N-1.0)*phi[1]*rN;  //{1A,(N-1)C} -> {NC} 

  //AG
  m[N+3]  [0]   = (N-1.0)*phi[0]*rN;  //{(N-1)A,1G} -> {NA} 
  m[2*N+1] [2]  = (N-1.0)*phi[2]*rN;  //{1A,(N-1)G} -> {NG} 

  //AT
  m[2*N+2][0]   = (N-1.0)*phi[0]*rN;  //{(N-1)A,1T} -> {NA} 
  m[3*N]   [3]  = (N-1.0)*phi[3]*rN;  //{1A,(N-1)T} -> {NT} 

  //CG
  m[3*N+1][1]   = (N-1.0)*phi[1]*rN;  //{(N-1)C,1G} -> {NC} 
  m[4*N-1] [2]  = (N-1.0)*phi[2]*rN;  //{1C,(N-1)G} -> {NG}

  //CT
  m[4*N]  [1]   = (N-1.0)*phi[1]*rN;  //{(N-1)C,1T} -> {NC} 
  m[5*N-2] [3]  = (N-1.0)*phi[3]*rN;  //{1C,(N-1)T} -> {NT} 

  //GT
  m[5*N-1][2]   = (N-1.0)*phi[2]*rN;  //{(N-1)G,1T} -> {NG} 
  m[6*N-3] [3]  = (N-1.0)*phi[3]*rN;  //{1G,(N-1)T} -> {NT} 

  //the pomo rate matrix is entirely defined by fixations and mutations if N=2
  if (N>2) {

    //frequency shifts from singletons
    //AC
    m[4]     [5]     = (N-1.0)*phi[1]*rN;  //{(N-1)A,1C} -> {(N-2)A,2C}
    m[N+2]   [N+1]   = (N-1.0)*phi[0]*rN;  //{1A,(N-1)C} -> {2A,(N-2)C}

    //AG
    m[N+3]   [N+4]   = (N-1.0)*phi[2]*rN;  //{(N-1)A,1G} -> {(N-2)A,2G}
    m[2*N+1] [2*N]   = (N-1.0)*phi[0]*rN;  //{1A,(N-1)G} -> {2A,(N-2)G}

    //AT
    m[2*N+2] [2*N+3] = (N-1.0)*phi[3]*rN;  //{(N-1)A,1T} -> {(N-2)A,2T}
    m[3*N]   [3*N-1] = (N-1.0)*phi[0]*rN;  //{1A,(N-1)T} -> {2A,(N-2)T}

    //CG
    m[3*N+1] [3*N+2] = (N-1.0)*phi[2]*rN;  //{(N-1)C,1G} -> {(N-2)C,2G}
    m[4*N-1] [4*N-2] = (N-1.0)*phi[1]*rN;  //{1C,(N-1)G} -> {2C,(N-2)G}

    //CT
    m[4*N]   [4*N+1] = (N-1.0)*phi[3]*rN;  //{(N-1)C,1T} -> {(N-2)C,2T}
    m[5*N-2] [5*N-3] = (N-1.0)*phi[1]*rN;  //{1C,(N-1)T} -> {2C,(N-2)T}

    //GT
    m[5*N-1] [5*N]   = (N-1.0)*phi[3]*rN;  //{(N-1)G,1T} -> {(N-2)G,2T}
    m[6*N-3] [6*N-4] = (N-1.0)*phi[2]*rN;  //{1G,(N-1)T} -> {2G,(N-2)T}


    //frequency shifts for all the other polymorphic states
    if (N>3) {

      //polymorphic states are populated in two fronts, thus the need for the middle frequency
      int S = N/2+1; 

      for (int n=2; n<S; n++){

        //populates the first half of the polymorphic edges
        //AC
        m[n+3]     [n+4]      = n*(N-n)*phi[1]*rN; //{nA,(N-n)C} -> {(n-1)A,(N-n+1)C}
        m[n+3]     [n+2]      = n*(N-n)*phi[0]*rN; //{nA,(N-n)C} -> {(n+1)A,(N-n-1)C}

        //AG
        m[N+n+2]   [N+n+3]    = n*(N-n)*phi[2]*rN; //{nA,(N-n)G} -> {(n-1)A,(N-n+1)G}
        m[N+n+2]   [N+n+1]    = n*(N-n)*phi[0]*rN; //{nA,(N-n)G} -> {(n+1)A,(N-n-1)G}

        //AT
        m[2*N+n+1] [2*N+n+2]  = n*(N-n)*phi[3]*rN; //{nA,(N-n)T} -> {(n-1)A,(N-n+1)T}
        m[2*N+n+1] [2*N+n]    = n*(N-n)*phi[0]*rN; //{nA,(N-n)T} -> {(n+1)A,(N-n-1)T}

        //CG
        m[3*N+n]   [3*N+n+1]  = n*(N-n)*phi[2]*rN; //{nC,(N-n)G} -> {(n-1)C,(N-n+1)G}
        m[3*N+n]   [3*N+n-1]  = n*(N-n)*phi[1]*rN; //{nC,(N-n)G} -> {(n+1)C,(N-n-1)G}

        //CT
        m[4*N+n-1] [4*N+n]    = n*(N-n)*phi[3]*rN; //{nC,(N-n)T} -> {(n-1)C,(N-n+1)T}
        m[4*N+n-1] [4*N+n-2]  = n*(N-n)*phi[1]*rN; //{nC,(N-n)T} -> {(n+1)C,(N-n-1)T}

        //GT
        m[5*N+n-2] [5*N+n-1]  = n*(N-n)*phi[3]*rN; //{nG,(N-n)T} -> {(n-1)G,(N-n+1)T}
        m[5*N+n-2] [5*N+n-3]  = n*(N-n)*phi[2]*rN; //{nG,(N-n)T} -> {(n+1)G,(N-n-1)T}


        //populates the second half of the polymorphic edges
        //AC
        m[N-n+3] [N-n+2]     = (N-n)*n*phi[0]*rN; //{(N-n)A,nC} -> {(N-n+1)A,(n-1)C}
        m[N-n+3] [N-n+4]     = (N-n)*n*phi[1]*rN; //{(N-n)A,nC} -> {(N-n-1)A,(n+1)C}

        //AG
        m[2*N-n+2] [2*N-n+1] = (N-n)*n*phi[0]*rN; //{(N-n)A,nG} -> {(N-n+1)A,(n-1)G}
        m[2*N-n+2] [2*N-n+3] = (N-n)*n*phi[2]*rN; //{(N-n)A,nG} -> {(N-n-1)A,(n+1)G}

        //AT
        m[3*N-n+1] [3*N-n]   = (N-n)*n*phi[0]*rN; //{(N-n)A,nT} -> {(N-n+1)A,(n-1)T}
        m[3*N-n+1] [3*N-n+2] = (N-n)*n*phi[3]*rN; //{(N-n)A,nT} -> {(N-n-1)A,(n+1)T}

        //CG
        m[4*N-n] [4*N-n-1]   = (N-n)*n*phi[1]*rN; //{(N-n)C,nG} -> {(N-n+1)C,(n-1)G}
        m[4*N-n] [4*N-n+1]   = (N-n)*n*phi[2]*rN; //{(N-n)C,nG} -> {(N-n-1)C,(n+1)G}

        //CT
        m[5*N-n-1] [5*N-n-2] = (N-n)*n*phi[1]*rN; //{(N-n)C,nT} -> {(N-n+1)C,(n-1)T}
        m[5*N-n-1] [5*N-n]   = (N-n)*n*phi[3]*rN; //{(N-n)C,nT} -> {(N-n-1)C,(n+1)T}

        //GT
        m[6*N-n-2] [6*N-n-3] = (N-n)*n*phi[2]*rN; //{(N-n)G,nT} -> {(N-n+1)G,(n-1)T}
        m[6*N-n-2] [6*N-n-1] = (N-n)*n*phi[3]*rN; //{(N-n)G,nT} -> {(N-n-1)G,(n+1)T}

      }

    }

  }


  // set flags
  needs_update = true;

}



/** Calculate the transition probabilities for the real case */
void RateMatrix_PoMo4N::tiProbsEigens(double t, TransitionProbabilityMatrix& P) const
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
void RateMatrix_PoMo4N::tiProbsComplexEigens(double t, TransitionProbabilityMatrix& P) const
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



void RateMatrix_PoMo4N::setN( long & ni )
{
    N = ni;
    
    // set flags
    needs_update = true;
    
}


void RateMatrix_PoMo4N::setMu( const std::vector<double> &m )
{
    mu = m;
    
    // set flags
    needs_update = true;
}

void RateMatrix_PoMo4N::setPhi( const std::vector<double> &f )
{
    phi = f;
    
    // set flags
    needs_update = true;
}



/** Update the eigen system */
void RateMatrix_PoMo4N::updateEigenSystem(void)
{
    
    eigen_system->update();
    calculateCijk();
    
}


void RateMatrix_PoMo4N::update( void )
{
    
    if ( needs_update )
    {
        // compute the off-diagonal values
        computeOffDiagonal();
        
        // set the diagonal values
        setDiagonal();
        
        // rescale
        //rescaleToAverageRate( 3.0 );
        
        // now update the eigensystem
        updateEigenSystem();
        
        // clean flags
        needs_update = false;
    }
}



