#include "RateMatrix_PoMoKNrecurrentMutations.h"

#include <cassert>
#include <cstddef>

#include "RbException.h"
#include "RbMathCombinatorialFunctions.h"
#include "Cloneable.h"
#include "MatrixReal.h"
#include "TransitionProbabilityMatrix.h"
#include "RbVector.h"
#include "RbVectorImpl.h"

using namespace RevBayesCore;

/*
RateMatrix_PoMoKNrecurrentMutations::RateMatrix_PoMoKNrecurrentMutations(size_t num_states) : 
AbstractRateMatrix( num_states ), 
K( in_k ),
N( in_n ),
mu( in_nmr, 0.01 ),
phi( in_k, 1.0 )
{
    update();
}
*/

/** Construct rate matrix with n states, an exchangeability matrix, a simplex of equilibrium frequencies, and a virtual
 * population size */
RateMatrix_PoMoKNrecurrentMutations::RateMatrix_PoMoKNrecurrentMutations(long ns, long na, double nv, size_t n_mr, bool rm)  :
AbstractRateMatrix( ns ),
S( ns ),
K( na ),
V( nv ),
mu( n_mr, 0.01 ),
phi( na, 1.0 ),
R( rm )
{
    
    MatrixReal& m = *the_rate_matrix;

    // filling the matrix with zeros
    for (int i=0; i< ns; i++)
    {
        for (int j=0; j< ns; j++)
        {
            m[i][j] = 0.0;
        }
    }
    
    update();
}



/** Destructor */
RateMatrix_PoMoKNrecurrentMutations::~RateMatrix_PoMoKNrecurrentMutations(void)
{

}



double RateMatrix_PoMoKNrecurrentMutations::averageRate(void) const
{
    return 1.0;
}





void RateMatrix_PoMoKNrecurrentMutations::buildRateMatrix(void)
{

    /*  
     INFORMATION ABOUT THE PoMoKNrecurrentMutations RATE MATRIX

     It includes both fixed and polymorphic sites, but only biallelic states are considered.

     The pomo rate matrices defined here first list the fixed states {Va0}, {Va1} ...,
     these occupying positions 0:(K-1), and then polymorphic states.

     K alleles comprise (K*K-K)/2 pairwise combinations of alleles.
     This is the number of edges in the pomo state-space. Each edge comprises V-1 polymorphic states, 
     each of which represents a state of allelic frequencies in the population.
     Example for a random allele pair aiaj: {(V-1)ai,1aj}, {(V-2)ai,2aj}, ..., {1ai,(V-1)aj}

     The polymorphic edges are listed in the following order a0a1, a0a2, a0a3, ..., a(K-2)aK-1
     We say the a0a1 polymorphic states sit at edge 0. Thus, {(V-1)a0,1a1} sits at position K and 
     {1a0,(V-1)a1} sits V-2 positions further (i.e., K+V-2).

     More generally, the polymorphic states of the allele pair sitting at edge E occupy the positions
     [K+E*V-E]:[K+(E+1)*(V-1)-1].
     */

    MatrixReal& m = *the_rate_matrix;

    // populate diagonal elements with 0.0
    // otherwise diagonal elements get super negative (suspect 3 times neg than they should)
    // ask sebastian
    for (int i=0; i< S; i++)
    {
        m[i][i] = 0.0;
    }

    // Notes
    // K: number of alleles
    // V: virtual population size
    // mu : vector of mutation rates
    // phi: vector of fitness coefficients
    // R: boolean for recurrent mutations
    // std::cout << "(K,V,mu[0],phi[0],R)=" << " (" << K << "," << V << "," << mu[0] << "," << phi[0] << "," << R << ")\n";
    double rm = R*1.0;

    // first edge
    size_t E = 0;
     
    // these for loops go through the (K*K-K) edges of the pomo state-space
    // their represent all the possible pairwise combinations of alleles: a0a1, a1a2, ..., aK-2aK-1
    
    double mutation_i, mutation_j, fixation_i, fixation_j, freq_shift_i, freq_shift_j;
    for (int i=0; i<K; i++)
    {
        for (int j=(i+1); j<K; j++)
        {

            // mutations
            mutation_i = V * mu[2*E+1];
            mutation_j = V * mu[2*E];

            m[i][K+E*V-E]         = mutation_j ;   //{Vai} -> {(V-1)ai,1aj}
            m[j][K+(E+1)*(V-1)-1] = mutation_i ;    //{Vaj} -> {1ai,(V-1)aj}

            // diagonal elements
            m[i][i]  -= mutation_j;
            m[j][j]  -= mutation_i;

            // fixations
            fixation_i = (V-1.0)*phi[i]/(phi[j] + (V-1.0)*phi[i]) + rm*mu[2*E+1];
            fixation_j = (V-1.0)*phi[j]/(phi[i] + (V-1.0)*phi[j]) + rm*mu[2*E];

            m[K+E*V-E]        [i]  = fixation_i;  //{(V-1)ai,1aj} -> {Vai}
            m[K+(E+1)*(V-1)-1][j]  = fixation_j;  //{1ai,(V-1)aj} -> {Vaj}

            // diagonal elements
            m[K+E*V-E]        [K+E*V-E]         -= fixation_i;
            m[K+(E+1)*(V-1)-1][K+(E+1)*(V-1)-1] -= fixation_j;

            // the pomo rate matrix is entirely defined by fixations and mutations if V=2
            if (V>2)
            {

                // frequency shifts for all the other polymorphic states
                for (int v=1; v<(V-1); v++)
                {
                    // populates the rates from 1 to N-1 and N-1 to 1 at the same time:
                    freq_shift_i = v*(V-v)*phi[i]/(v*phi[i] + (V-v)*phi[j]) + rm*(V-v)*mu[2*E+1];
                    freq_shift_j = v*(V-v)*phi[j]/(v*phi[j] + (V-v)*phi[i]) + rm*(V-v)*mu[2*E];
                    
                    m[K+E*V-E-1+v]    [K+E*V-E+v]          = freq_shift_j;  //{(V-1)ai,1aj} -> {(V-2)ai,2aj}
                    m[K+(E+1)*(V-1)-v][K+(E+1)*(V-1)-v-1]  = freq_shift_i;  //{1ai,(V-1)aj} -> {2ai,(V-2)aj}

                    // diagonal elements
                    m[K+E*V-E-1+v]    [K+E*V-E-1+v]     -=  freq_shift_j;
                    m[K+(E+1)*(V-1)-v][K+(E+1)*(V-1)-v] -=  freq_shift_i;  

                }
            
            }

            //update edge
            E += 1; 

        }
    }

}


/** Calculate the transition probabilities */
void RateMatrix_PoMoKNrecurrentMutations::calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const
{

  // We use repeated squaring to quickly obtain exponentials, as in Poujol and Lartillot, Bioinformatics 2014.
  // Mayrose et al. 2010 also used this method for chromosome evolution (named the squaring and scaling method in Moler and Van Loan 2003).
  double t = rate * (startAge - endAge);
  exponentiateMatrixByScalingAndSquaring(t, P );

  return;
}


RateMatrix_PoMoKNrecurrentMutations* RateMatrix_PoMoKNrecurrentMutations::clone( void ) const
{
    return new RateMatrix_PoMoKNrecurrentMutations( *this );
}


std::vector<double> RateMatrix_PoMoKNrecurrentMutations::getStationaryFrequencies( void ) const
{
//  return stationaryVector;
    return calculateStationaryFrequencies();
}



void RateMatrix_PoMoKNrecurrentMutations::setNumberOfAlleles( long na )
{
    K = na;
    
    // set flags
    needs_update = true;
    
}


void RateMatrix_PoMoKNrecurrentMutations::setVirtualPopulationSize( long ni )
{
    V = ni;
    
    // set flags
    needs_update = true;
    
}



void RateMatrix_PoMoKNrecurrentMutations::setMu( const std::vector<double> &m )
{
    mu = m;
    
    // set flags
    needs_update = true;
}

void RateMatrix_PoMoKNrecurrentMutations::setPhi( const std::vector<double> &f )
{
    phi = f;
    
    // set flags
    needs_update = true;
}


void RateMatrix_PoMoKNrecurrentMutations::setRecurrentMutations( bool &r )
{
    R = r;
    
    // set flags
    needs_update = true;
}




void RateMatrix_PoMoKNrecurrentMutations::update( void )
{

    if ( needs_update )
    {

        buildRateMatrix();

        // rescale: not useful, same loglk.
        //rescaleToAverageRate( 1.0 );

        // clean flags
        needs_update = false;
    }
}
