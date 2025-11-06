//
//  RateMatrix_ReversiblePoMo.cpp
//
//  Created by Bastien Boussau on 4/4/14.
//  Copyright (c) 2014 Bastien Boussau. all rights reserved.
//

#include "RateMatrix_ReversiblePoMo.h"

#include <cassert>
#include <cstddef>

#include "RbException.h"
#include "Cloneable.h"
#include "MatrixReal.h"
#include "RbVector.h"
#include "RbVectorImpl.h"

using namespace RevBayesCore;

/** Construct rate matrix with n states */
RateMatrix_ReversiblePoMo::RateMatrix_ReversiblePoMo(size_t n) : AbstractRateMatrix( n ), rho( 6, 1.0 ), pi(4, 0.25), N( 10 ), matrixSize( n )
{
    // theEigenSystem       = new EigenSystem(the_rate_matrix);
    // c_ijk.resize(num_states * num_states * num_states);
    // cc_ijk.resize(num_states * num_states * num_states);
    update();
}

/** Construct rate matrix with n states and a matrix of mutation rates */
// RateMatrix_ReversiblePoMo::RateMatrix_ReversiblePoMo(size_t n,  const RateGenerator &Qmut) : AbstractRateMatrix( n ), Q_mut( static_cast< const Ratematrix_GtR& >( Qmut ) ), N( 10 ), matrixSize( n )
// {
//
//     stationary_freqs = Q_mut.getstationaryFrequencies();
//     update();
// }
//
// /** Construct rate matrix with n states, a matrix of mutation rates, and a virtual population size */
// RateMatrix_ReversiblePoMo::RateMatrix_ReversiblePoMo(const size_t n,  const RateGenerator &Qmut, const size_t vps )  : AbstractRateMatrix( n ), Q_mut( static_cast< const Ratematrix_GtR& > ( Qmut ) ), N( vps ), matrixSize( n )
// {
//   stationary_freqs = Q_mut.getstationaryFrequencies();
//   update();
// }


/** Construct rate matrix with n states, an exchangeability matrix, a simplex of equilibrium frequencies, and a virtual population size */
RateMatrix_ReversiblePoMo::RateMatrix_ReversiblePoMo(const size_t n,  const std::vector<double> &rh, const Simplex p, const size_t vps  )  : AbstractRateMatrix( n ), rho(rh), pi(p), N( vps ), matrixSize( n )
{

  // theEigenSystem       = new EigenSystem(the_rate_matrix);
  // c_ijk.resize(num_states * num_states * num_states);
  // cc_ijk.resize(num_states * num_states * num_states);

  update();

}

/** Copy constructor */
RateMatrix_ReversiblePoMo::RateMatrix_ReversiblePoMo(const RateMatrix_ReversiblePoMo& m) : AbstractRateMatrix( m ), rho(m.rho), pi(m.pi), N( m.N ), matrixSize( m.matrixSize )
{

    // theEigenSystem       = new EigenSystem( *m.theEigenSystem );
    // c_ijk                = m.c_ijk;
    // cc_ijk               = m.cc_ijk;
    //
    // theEigenSystem->setRateMatrixPtr(the_rate_matrix);
}


RateMatrix_ReversiblePoMo& RateMatrix_ReversiblePoMo::operator=(const RateMatrix_ReversiblePoMo &r)
{

  if (this != &r)
  {
    AbstractRateMatrix::operator=( r );
    rho = r.rho;
    pi = r.pi;
    N = r.N;
    matrixSize = r.matrixSize ;
    // delete theEigenSystem;
    //
    // theEigenSystem       = new EigenSystem( *r.theEigenSystem );
    // c_ijk                = r.c_ijk;
    // cc_ijk               = r.cc_ijk;
    //
    // theEigenSystem->setRateMatrixPtr(the_rate_matrix);
  }

  return *this;
}




/** Destructor */
RateMatrix_ReversiblePoMo::~RateMatrix_ReversiblePoMo(void)
{
    // delete theEigenSystem;

}


double RateMatrix_ReversiblePoMo::averageRate(void) const
{
    return 1.0;
}


// Code modified from IQ-tree:
void RateMatrix_ReversiblePoMo::decomposeState(int state, int &i, int &nt1, int &nt2)
{
    if (state < 4) {
        // Fixed a, C, G or t
        i = N;
        nt1 = state;
        nt2 = -1; // -1 for unknown nt
    } else if (state < 4+(N-1)) {
        // (ia,N-iC)
        i = state-3;
        nt1 = 0; // a
        nt2 = 1; // C
    } else if (state < 4+2*(N-1)) {
        // (ia,N-iG)
        i = state-3-(N-1);
        nt1 = 0; // a
        nt2 = 2; // G
    } else if (state < 4+3*(N-1)) {
        // (ia,N-it)
        i = state-3-2*(N-1);
        nt1 = 0; // a
        nt2 = 3; // t
    } else if (state < 4+4*(N-1)) {
        // (iC,N-iG)
        i = state-3-3*(N-1);
        nt1 = 1; // C
        nt2 = 2; // G
    } else if (state < 4+5*(N-1)) {
        // (iC,N-it)
        i = state-3-4*(N-1);
        nt1 = 1; // C
        nt2 = 3; // t
    } else if (state < 4+6*(N-1)) {
        // (iG,N-it)
        i = state-3-5*(N-1);
        nt1 = 2; // G
        nt2 = 3; // t
    } else {
      throw RbException("state exceeds limit");
    }
}

// // Code modified from IQ-tree:
// // Here, we return directly m_xy * \pi_xy, whereas the IQ-tree code returned only m_xy (I think)
// double RateMatrix_ReversiblePoMo::mutCoeff(int nt1, int nt2) {
//     assert(nt1!=nt2 && nt1<4 && nt2<4);
//     // if (nt2 < nt1) {
//     //     int tmp=nt1;
//     //     nt1=nt2;
//     //     nt2=tmp;
//     // }
//     //if (nt1==0)
//
//
//
//
//     return Q_mut.getRate(nt1, nt2);//mutation_prob[nt2-1];
//     //if (nt1==1) return mutation_prob[nt2+1];
//     //if (nt1==2) return mutation_prob[5];
//     //assert(0);
// }
//


double RateMatrix_ReversiblePoMo::mutCoeff(int nt1, int nt2) {

  // Inspired by code from Rui in revbayescode.txt
  // Order: aCGt
  // r[1][N+3]   <- abs(pi[2]*rho[1]) //aC?
  // r[1][2*N+2] <- abs(pi[3]*rho[2]) //aG
  // r[1][3*N+1] <- abs(pi[4]*rho[3]) //at
  //
  // r[2][5]     <- abs(pi[1]*rho[1]) //Ca
  // r[2][4*N]   <- abs(pi[3]*rho[4]) //CG
  // r[2][5*N-1] <- abs(pi[4]*rho[5]) //Ct
  //
  // r[3][N+4]   <- abs(pi[1]*rho[2]) //Ga
  // r[3][3*N+2] <- abs(pi[2]*rho[4]) //GC
  // r[3][6*N-2] <- abs(pi[4]*rho[6]) //Gt
  //
  // r[4][2*N+3] <- abs(pi[1]*rho[3]) //ta
  // r[4][4*N+1] <- abs(pi[2]*rho[5]) //tC
  // r[4][5*N]   <- abs(pi[3]*rho[6]) //tG
  double r = 0.0;
  if (nt1 == 0) {
    if (nt2 == 1) {
      r = pi[1]*rho[0];
    }
    else if (nt2 == 2) {
      r = pi[2]*rho[1];
    }
    else if (nt2 == 3) {
      r = pi[3]*rho[2];
    }
  }
  else if (nt1 == 1 ) {
    if (nt2 == 0) {
      r = pi[0]*rho[0];
    }
    else if (nt2 == 2) {
      r = pi[2]*rho[3];
    }
    else if (nt2 == 3) {
      r = pi[3]*rho[4];
    }

  }
  else if (nt1 == 2 ) {
    if (nt2 == 0) {
      r = pi[0]*rho[1];
    }
    else if (nt2 == 1) {
      r = pi[1]*rho[3];
    }
    else if (nt2 == 3) {
      r = pi[3]*rho[5];
    }

  }
  else if (nt1 == 3 ) {
    if (nt2 == 0) {
      r = pi[0]*rho[2];
    }
    else if (nt2 == 1) {
      r = pi[1]*rho[4];
    }
    else if (nt2 == 2) {
      r = pi[2]*rho[5];
    }

  }
//  assert(r != 0.0);
  return r;

}





// Code modified from IQ-tree:
double RateMatrix_ReversiblePoMo::computeProbBoundaryMutation(int state1, int state2)
{
    // the transition rate to the same state will be calculated by
    // (row_sum = 0).
    assert(state1 != state2);

    // Both states are decomposed into the abundance of the first
    // allele as well as the nucleotide of the first and the second
    // allele.
    int i1=0, i2=0, nt1=-1, nt2=-1, nt3=-1, nt4=-1;
    decomposeState(state1, i1, nt1, nt2);
    decomposeState(state2, i2, nt3, nt4);

    // Either the first nucleotides match or the first of state 1 with
    // the second of state 2 or the first of state 2 with the second
    // of state 1.  additionally, we have to consider fixed states as
    // special cases.
    if (nt1 == nt3 && (nt2==nt4 || nt2==-1 || nt4 == -1))
    {
        assert(i1 != i2); // because state1 != state2
        if (i1+1==i2)
            // e.g.: 2a8C -> 3a7C or 9a1C -> 10a
            // Changed Dom tue sep 29 13:31:02 CEst 2015
            // return double(i1*(N-i1)) / double(N*N);
            return double(i1*(N-i1)) / double(N);
        else if (i1-1 == i2)
            // e.g.: 3a7C -> 2a8C or 10a -> 9a1C
            if (nt2 == -1)
                // e.g. 10a -> 9a1C
                // return mutCoeff(nt1,nt4) * state_freq[nt4];
                return mutCoeff(nt1,nt4) ;//* freq_fixed_states[nt4];
            else
                // e.g. 9a1C -> 8a2C
                // Changed Dom tue sep 29 13:30:43 CEst 2015
                // return double(i1*(N-i1)) / double(N*N);
                return double(i1*(N-i1)) / double(N);
        else
            return 0.0;
    }
    else if (nt1 == nt4 && nt2 == -1 && i2 == 1)
    {
        // e.g.: 10G -> 1a9G
        //return mutCoeff(nt1,nt3) * state_freq[nt3];
        return mutCoeff(nt1,nt3) ;//* freq_fixed_states[nt3];
    }
    else if (nt2 == nt3  && i1 == 1 && nt4 == -1)
    {
        // E.g.: 1a9G -> 10G
        // Changed Dom tue sep 29 13:30:25 CEst 2015
        // return double(i1*(N-i1)) / double(N*N);
        return double(i1*(N-i1)) / double(N);
    }
    else
    {
        // 0 for all other transitions
        return 0.0;
    }
}





void RateMatrix_ReversiblePoMo::buildRateMatrix(void)
{
    // the destination state is in column, the starting state in rows.
//    int state1, state2;

    // activate this if frequencies of fixed states sum up to 1.0.
    // updateFreqFixedstate();

    //  computestateFreq();

    std::vector< double > stationaryFrequencies = getStationaryFrequencies();

    // Code modified from IQ-tree
    double tot_sum = 0.0;
    // Loop over rows (transition starting from state1).
    for (size_t state1 = 0; state1 < num_states; state1++)
    {
        double row_sum = 0.0;
        // Loop over columns in row state1 (transition to state2).
        for (size_t state2 = 0; state2 < num_states; state2++)
        {
            if (state2 != state1)
            {
                (*the_rate_matrix)[state1][state2] = computeProbBoundaryMutation(state1, state2);
                row_sum += (*the_rate_matrix)[state1][state2];
            }
        }
        tot_sum += stationaryFrequencies[state1]*row_sum;
        (*the_rate_matrix)[state1][state1] = -(row_sum);
    }
        
    // Normalize rate matrix such that one event happens per unit time.
    for (size_t state1 = 0; state1 < num_states; state1++)
    {
        for (size_t state2 = 0; state2 < num_states; state2++)
        {
            (*the_rate_matrix)[state1][state2] /= tot_sum;
        }
    }
    
}


/** Calculate the transition probabilities */
void RateMatrix_ReversiblePoMo::calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const
{


  // We use repeated squaring to quickly obtain exponentials, as in Poujol and Lartillot, Bioinformatics 2014.
      // Mayrose et al. 2010 also used this method for chromosome evolution (named the squaring and scaling method in Moler and Van Loan 2003).
  double t = rate * (startAge - endAge);
  exponentiateMatrixByScalingAndSquaring(t, P );

  //std::cout << "Transition probability matrix on branch of length "<<t<< " : " << P << std::endl;

    return;
}


RateMatrix_ReversiblePoMo* RateMatrix_ReversiblePoMo::clone( void ) const
{
    return new RateMatrix_ReversiblePoMo( *this );
}

std::vector<double> RateMatrix_ReversiblePoMo::getStationaryFrequencies( void ) const
{
  //indexes for the polymorphic states: a=0, C=1, G=2 and t=3
  // C++11:
  // std::vector<size_t> i1 { 0,0,0,1,1,2 };
  // std::vector<size_t> i2 { 1,2,3,2,3,3 };
  std::vector<size_t> i1 (6, 0);
  i1[3] = 1;
  i1[4] = 1;
  i1[5] = 2;

  std::vector<size_t> i2 (6, 3);
  i2[0] = 1;
  i2[1] = 2;
  i2[3] = 2;

  size_t numElements = 4+6*(N-1);
  std::vector<double> stationaryVector ( numElements ,0.0 );
  double sum = 0.0;

  //stationary frequencies of the monomorphic states
  //first 4 elements of the root distribution are proportional to pi
  for (size_t i = 0; i < 4; ++i) {
    stationaryVector[i] = pi[i];
    sum = sum + stationaryVector[i];

  }

  //stationary frequencies of the polymorphic states
  size_t ind = 4;
  for (size_t i = 0; i < 6; ++i) {
    for (size_t j = 1; j < N; ++j) {
      stationaryVector[ind] = rho[i]*pi[i1[i]]*pi[i2[i]]*N/(j*(N-j));
      sum += stationaryVector[ind];
      ind = ind+1;
    }
  }

  //normalization
  for (size_t i = 0; i < numElements; ++i )
  {
    stationaryVector[i] = stationaryVector[i]/sum;
  }

  return stationaryVector;
}

void RateMatrix_ReversiblePoMo::setExchangeabilityRates(const std::vector<double>& r)
{
    rho = r;
}



void RateMatrix_ReversiblePoMo::setStationaryFrequencies(const Simplex& s)
{
    pi = s;
}



void RateMatrix_ReversiblePoMo::update( void )
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
