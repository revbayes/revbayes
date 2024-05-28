#include "RateMatrix_PoMo.h"

#include "MatrixReal.h"
#include "RbException.h"
#include "RbMathCombinatorialFunctions.h"
#include "TransitionProbabilityMatrix.h"
#include "Cloneable.h"
#include "RbVector.h"
#include "RbVectorImpl.h"

using namespace RevBayesCore;

/** Construct rate matrix with n states, virtual population size, mutation rates, selection coefficients */
RateMatrix_PoMo::RateMatrix_PoMo(size_t n, size_t vps, const std::vector<double> &mr, const std::vector<double> &sc) : AbstractRateMatrix( n + size_t(RbMath::kchoose2(n))*(vps-1) ),
    N( vps ),
    matrix_size( n + size_t(RbMath::kchoose2(n))*(vps-1) ),
    num_raw_states( n )
{
    std::vector<double> temp (n, 0.0);
    for (size_t i = 0; i<n ; ++i)
    {
        mu.push_back(temp);
        s.push_back(1.0);
    }
    
    if ( mr.size() > 0 )
    {
        setMutationRates(mr);
    }
    if ( sc.size() > 0 )
    {
        setSelectionCoefficients(sc);
    }
    
    update();
}

/** Construct rate matrix with n states, a matrix of mutation rates, and a vector of selection coefficients */
RateMatrix_PoMo::RateMatrix_PoMo(size_t n, size_t vps, const RateGenerator &mm, const std::vector<double> sc) : AbstractRateMatrix( n + size_t(RbMath::kchoose2(n))*(vps-1) ),
    N( vps ),
    matrix_size( n + size_t(RbMath::kchoose2(n))*(vps-1) ),
    num_raw_states( n )
{
    std::vector<double> temp (n, 0.0);
    for (size_t i = 0; i<n ; ++i)
    {
        mu.push_back(temp);
        s.push_back(1.0);
    }
    setMutationRates(mm);
    setSelectionCoefficients(sc);
    
    update();
}


/** Destructor */
RateMatrix_PoMo::~RateMatrix_PoMo(void)
{
    
}


double RateMatrix_PoMo::averageRate(void) const
{
    return 1.0;
}

void RateMatrix_PoMo::buildRateMatrix(void) 
{
    
    // compute auxilliary variables
    double N2 = 1.0;//(double) (N*N);
    int Nminus1 = (int)N-1;
    double Nminus1d = (double) Nminus1;
    for (size_t i = 0 ; i < num_raw_states; i++)
    {
        mu[i][i] = 0.0;
    }
    
    // calculate the transition probabilities
    for (size_t i=0; i< matrix_size; i++)
    {
        //The first 4 states are the monomorphic states; we can't directly change from one into another one
        for (size_t j=0; j< matrix_size; j++)
        {
            (*the_rate_matrix)[i][j] = 0.0;
        }
        
    }
    
    // Change from a monomorphic into a polymorphic state
    //(i.e. the first 4 lines in the matrix )
    //The 4..4+N-1 states are the AC matrix
    //Only 2 entries can differ from 0, (N-1)A and (N-1)C
    //(N-1)A can only come from monomorphic state A, i.e. i=0
    //(N-1)A is at the end of the submatrix, j=4+N-1
    (*the_rate_matrix)[0][4 + Nminus1 - 1] = N2 * mu[0][1];
    //(N-1)C can only come from monomorphic state C, i.e. i=1
    //(N-1)C is at the begining of the submatrix, j=4+N
    (*the_rate_matrix)[1][4] = N2 * mu[1][0];
    
    //The 4+Nminus1..4+2Nminus1 states are the AG matrix
    //Only 2 entries can differ from 0, (N-1)A and (N-1)G
    //(N-1)A can only come from monomorphic state A, i.e. i=0
    //(N-1)A is at the end of the submatrix, j=4+2*Nminus1
    (*the_rate_matrix)[0][4 + 2*Nminus1 -1] = N2 * mu[0][2];
    //(N-1)G can only come from monomorphic state G, i.e. i=1
    //(N-1)G is at the begining of the submatrix, j=4+N
    (*the_rate_matrix)[2][4 + Nminus1] = N2 * mu[2][0];

    //The 4+2Nminus1..4+3Nminus1 states are the AT matrix
    //Only 2 entries can differ from 0, (N-1)A and (N-1)T
    //(N-1)A can only come from monomorphic state A, i.e. i=0
    //(N-1)A is at the end of the submatrix, j=4+3*Nminus1
    (*the_rate_matrix)[0][4 + 3*Nminus1 -1] = N2 * mu[0][3];
    //(N-1)T can only come from monomorphic state T, i.e. i=1
    //(N-1)T is at the begining of the submatrix, j=4+N
    (*the_rate_matrix)[3][4 + 2*Nminus1] = N2 * mu[3][0];

    //The 4+3Nminus1..4+4Nminus1 states are the CG matrix
    //Only 2 entries can differ from 0, (N-1)C and (N-1)G
    //(N-1)C can only come from monomorphic state C, i.e. i=0
    //(N-1)C is at the end of the submatrix, j=4+4*Nminus1
    (*the_rate_matrix)[1][4 + 4*Nminus1 -1] = N2 * mu[1][2];
    //(N-1)G can only come from monomorphic state G, i.e. i=1
    //(N-1)G is at the begining of the submatrix, j=4+N
    (*the_rate_matrix)[2][4 + 3*Nminus1] = N2 * mu[2][1];

    //The 4+4Nminus1..4+5Nminus1 states are the CT matrix
    //Only 2 entries can differ from 0, (N-1)C and (N-1)T
    //(N-1)C can only come from monomorphic state C, i.e. i=0
    //(N-1)C is at the end of the submatrix, j=4+5*Nminus1
    (*the_rate_matrix)[1][4 + 5*Nminus1 -1] = N2 * mu[1][3];
    //(N-1)T can only come from monomorphic state T, i.e. i=1
    //(N-1)T is at the begining of the submatrix, j=4+N
    (*the_rate_matrix)[3][4 + 4*Nminus1] = N2 * mu[3][1];

    //The 4+5Nminus1..4+6Nminus1 states are the GT matrix
    //Only 2 entries can differ from 0, (N-1)G and (N-1)T
    //(N-1)G can only come from monomorphic state G, i.e. i=0
    //(N-1)G is at the end of the submatrix, j=4+5*Nminus1
    (*the_rate_matrix)[2][4 + 6*Nminus1 -1] = N2 * mu[2][3];
    //(N-1)T can only come from monomorphic state T, i.e. i=1
    //(N-1)T is at the begining of the submatrix, j=4+N
    (*the_rate_matrix)[3][4 + 5*Nminus1] = N2 * mu[3][2];


    //Now we move from a polymorphic state to a monomorphic state
    //(i.e. the first four columns in the matrix)
    //The [4..4+Nminus1[ states are the AC matrix
    //Only 2 entries can differ from 0, (N-1)A going to mono A and (N-1)C going to mono C
    //(N-1)A can only go to monomorphic state A, i.e. j=0
    //(N-1)A is at the end of the submatrix, i=4+N-1
    
    (*the_rate_matrix)[4 + Nminus1 - 1][0] = computeEntryFromMoranProcessWithSelection(0, 1, Nminus1d);
//    double temp = (N-1)*(1+s[0]-s[1]);
  //  (*the_rate_matrix)[4 + Nminus1 - 1][0] = temp / ( temp + 1) * (1) / N;
    //(N-1)C can only go to monomorphic state C, i.e. j=1
    //(N-1)C is at the begining of the submatrix, i=4
   // temp = (N-1)*(1+s[1]-s[0]);
    (*the_rate_matrix)[4 ][1] = computeEntryFromMoranProcessWithSelection(1, 0, Nminus1d);
  //  (*the_rate_matrix)[4 ][1] = temp / ( temp + 1) * (1) / N;

    //The 4+Nminus1..4+2Nminus1 states are the AG matrix
    //Only 2 entries can differ from 0, (N-1)A going to mono A and (N-1)G going to mono G
    //(N-1)A can only go to monomorphic state A, i.e. j=0
    //(N-1)A is at the end of the submatrix, i=4+2*(Nminus1)
  //  temp = (N-1)*(1+s[0]-s[2]);
    (*the_rate_matrix)[4 + 2*Nminus1 - 1][0] = computeEntryFromMoranProcessWithSelection(0, 2, Nminus1d);
  //  (*the_rate_matrix)[4 + 2*Nminus1 - 1][0] = temp / ( temp + 1) * (1) / N;
    //(N-1)G can only go to monomorphic state G, i.e. j=2
    //(N-1)G is at the begining of the submatrix, i=4
  //  temp = (N-1)*(1+s[2]-s[0]);
    (*the_rate_matrix)[4 + Nminus1 ][2] = computeEntryFromMoranProcessWithSelection(2, 0, Nminus1d);
   // (*the_rate_matrix)[4 + Nminus1 ][2] = temp / ( temp + 1) * (1) / N;

    //The 4+2Nminus1..4+3Nminus1 states are the AT matrix
    //Only 2 entries can differ from 0, (N-1)A going to mono A and (N-1)T going to mono T
    //(N-1)A can only go to monomorphic state A, i.e. j=0
    //(N-1)A is at the end of the submatrix, i=4+3*(Nminus1)
//    temp = (N-1)*(1+s[0]-s[3]);
    (*the_rate_matrix)[4 + 3*Nminus1 - 1][0] = computeEntryFromMoranProcessWithSelection(0, 3, Nminus1d);
   // (*the_rate_matrix)[4 + 3*Nminus1 - 1][0] = temp / ( temp + 1) * (1) / N;
    //(N-1)T can only go to monomorphic state T, i.e. j=3
    //(N-1)T is at the begining of the submatrix, i=4+ 2*Nminus1
   // temp = (N-1)*(1+s[3]-s[0]);
    (*the_rate_matrix)[4 + 2*Nminus1 ][3] = computeEntryFromMoranProcessWithSelection(3, 0, Nminus1d);
   // (*the_rate_matrix)[4 + 2*Nminus1 ][3] = temp / ( temp + 1) * (1) / N;

    //The 4+3Nminus1..4+4Nminus1 states are the CG matrix
    //Only 2 entries can differ from 0, (N-1)C going to mono C and (N-1)G going to mono G
    //(N-1)C can only go to monomorphic state C, i.e. j=1
    //(N-1)C is at the end of the submatrix, i=4+4*(Nminus1)
  //  temp = (N-1)*(1+s[1]-s[2]);
    (*the_rate_matrix)[4 + 4*Nminus1 - 1][1] = computeEntryFromMoranProcessWithSelection(1, 2, Nminus1d);
//    (*the_rate_matrix)[4 + 4*Nminus1 - 1][1] = temp / ( temp + 1) * (1) / N;
    //(N-1)G can only go to monomorphic state G, i.e. j=2
    //(N-1)G is at the begining of the submatrix, i=4
 //   temp = (N-1)*(1+s[2]-s[1]);
     (*the_rate_matrix)[4 + 3*Nminus1 ][2] = computeEntryFromMoranProcessWithSelection(2, 1, Nminus1d);
  //  (*the_rate_matrix)[4 + 3*Nminus1 ][2] = temp / ( temp + 1) * (1) / N;
    
    //The 4+4Nminus1..4+5Nminus1 states are the CT matrix
    //Only 2 entries can differ from 0, (N-1)C going to mono C and (N-1)T going to mono T
    //(N-1)C can only go to monomorphic state C, i.e. j=1
    //(N-1)C is at the end of the submatrix, i=4+5*(Nminus1)
 //   temp = (N-1)*(1+s[1]-s[3]);
   (*the_rate_matrix)[4 + 5*Nminus1 - 1][1] =computeEntryFromMoranProcessWithSelection(1, 3, Nminus1d);
  //  (*the_rate_matrix)[4 + 5*Nminus1 - 1][1] = temp / ( temp + 1) * (1) / N;
    //(N-1)T can only go to monomorphic state T, i.e. j=3
    //(N-1)T is at the begining of the submatrix, i=4
 //   temp = (N-1)*(1+s[3]-s[1]);
    (*the_rate_matrix)[4 + 4*Nminus1 ][3] =computeEntryFromMoranProcessWithSelection(3, 1, Nminus1d);
  //  (*the_rate_matrix)[4 + 4*Nminus1 ][3] = temp / ( temp + 1) * (1) / N;

    //The 4+5Nminus1..4+6Nminus1 states are the GT matrix
    //Only 2 entries can differ from 0, (N-1)G going to mono G and (N-1)T going to mono T
    //(N-1)G can only go to monomorphic state G, i.e. j=2
    //(N-1)G is at the end of the submatrix, i=4+6*(Nminus1)
  //  temp = (N-1)*(1+s[2]-s[3]);
    (*the_rate_matrix)[4 + 6*Nminus1 - 1][2] = computeEntryFromMoranProcessWithSelection(2, 3, Nminus1d);
  //  (*the_rate_matrix)[4 + 6*Nminus1 - 1][2] = temp / ( temp + 1) * (1) / N;
    //(N-1)T can only go to monomorphic state T, i.e. j=3
    //(N-1)T is at the begining of the submatrix, i=4
  //  temp = (N-1)*(1+s[3]-s[2]);
    (*the_rate_matrix)[4 + 5*Nminus1 ][3] = computeEntryFromMoranProcessWithSelection(3, 2, Nminus1d);
   // (*the_rate_matrix)[4 + 5*Nminus1 ][3] = temp / ( temp + 1) * (1) / N;

    
    //Now we need to fill the rest of the matrix, i.e. the B matrices along the diagonal.
    //In these B matrices, again most cells = 0.
    //The diagonal is such that it's 0 - (sum of the cells in the line)
    
    for (size_t k = 0; k <= 5; k++)
    {
        //Definition of the fitnesses
        double f1, f2;
        if (k<3)
        {
            f1 = s[0];
        }
        else if (k<5)
        {
            f1 = s[1];
        }
        else
        {
            f1 = s[2];
        }
        
        if (k==0)
        {
            f2 = s[1];
        }
        else if (k==1 || k==3)
        {
            f2 = s[2];
        }
        else
        {
            f2 = s[3];
        }
        
        for (size_t i = 1; i <= N-2 ; ++i)
        {
            size_t j = i+1;
            (*the_rate_matrix)[3+i+Nminus1*k][3+j+Nminus1*k] = (f1*i/(f1*i + f2*(N-i)) * (N-i)/N);
            (*the_rate_matrix)[3+j+Nminus1*k][3+i+Nminus1*k] = (f2*j/(f2*j + f1*(N-j)) * (N-j)/N);
        }

    }
    
    
    //In the first 4 rows/columns, the diagonal is defined such that the sum by line is 1.
   /* double sum = 0.0;
    for (size_t i=0; i< matrix_size; i++)
    {
        sum += (*the_rate_matrix)[0][i];
    }
    (*the_rate_matrix)[0][0] = 0-sum;
    
    sum = 0.0;
    for (size_t i=0; i< matrix_size; i++)
    {
        sum += (*the_rate_matrix)[1][i];
    }
    (*the_rate_matrix)[1][1] = 0-sum;
    
    sum = 0.0;
    for (size_t i=0; i< matrix_size; i++)
    {
        sum += (*the_rate_matrix)[2][i];
    }
    (*the_rate_matrix)[2][2] = 0-sum;
    
    sum = 0.0;
    for (size_t i=0; i< matrix_size; i++)
    {
        sum += (*the_rate_matrix)[3][i];
    }
    (*the_rate_matrix)[3][3] = 0-sum;*/
    
    /*
    for (size_t i=0; i< matrix_size; i++)
    {
        for (size_t j=0; j< matrix_size; j++)
        {
        (*the_rate_matrix)[i][j] *= (double) N;
        }
    }
    */
    
    // set the diagonal values
    setDiagonal();
    
    
    
    // rescale
    //rescaleToAverageRate( 1.0 );

    
    //Then we remove the identity matrix
  /*  for (size_t i=0; i< matrix_size; i++)
    {
        (*the_rate_matrix)[i][i] = (*the_rate_matrix)[i][i] - 1 ;
    }*/
}


double RateMatrix_PoMo::computeEntryFromMoranProcessWithSelection(size_t state1, size_t state2, double& count1)
{
    // We always assume state1 with count1 is increasing
    double count2 = (double)N-count1;
    
    // One of state2 alleles is chosen for disappearance
    double result = count2/(double)N; // 1/count2;
    
    // One of state1 alleles is chosen for replication
    result *= s[state1]*count1 / ( s[state2]*count2 + s[state1]*count1) ;
    return result;
}


/** Calculate the transition probabilities */
void RateMatrix_PoMo::calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const
{
    
    //Now the instantaneous rate matrix has been filled up entirely.
    //We use repeated squaring to quickly obtain exponentials, as in Poujol and Lartillot, Bioinformatics 2014.
    double t = rate * (startAge - endAge);
    computeExponentialMatrixByRepeatedSquaring(t, P);
    
    return;
}

void RateMatrix_PoMo::computeExponentialMatrixByRepeatedSquaring(double t,  TransitionProbabilityMatrix& P ) const
{
    //We use repeated squaring to quickly obtain exponentials, as in Poujol and Lartillot, Bioinformatics 2014.
    //Ideally one should dynamically decide how many squarings are necessary. 
    //For the moment, we arbitrarily do 10 such squarings, as it seems to perform well in practice (N. Lartillot, personal communication).
    //first, multiply the matrix by the right scalar
    //2^10 = 1024
    double tOver2s = t/(1024);
    for ( size_t i = 0; i < matrix_size; i++ )
    {
        for ( size_t j = 0; j < matrix_size; j++ )
        {
            P[i][j] = (*the_rate_matrix)[i][j] * tOver2s; 
        }
    }
    
    //Add the identity matrix:
     for ( size_t i = 0; i < matrix_size; i++ )
     {
         P[i][i] += 1;
     }
     //Now we can do the multiplications
     TransitionProbabilityMatrix P2 (matrix_size);
     squareMatrix (P, P2); //P2 at power 2
     squareMatrix (P2, P); //P at power 4
     squareMatrix (P, P2); //P2 at power 8
     squareMatrix (P2, P); //P at power 16
     squareMatrix (P, P2); //P2 at power 32
     squareMatrix (P2, P); //P at power 64
     squareMatrix (P, P2); //P2 at power 128
     squareMatrix (P2, P); //P at power 256
     squareMatrix (P, P2); //P2 at power 512
     squareMatrix (P2, P); //P at power 1024

     return;
}

inline void RateMatrix_PoMo::squareMatrix( TransitionProbabilityMatrix& P,  TransitionProbabilityMatrix& P2) const
{
    //Could probably use boost::ublas here, for the moment we do it ourselves.
    for ( size_t i = 0; i < matrix_size; i++ )
    {
        for ( size_t j = 0; j < matrix_size; j++ )
        {
            P2.getElement ( i, j ) = 0;
            for ( size_t k = 0; k < matrix_size; k++ )
            {
                P2.getElement ( i, j ) += P.getElement ( i, k ) * P.getElement ( k, j );
                
            }
        }
    }
}



RateMatrix_PoMo* RateMatrix_PoMo::clone( void ) const
{
    return new RateMatrix_PoMo( *this );
}

std::vector<double> RateMatrix_PoMo::getStationaryFrequencies( void ) const
{
    
    return stationary_freqs;
}


void RateMatrix_PoMo::update( void )
{
    
    if ( needs_update )
    {
        buildRateMatrix();
        // clean flags
        needs_update = false;
    }
}


void RateMatrix_PoMo::setMutationRates(const std::vector<double>& mr)
{

    size_t index = 0;
    for (size_t i=0; i<num_raw_states; ++i)
    {
        for (size_t j=0; j<num_raw_states; ++j)
        {
            if ( i!=j )
            {
                mu[i][j] = mr[index];
                ++index;
            }
        }
    }
    
}


void RateMatrix_PoMo::setMutationRates(const RateGenerator& mm)
{
    
    double age = 0.0;
    double rate = 1.0;
    
    for (size_t i=0; i<num_raw_states; ++i)
    {
        for (size_t j=0; j<num_raw_states; ++j)
        {
            if ( i!=j )
            {
                mu[i][j] = mm.getRate(i,j,age,rate);
            }
        }
    }
    
}


void RateMatrix_PoMo::setSelectionCoefficients(const std::vector<double>& sc)
{
    s = sc;

}
