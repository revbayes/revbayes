#ifndef RateMatrix_revPoMo2Nrecurrent_H
#define RateMatrix_revPoMo2Nrecurrent_H

#include "TimeReversibleRateMatrix.h"
#include <complex>
#include <vector>


namespace RevBayesCore {
    
    class EigenSystem;
    class TransitionProbabilityMatrix;
    
    /*
     */

    class RateMatrix_revPoMo2Nrecurrent : public TimeReversibleRateMatrix {
        
    public:
        RateMatrix_revPoMo2Nrecurrent(  long num_states,  long in_n );                                                                                            //!< Construct rate matrix with n states
        RateMatrix_revPoMo2Nrecurrent(const RateMatrix_revPoMo2Nrecurrent& m);                                                                  //!< Copy constructor
        virtual                                                ~RateMatrix_revPoMo2Nrecurrent(void);                                                              //!< Destructor
        
        // overloaded operators
        RateMatrix_revPoMo2Nrecurrent&                                   operator=(const RateMatrix_revPoMo2Nrecurrent& r);
        
        void                                                    calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;    //!< Calculate the transition matrix
        RateMatrix_revPoMo2Nrecurrent*                                   clone(void) const;

        void                                                    setN( long &ni );
        void                                                    setPi(  const Simplex &bf );
        void                                                    setRho( double &ex );

        void                                                    update(void);
        
    private:
        void                                                    calculateCijk(void);                                                                //!< Do precalculations on eigenvectors and their inverse
        void                                                    computeOffDiagonal( void );
        void                                                    tiProbsEigens(double t, TransitionProbabilityMatrix& P) const;                      //!< Calculate transition probabilities for real case
        void                                                    tiProbsComplexEigens(double t, TransitionProbabilityMatrix& P) const;               //!< Calculate transition probabilities for complex case
        void                                                    updateEigenSystem(void);                                                            //!< Update the system of eigenvalues and eigenvectors
        
        EigenSystem*                                            eigen_system;                                                                       //!< Holds the eigen system
        std::vector<double>                                     c_ijk;                                                                              //!< Vector of precalculated product of eigenvectors and their inverse
        std::vector<std::complex<double> >                      cc_ijk;                                                                             //!< Vector of precalculated product of eigenvectors and thier inverse for complex case
        
        long                                                    N;
        Simplex                                                 pi; 
        double                                                  rho;   
                
    };
    
}

#endif


