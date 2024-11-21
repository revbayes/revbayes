#ifndef RateMatrix_revPoMo2N_H
#define RateMatrix_revPoMo2N_H

#include "TimeReversibleRateMatrix.h"
#include <complex>
#include <cstdint>
#include <vector>


namespace RevBayesCore {
    
    class EigenSystem;
    class TransitionProbabilityMatrix;
    
    /*
     */

    class RateMatrix_revPoMo2N : public TimeReversibleRateMatrix {
        
    public:
        RateMatrix_revPoMo2N(  std::int64_t num_states, std::int64_t in_n );                                                                                            //!< Construct rate matrix with n states
        RateMatrix_revPoMo2N(const RateMatrix_revPoMo2N& m);                                                                  //!< Copy constructor
        virtual                             ~RateMatrix_revPoMo2N(void);                                                              //!< Destructor
        
        // overloaded operators
        RateMatrix_revPoMo2N&                                   operator=(const RateMatrix_revPoMo2N& r);
        
        // RateMatrix functions
        void                                                    calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;    //!< Calculate the transition matrix
        RateMatrix_revPoMo2N*                                   clone(void) const;

        void                                                    setN( std::int64_t &ni );
        void                                                    setPi(  const Simplex &bf );
        void                                                    setRho( double &ex );
        void                                                    setPhi( const std::vector<double> &f );
        std::vector<double>                                     getStationaryFrequencies( void ) const;

        void                                                    update(void);
        
    private:
        void                                                    calculateCijk(void);                                                                //!< Do precalculations on eigenvectors and their inverse
        double                                                  calculateReciprocalExpectedDivergence( void );
        void                                                    computeOffDiagonal( void );
        void                                                    tiProbsEigens(double t, TransitionProbabilityMatrix& P) const;                      //!< Calculate transition probabilities for real case
        void                                                    tiProbsComplexEigens(double t, TransitionProbabilityMatrix& P) const;               //!< Calculate transition probabilities for complex case
        void                                                    updateEigenSystem(void);                                                            //!< Update the system of eigenvalues and eigenvectors
        
        EigenSystem*                                            eigen_system;                                                                       //!< Holds the eigen system
        std::vector<double>                                     c_ijk;                                                                              //!< Vector of precalculated product of eigenvectors and their inverse
        std::vector<std::complex<double> >                      cc_ijk;                                                                             //!< Vector of precalculated product of eigenvectors and thier inverse for complex case
        
        std::int64_t                                                    N;
        Simplex                                                 pi; 
        double                                                  rho;   
        std::vector<double>                                     phi;                                                                           //!< Vector of precalculated product of eigenvectors and their inverse
                
    };
    
}

#endif


