#ifndef RateMatrix_revPoMoBalance4N_H
#define RateMatrix_revPoMoBalance4N_H

#include "TimeReversibleRateMatrix.h"
#include <complex>
#include <vector>


namespace RevBayesCore {
    
    class EigenSystem;
    class TransitionProbabilityMatrix;
    
    /*
     */

    class RateMatrix_revPoMoBalance4N : public TimeReversibleRateMatrix {
        
    public:
        RateMatrix_revPoMoBalance4N( std::int64_t ss, std::int64_t in_n );                                                                                            //!< Construct rate matrix with n states
        RateMatrix_revPoMoBalance4N(const RateMatrix_revPoMoBalance4N& m);                                                                  //!< Copy constructor
        virtual                             ~RateMatrix_revPoMoBalance4N(void);                                                              //!< Destructor
        
        // overloaded operators
        RateMatrix_revPoMoBalance4N&                                 operator=(const RateMatrix_revPoMoBalance4N& r);
        
        // RateMatrix functions
        void                                                    calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;    //!< Calculate the transition matrix
        RateMatrix_revPoMoBalance4N*                                 clone(void) const;

        void                                                    setN( std::int64_t n );
        void                                                    setPi(const std::vector<double> &p );
        void                                                    setRho( const std::vector<double> &r );
        void                                                    setPhi( const std::vector<double> &s );
        void                                                    setBeta( const std::vector<double> &b );
        void                                                    setB( const std::vector<std::int64_t> &Bf );

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
        
        std::int64_t                                                     N;
        std::vector<double>                                     pi;
        std::vector<double>                                     rho;
        std::vector<double>                                     phi;   
        std::vector<double>                                     beta;                                                                //!< Vector of precalculated product of eigenvectors and their inverse
        std::vector<std::int64_t>                                        B;
    };
    
}

#endif