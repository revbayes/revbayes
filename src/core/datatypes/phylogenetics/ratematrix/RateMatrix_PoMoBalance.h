#ifndef RateMatrix_PoMoBalance_H
#define RateMatrix_PoMoBalance_H

#include "TimeReversibleRateMatrix.h"
#include <complex>
#include <vector>


namespace RevBayesCore {
    
    class EigenSystem;
    class TransitionProbabilityMatrix;
    
    /*
     */

    class RateMatrix_PoMoBalance : public TimeReversibleRateMatrix {
        
    public:
        RateMatrix_PoMoBalance( size_t ss, double in_n );                                                                                            //!< Construct rate matrix with n states
        RateMatrix_PoMoBalance(const RateMatrix_PoMoBalance& m);                                                                  //!< Copy constructor
        virtual                             ~RateMatrix_PoMoBalance(void);                                                              //!< Destructor
        
        // overloaded operators
        RateMatrix_PoMoBalance&             operator=(const RateMatrix_PoMoBalance& r);
        
        // RateMatrix functions
        virtual RateMatrix_PoMoBalance&                   assign(const Assignable &m);                                                                                            //!< Assign operation that can be called on a base class instance.
        void                                                    calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;    //!< Calculate the transition matrix
        RateMatrix_PoMoBalance*                           clone(void) const;

        void                                                    setN( double n );
        void                                                    setPi(const std::vector<double> &p );
        void                                                    setRho( const std::vector<double> &r );
        void                                                    setSigma( const std::vector<double> &s );
        void                                                    setBeta( double b );

        void                                                    update(void);
        
    private:
        void                                                    calculateCijk(void);                                                                //!< Do precalculations on eigenvectors and their inverse
        void                                                    computeOffDiagonal( void );
        double                                                  midpoint_beta( double n, int indicator ) const;
        void                                                    tiProbsEigens(double t, TransitionProbabilityMatrix& P) const;                      //!< Calculate transition probabilities for real case
        void                                                    tiProbsComplexEigens(double t, TransitionProbabilityMatrix& P) const;               //!< Calculate transition probabilities for complex case
        void                                                    updateEigenSystem(void);                                                            //!< Update the system of eigenvalues and eigenvectors
        
        EigenSystem*                                            eigen_system;                                                                       //!< Holds the eigen system
        std::vector<double>                                     c_ijk;                                                                              //!< Vector of precalculated product of eigenvectors and their inverse
        std::vector<std::complex<double> >                      cc_ijk;                                                                             //!< Vector of precalculated product of eigenvectors and thier inverse for complex case
        
        double                                                  N;
        std::vector<double>                                     pi;
        std::vector<double>                                     rho;
        std::vector<double>                                     sigma;   
        double                                                  beta;                                                                //!< Vector of precalculated product of eigenvectors and their inverse
                
    };
    
}

#endif