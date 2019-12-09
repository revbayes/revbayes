#ifndef RateMatrix_PoMoThree_H
#define RateMatrix_PoMoThree_H

#include "TimeReversibleRateMatrix.h"
#include <complex>
#include <vector>


namespace RevBayesCore {
    
    class EigenSystem;
    class TransitionProbabilityMatrix;
    
    
    
    class RateMatrix_PoMoThree : public TimeReversibleRateMatrix {
        
    public:
        RateMatrix_PoMoThree(void);                                                                                               //!< Construct rate matrix with n states
        RateMatrix_PoMoThree(const RateMatrix_PoMoThree& m);                                                                                //!< Copy constructor
        virtual                             ~RateMatrix_PoMoThree(void);                                                              //!< Destructor
        
        // overloaded operators
        RateMatrix_PoMoThree&                          operator=(const RateMatrix_PoMoThree& r);
        
        // RateMatrix functions
        virtual RateMatrix_PoMoThree&                  assign(const Assignable &m);                                                                                            //!< Assign operation that can be called on a base class instance.
        void                                           calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;    //!< Calculate the transition matrix
        RateMatrix_PoMoThree*                          clone(void) const;

		void 										   setPopulationSize(double ps);
		void 										   setEquilibriumFrequencies( const std::vector<double> &f );
		void 										   setExchangeabilities( const std::vector<double> &r );
		void 										   setSelectionCoefficients( const std::vector<double> &s );
        void                                           update(void);
        
    private:
        void                                            calculateCijk(void);                                                                //!< Do precalculations on eigenvectors and their inverse
        void                                            computeOffDiagonal( void );
		std::vector<double> 						    transformSigma( void ) const ;
		std::vector<double> 							transformRho( void ) const ;
        std::vector<double>                             getStationaryFrequencies( void ) const ;
        double                                          getExpectedNumberEvents( std::vector<double> rho2, std::vector<double> sigma2 ) const;


        void                                            tiProbsEigens(double t, TransitionProbabilityMatrix& P) const;                      //!< Calculate transition probabilities for real case
        void                                            tiProbsComplexEigens(double t, TransitionProbabilityMatrix& P) const;               //!< Calculate transition probabilities for complex case
        void                                            updateEigenSystem(void);                                                            //!< Update the system of eigenvalues and eigenvectors
        
        EigenSystem*                                    eigen_system;                                                                       //!< Holds the eigen system
        std::vector<double>                             c_ijk;                                                                              //!< Vector of precalculated product of eigenvectors and their inverse
        std::vector<std::complex<double> >              cc_ijk;                                                                             //!< Vector of precalculated product of eigenvectors and thier inverse for complex case
        
        double                                          N;
        std::vector<double>                             pi;                                                                                 //!< Vector of precalculated product of eigenvectors and their inverse
        std::vector<double>                             rho;                                                                                //!< Vector of precalculated product of eigenvectors and their inverse
        std::vector<double>                             sigma;                                                                              //!< Vector of precalculated product of eigenvectors and their inverse
        
    };
    
}

#endif


