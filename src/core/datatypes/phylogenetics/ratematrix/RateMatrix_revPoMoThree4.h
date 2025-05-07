#ifndef RateMatrix_revPoMoThree4_H
#define RateMatrix_revPoMoThree4_H

#include "TimeReversibleRateMatrix.h"
#include <complex>
#include <vector>


namespace RevBayesCore {
    
    class EigenSystem;
    class TransitionProbabilityMatrix;
    
    
    
    class RateMatrix_revPoMoThree4 : public TimeReversibleRateMatrix {
        
    public:
        RateMatrix_revPoMoThree4(void);                                                                                                       //!< Construct rate matrix with n states
        RateMatrix_revPoMoThree4(const RateMatrix_revPoMoThree4& m);                                                                                //!< Copy constructor
        virtual                             ~RateMatrix_revPoMoThree4(void);                                                                  //!< Destructor
        
        // overloaded operators
        RateMatrix_revPoMoThree4&                     operator=(const RateMatrix_revPoMoThree4& r);
        
        // RateMatrix functions
        void                                           calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;    //!< Calculate the transition matrix
        RateMatrix_revPoMoThree4*                     clone(void) const;

		
		void 										   setPi( const std::vector<double> &f );
		void 										   setRho( const std::vector<double> &r );
        void                                           setPhi( const std::vector<double> &fc );

        void                                           update(void);
        
    private:
        void                                            calculateCijk(void);                                                                //!< Do precalculations on eigenvectors and their inverse
        void                                            computeOffDiagonal( void );

        void                                            tiProbsEigens(double t, TransitionProbabilityMatrix& P) const;                      //!< Calculate transition probabilities for real case
        void                                            tiProbsComplexEigens(double t, TransitionProbabilityMatrix& P) const;               //!< Calculate transition probabilities for complex case
        void                                            updateEigenSystem(void);                                                            //!< Update the system of eigenvalues and eigenvectors
        
        EigenSystem*                                    eigen_system;                                                                       //!< Holds the eigen system
        std::vector<double>                             c_ijk;                                                                              //!< Vector of precalculated product of eigenvectors and their inverse
        std::vector<std::complex<double> >              cc_ijk;                                                                             //!< Vector of precalculated product of eigenvectors and thier inverse for complex case
        
        std::vector<double>                             pi;                                                                                 //!< Vector of precalculated product of eigenvectors and their inverse
        std::vector<double>                             rho;                                                                                //!< Vector of precalculated product of eigenvectors and their inverse
        std::vector<double>                             phi;                                                                                //!< Vector of precalculated product of eigenvectors and their inverse
       
    };
    
}

#endif


