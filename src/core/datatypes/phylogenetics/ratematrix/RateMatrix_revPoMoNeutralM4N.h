#ifndef RateMatrix_revPoMoNeutralM4N_H
#define RateMatrix_revPoMoNeutralM4N_H

#include "TimeReversibleRateMatrix.h"
#include <complex>
#include <vector>


namespace RevBayesCore {
    
    class EigenSystem;
    class TransitionProbabilityMatrix;
    
    
    
    class RateMatrix_revPoMoNeutralM4N : public TimeReversibleRateMatrix {
        
    public:
        RateMatrix_revPoMoNeutralM4N( std::int64_t num_states, std::int64_t in_m );                                                                                                       //!< Construct rate matrix with n states
        RateMatrix_revPoMoNeutralM4N(const RateMatrix_revPoMoNeutralM4N& m);                                                                                //!< Copy constructor
        virtual                             ~RateMatrix_revPoMoNeutralM4N(void);                                                                  //!< Destructor
        
        // overloaded operators
        RateMatrix_revPoMoNeutralM4N&                            operator=(const RateMatrix_revPoMoNeutralM4N& r);
        
        // RateMatrix functions
        void                                           calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;    //!< Calculate the transition matrix
        RateMatrix_revPoMoNeutralM4N*                       clone(void) const;

		
		void 										   setN(std::int64_t ps);
        void                                           setM(std::int64_t vps);
		void 										   setPi( const std::vector<double> &f );
		void 										   setRho( const std::vector<double> &r );
        std::vector<double>                            getStationaryFrequencies( void ) const;

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
        
        std::int64_t                                            N;
        std::int64_t                                            M;
        std::vector<double>                             pi;                                                                                 //!< Vector of precalculated product of eigenvectors and their inverse
        std::vector<double>                             rho;                                                                                //!< Vector of precalculated product of eigenvectors and their inverse
        
    };
    
}

#endif


