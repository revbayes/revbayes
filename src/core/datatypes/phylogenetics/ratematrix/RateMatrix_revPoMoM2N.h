#ifndef RateMatrix_revPoMoM2N_H
#define RateMatrix_revPoMoM2N_H

#include "TimeReversibleRateMatrix.h"
#include <complex>
#include <vector>


namespace RevBayesCore {
    
    class EigenSystem;
    class TransitionProbabilityMatrix;
    
    
    
    class RateMatrix_revPoMoM2N : public TimeReversibleRateMatrix {
        
    public:
        RateMatrix_revPoMoM2N(size_t ss );                                                                                                       //!< Construct rate matrix with n states
        RateMatrix_revPoMoM2N(const RateMatrix_revPoMoM2N& m);                                                                                //!< Copy constructor
        virtual                             ~RateMatrix_revPoMoM2N(void);                                                                  //!< Destructor
        
        // overloaded operators
        RateMatrix_revPoMoM2N&                            operator=(const RateMatrix_revPoMoM2N& r);
        
        // RateMatrix functions
        virtual RateMatrix_revPoMoM2N&               assign(const Assignable &m);                                                                                            //!< Assign operation that can be called on a base class instance.
        void                                           calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;    //!< Calculate the transition matrix
        RateMatrix_revPoMoM2N*                       clone(void) const;

		
		void 										   setN(long ps);
		void 										   setMu( const std::vector<double> &r );
        void                                           setM(long vps);
        void                                           setGen(double g);

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
        
        long                                            N;
        std::vector<double>                             mu;                                                                                //!< Vector of precalculated product of eigenvectors and their inverse
        long                                            M;
        double                                          gen;

    };
    
}

#endif


