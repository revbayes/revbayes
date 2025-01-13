#ifndef RateMatrix_FlowT2Populations_H
#define RateMatrix_FlowT2Populations_H

#include "TimeReversibleRateMatrix.h"
#include <complex>
#include <vector>


namespace RevBayesCore {
    
    class EigenSystem;
    class TransitionProbabilityMatrix;
    
    /*
     */

    class RateMatrix_FlowT2Populations : public TimeReversibleRateMatrix {
        
    public:
        RateMatrix_FlowT2Populations( size_t ss, const std::vector<double>& in_nu );                                                                                            //!< Construct rate matrix with n states
        RateMatrix_FlowT2Populations(const RateMatrix_FlowT2Populations& m);                                                                  //!< Copy constructor
        virtual                             ~RateMatrix_FlowT2Populations(void);                                                              //!< Destructor
        
        // overloaded operators
        RateMatrix_FlowT2Populations&             operator=(const RateMatrix_FlowT2Populations& r);
        
        // RateMatrix functions
        void                                                    calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;    //!< Calculate the transition matrix
        RateMatrix_FlowT2Populations*                           clone(void) const;

        void                                                    setNu(const std::vector<double> &n );
        void                                                    setMu(const std::vector<double> &m );
        void                                                    setLambda( const std::vector<double> &l );
        void                                                    setSigma( const std::vector<double> &s );

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
        
        std::vector<double>                                     nu;
        std::vector<double>                                     mu;
        std::vector<double>                                     lambda;
        std::vector<double>                                     sigma;                                                                           //!< Vector of precalculated product of eigenvectors and their inverse
                
    };
    
}

#endif


