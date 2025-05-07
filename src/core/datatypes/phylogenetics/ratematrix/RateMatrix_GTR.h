#ifndef RateMatrix_GTR_H
#define RateMatrix_GTR_H

#include "TimeReversibleRateMatrix.h"
#include "EigenSystem.h"
#include <complex>
#include <vector>


namespace RevBayesCore {
    
//    class EigenSystem;
    class TransitionProbabilityMatrix;
    
    /**
     * @brief Generalised Time Reversible Model (GTR)rate matrix class (Tavare 1986).
     *
     * This class implements the special GTR rate matrix with the known analytical solution for the transition probabilities.
     * The GTR matrix has 6 transition probabilities and the four base frequency parameters.
     *
     */
    class RateMatrix_GTR : public TimeReversibleRateMatrix {
        
    public:
        RateMatrix_GTR(size_t n);                                                                                               //!< Construct rate matrix with n states
        RateMatrix_GTR(const RateMatrix_GTR& m);                                                                                //!< Copy constructor
        virtual                             ~RateMatrix_GTR(void);                                                              //!< Destructor
        
        // overloaded operators
        RateMatrix_GTR&                     operator=(const RateMatrix_GTR& r);
        
        // RateMatrix functions
        void                                calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;    //!< Calculate the transition matrix
        RateMatrix_GTR*                     clone(void) const;
        void                                update(void);
        virtual void                        initFromString( const std::string &s );                                             //!< Serialize (resurrect) the object from a string value
        
    private:
        void                                calculateCijk(void);                                                                //!< Do precalculations on eigenvectors and their inverse
        void                                tiProbsEigens(double t, TransitionProbabilityMatrix& P) const;                      //!< Calculate transition probabilities for real case
        void                                tiProbsComplexEigens(double t, TransitionProbabilityMatrix& P) const;               //!< Calculate transition probabilities for complex case
        void                                updateEigenSystem(void);                                                            //!< Update the system of eigenvalues and eigenvectors
        
        EigenSystem*                        theEigenSystem;                                                                     //!< Holds the eigen system
        std::vector<double>                 c_ijk;                                                                              //!< Vector of precalculated product of eigenvectors and their inverse
        std::vector<std::complex<double> >  cc_ijk;                                                                             //!< Vector of precalculated product of eigenvectors and thier inverse for complex case
        
                
    };
        
}

#endif

