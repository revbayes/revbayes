#ifndef RateMatrix_Senca_H
#define RateMatrix_Senca_H

#include "TimeReversibleRateMatrix.h"
#include <complex>
#include <vector>


namespace RevBayesCore {
    
    class EigenSystem;
    class TransitionProbabilityMatrix;
    
    
    
    /**
     * @brief Codon rate matrix class with different rates for synonymous and non-synonymous substitutions.
     *
     * This class implements the special codon rate matrix.
     * The resulting rate matrix is computed by:
     *
     *		0: codon changes in more than one codon position (or stop codons)
     *		1: synonymous transition
     *		2: synonymous transversion
     *		3: non-synonymous transition
     *		4: non-synonymous transversion
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2015-07-15, version 1.0
     */
    class RateMatrix_Senca : public TimeReversibleRateMatrix {
        
    public:
        RateMatrix_Senca(size_t n);                                                                                               //!< Construct rate matrix with n states
        RateMatrix_Senca(const RateMatrix_Senca& m);                                                                                //!< Copy constructor
        virtual                             ~RateMatrix_Senca(void);                                                              //!< Destructor
        
        // overloaded operators
        RateMatrix_Senca&                   operator=(const RateMatrix_Senca& r);
        
        // RateMatrix functions
        void                                calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;    //!< Calculate the transition matrix
        RateMatrix_Senca*                   clone(void) const;
        void                                update(void);
        
    private:
        void                                calculateCijk(void);                                                                //!< Do precalculations on eigenvectors and their inverse
        double                              computePreferenceRatio(const CodonState &from, const CodonState &to) const;
        void                                computeOffDiagonal( void );
        void                                computeStateFrequencies( void );
        void                                tiProbsEigens(double t, TransitionProbabilityMatrix& P) const;                      //!< Calculate transition probabilities for real case
        void                                tiProbsComplexEigens(double t, TransitionProbabilityMatrix& P) const;               //!< Calculate transition probabilities for complex case
        void                                updateEigenSystem(void);                                                            //!< Update the system of eigenvalues and eigenvectors
        
        EigenSystem*                        theEigenSystem;                                                                     //!< Holds the eigen system
        std::vector<double>                 c_ijk;                                                                              //!< Vector of precalculated product of eigenvectors and their inverse
        std::vector<std::complex<double> >  cc_ijk;                                                                             //!< Vector of precalculated product of eigenvectors and thier inverse for complex case
        
        double                              kappa;
        double                              omega;
        
        std::vector<double>                 codon_freqs;
        
    };
    
}

#endif

