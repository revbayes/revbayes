#ifndef RateMatrix_CodonSynonymousNonsynonymous_H
#define RateMatrix_CodonSynonymousNonsynonymous_H

#include "TimeReversibleRateMatrix.h"
#include <complex>
#include <vector>


namespace RevBayesCore {
    
    class EigenSystem;
    class TransitionProbabilityMatrix;
    
    
    
    /**
     * @brief Codon rate matrix class with different rates for synonymous and non-synonymous substitutions.
     *
     * This class implements the codon rate matrix with different dN/ dS rates.
     * The only parameter is the stationary codon frequencies.
     * The resulting rate matrix is computed by:
     *
     *		0: codon changes in more than one codon position (or stop codons)
     *		1: synonymous substitution * pi(target)
     *		2: non-synonymous substitution * pi(target)
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2015-07-15, version 1.0
     */
    class RateMatrix_CodonSynonymousNonsynonymous : public TimeReversibleRateMatrix {
        
    public:
        RateMatrix_CodonSynonymousNonsynonymous(void);                                                                                               //!< Construct rate matrix with n states
        RateMatrix_CodonSynonymousNonsynonymous(const RateMatrix_CodonSynonymousNonsynonymous& m);                                                                                //!< Copy constructor
        virtual                             ~RateMatrix_CodonSynonymousNonsynonymous(void);                                                              //!< Destructor
        
        // overloaded operators
        RateMatrix_CodonSynonymousNonsynonymous&                operator=(const RateMatrix_CodonSynonymousNonsynonymous& r);
        
        // RateMatrix functions
        void                                                    calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;    //!< Calculate the transition matrix
        RateMatrix_CodonSynonymousNonsynonymous*                clone(void) const;
        void                                                    setCodonFrequencies(const std::vector<double> &f);                                 //!< Set the nucleotide frequencies
        void                                                    setOmega(double o);
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
        
        double                                                  omega;
        std::vector<double>                                     codon_freqs; 
        
    };
    
}

#endif

