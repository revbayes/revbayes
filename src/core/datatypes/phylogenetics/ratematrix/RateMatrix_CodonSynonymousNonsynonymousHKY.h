#ifndef RateMatrix_CodonSynonymousNonsynonymousHKY_H
#define RateMatrix_CodonSynonymousNonsynonymousHKY_H

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
     * The JC has no parameter but can be applied to any number of states.
     * The resulting rate matrix is computed by:
     *
     *        0: codon changes in more than one codon position (or stop codons)
     *        1: synonymous transition
     *        2: synonymous transversion
     *        3: non-synonymous transition
     *        4: non-synonymous transversion
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2015-07-15, version 1.0
     */
    class RateMatrix_CodonSynonymousNonsynonymousHKY : public TimeReversibleRateMatrix {
        
    public:
        RateMatrix_CodonSynonymousNonsynonymousHKY(void);                                                                                               //!< Construct rate matrix with n states
        RateMatrix_CodonSynonymousNonsynonymousHKY(const RateMatrix_CodonSynonymousNonsynonymousHKY& m);                                                                                //!< Copy constructor
        virtual                             ~RateMatrix_CodonSynonymousNonsynonymousHKY(void);                                                              //!< Destructor
        
        // overloaded operators
        RateMatrix_CodonSynonymousNonsynonymousHKY&             operator=(const RateMatrix_CodonSynonymousNonsynonymousHKY& r);
        
        // RateMatrix functions
        void                                                    calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;    //!< Calculate the transition matrix
        RateMatrix_CodonSynonymousNonsynonymousHKY*             clone(void) const;
        void                                                    setNucleotideFrequencies(const std::vector<double> &f);                             //!< Set the nucleotide frequencies
        void                                                    setKappa(double k);
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
        
        double                                                  kappa;
        double                                                  omega;
        std::vector<double>                                     nucleotide_freqs;                                                                              //!< Vector of precalculated product of eigenvectors and their inverse
                
    };
    
}

#endif


