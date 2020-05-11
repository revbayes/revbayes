#ifndef RateMatrix_Blosum62_H
#define RateMatrix_Blosum62_H

#include "RateMatrix_Empirical.h"

/*
 * @file
 * This file contains the declaration of Blosum62 rate matrix, which is a
 * class that holds a rate matrix for a continuous-time Markov model
 *
 * @brief Declaration of RateMatrix_Blosum62.h
 *
 * This class implements the Blosum62 rate matrix. The Blosum62 matrix has empirically derived exchangeability rates and stationary frequencies for the 20 amino acids.
 *
 *
 */


namespace RevBayesCore {
    
    class RateMatrix_Blosum62 : public RateMatrix_Empirical {
        
    public:
        RateMatrix_Blosum62(void);                                                                                                  //!< Construct rate matrix with n states
        virtual                            ~RateMatrix_Blosum62(void);                                                              //!< Destructor
        
        // RateMatrix functions
        RateMatrix_Blosum62*                clone(void) const;
        
    private:
        
        
    };
    
}

#endif

