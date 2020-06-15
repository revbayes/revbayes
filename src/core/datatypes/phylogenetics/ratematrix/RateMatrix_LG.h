#ifndef RateMatrix_LG_H
#define RateMatrix_LG_H

#include "RateMatrix_Empirical.h"


namespace RevBayesCore {

/**
 * @brief
 * This file contains the declaration of RateMatrix_LG, which is
 * class that holds a rate matrix for a continuous-time Markov model.
 *
 * This is an empirically derived rate matrix for the 20 amino acids obtained from 3,912 Pfam alignments (Le and Gascuel 2008)
 *
 */
    class RateMatrix_LG : public RateMatrix_Empirical {
        
    public:
        RateMatrix_LG(void);                                                                                                //!< Construct rate matrix with n states
        virtual                             ~RateMatrix_LG(void);                                                               //!< Destructor
        
        // RateMatrix functions
        RateMatrix_LG*                     clone(void) const;
        
    private:
        
        
    };
    
}

#endif

