

#ifndef RateMatrix_Jones_H
#define RateMatrix_Jones_H

#include "RateMatrix_Empirical.h"


namespace RevBayesCore {
    

/**
 * @brief Declaration of RateMatrix_Jones
 *
 * This contains the declaration of RateMatrix_Jones, which is
 * class that holds a rate matrix for a continuous-time Markov model.
 *
 * This class implements the Jones rate matrix. The JTT92 matrix has empirically derived exchangeability rates and stationary frequencies for the 20 amino acids.
 */
    class RateMatrix_Jones : public RateMatrix_Empirical {
        
    public:
        RateMatrix_Jones(void);                                                                                                //!< Construct rate matrix with n states
        virtual                             ~RateMatrix_Jones(void);                                                               //!< Destructor
                
        // RateMatrix functions
        RateMatrix_Jones*                   clone(void) const;
        
    private:
        
        
    };
    
}

#endif

