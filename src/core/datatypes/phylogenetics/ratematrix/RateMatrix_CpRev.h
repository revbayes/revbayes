#ifndef RateMatrix_CpRev_H
#define RateMatrix_CpRev_H

#include "RateMatrix_Empirical.h"


namespace RevBayesCore {
    
    class TransitionProbabilityMatrix;
    
    class RateMatrix_CpRev : public RateMatrix_Empirical {
        
    public:
        RateMatrix_CpRev(void);                                                                                                 //!< Construct rate matrix with n states
        virtual                            ~RateMatrix_CpRev(void);                                                             //!< Destructor
        
        // RateMatrix functions
        RateMatrix_CpRev*                   clone(void) const;
        
    private:
        
        
    };
    
}

#endif

