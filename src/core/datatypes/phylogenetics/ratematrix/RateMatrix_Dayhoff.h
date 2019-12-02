#ifndef RateMatrix_Dayhoff_H
#define RateMatrix_Dayhoff_H

#include "RateMatrix_Empirical.h"


namespace RevBayesCore {
    
    class TransitionProbabilityMatrix;
    
    class RateMatrix_Dayhoff : public RateMatrix_Empirical {
        
    public:
        RateMatrix_Dayhoff(void);                                                                                                   //!< Construct rate matrix with n states
        virtual                             ~RateMatrix_Dayhoff(void);                                                              //!< Destructor
                
        // RateMatrix functions
        RateMatrix_Dayhoff*                 clone(void) const;
        
    private:
        
        
    };
    
}

#endif

