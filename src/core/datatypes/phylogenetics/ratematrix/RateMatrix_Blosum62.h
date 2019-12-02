#ifndef RateMatrix_Blosum62_H
#define RateMatrix_Blosum62_H

#include "RateMatrix_Empirical.h"


namespace RevBayesCore {
    
    class TransitionProbabilityMatrix;
    
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

