/**
 * @file
 * This file contains the declaration of RateMatrix_JC, which is
 * class that holds a rate matrix for a continuous-time Markov model.
 *
 * @brief Declaration of RateMatrix_JC
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date: 2012-12-05 12:20:37 -0800 (Wed, 05 Dec 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 *
 * $Id: RateMatrix.h 1901 2012-12-05 20:20:37Z boussau $
 */

#ifndef RateMatrix_MtRev24_H
#define RateMatrix_MtRev24_H

#include "RateMatrix_Empirical.h"


namespace RevBayesCore {
    
    class TransitionProbabilityMatrix;
    
    class RateMatrix_MtRev24 : public RateMatrix_Empirical {
        
    public:
        RateMatrix_MtRev24(void);                                                                                                //!< Construct rate matrix with n states
        virtual                             ~RateMatrix_MtRev24(void);                                                               //!< Destructor
        
        // RateMatrix functions
        RateMatrix_MtRev24*                 clone(void) const;
        
    private:
        
        
    };
    
}

#endif

