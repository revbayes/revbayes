#ifndef ProgressBar_H
#define ProgressBar_H

#include <cstddef>

#include "Parallelizable.h"

namespace RevBayesCore {
    
    
    /**
     * Show progress on terminal.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2016-10-07, version 1.0
     *
     */
    class ProgressBar : public Parallelizable {
        
    public:
        ProgressBar(size_t m, size_t o = 0);
        
        
        void                            finish(void);
        void                            start(void);
        void                            update(size_t i);
        
        
        
    private:
        
        size_t                          max;
        size_t                          num_stars;
        size_t                          offset;
        
        
        
    };
    
}
#endif
