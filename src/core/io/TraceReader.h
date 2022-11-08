#ifndef TraceReader_H
#define TraceReader_H

#include "Trace.h"
#include "RbFileManager.h"

#include <vector>
#include <string>

namespace RevBayesCore {
    
    
    /**
     * Reader for trace files.
     *
     * This reader is a reader of a trace files, e.g., tree-traces or stochastic variable traces.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-03-18, version 1.0
     *
     */
    class TraceReader {
        
    public:
//        TraceReader();
        
        std::vector<ModelTrace>             readStochasticVariableTrace( const path &fn, const std::string &delimiter );

        
    protected:
        
        
    };
    
}
#endif
