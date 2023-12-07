#ifndef TraceContinuousReader_H
#define TraceContinuousReader_H

#include "DelimitedDataReader.h"
#include "TraceNumeric.h"

namespace RevBayesCore {
    
    
    
    /**
     * Reader for continuous traces of MCMC samples.
     *
     * This reader is a simple file reader of a delimited file, e.g., by tab-stops.
     * In the first column should be the traces names and in the second column the data of each sample.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2015-03-02, version 1.0
     *
     */
    class TraceContinuousReader : public DelimitedDataReader {
        
    public:
        
        TraceContinuousReader(const path &fn, std::string d="");              //!< Constructor
        
        std::vector<TraceNumeric>&          getTraces(void);                            //!< Get the data.
        const std::vector<TraceNumeric>&    getTraces(void) const;                      //!< Get the data.
        
        
    protected:
        
        std::vector<TraceNumeric>           data;
    };
}

#endif
