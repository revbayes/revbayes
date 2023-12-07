#ifndef DelimitedDataReader_H
#define DelimitedDataReader_H

#include <cstddef>
#include <vector>
#include <iosfwd>
#include <string> // IWYU pragma: keep
#include "RbFileManager.h"

namespace RevBayesCore {
    
    
    /**
     * Reader for delimited files.
     *
     * This reader is a simple file reader of a delimited files, e.g., by tab-stops or ';'.
     * That is, all files are in table format.
     * This reader provides a simple base class for other specific reader that know what data type
     * is in each column and can covert this.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Michael Landis)
     * @since 2014-03-18, version 1.0
     *
     */
    class DelimitedDataReader {
        
    public:
        DelimitedDataReader(const path &fn, std::string d="", size_t ns=0);
        
        void                                                readData( size_t ls);
        const std::vector<std::vector<std::string> >&       getChars(void);
        const path&                                         getFilename(void);
        
    
    protected:
        
        // protected member only accessible for derived classes
        path                                                filename;
        std::string                                         delimiter;
        std::vector<std::vector<std::string> >              chars;
        
    };
    
}
#endif 
