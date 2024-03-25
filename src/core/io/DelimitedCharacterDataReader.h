#ifndef DelimitedCharacterDataReader_H
#define DelimitedCharacterDataReader_H

#include <cstddef>
#include <vector>
#include <iosfwd>

#include "DelimitedDataReader.h"

namespace RevBayesCore {
    
    
    /**
     * Reader for character data matrices from delimited files.
     *
     * This reader is a simple file reader of a delimited file, e.g., by tab-stops.
     * In the first row should be the column names.
     * In the first column we expect the species names.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Will Freyman & Sebastian Hoehna)
     * @since 2015-03-03, version 1.0
     *
     */
    class DelimitedCharacterDataReader : public DelimitedDataReader {
        
    public:
        
        DelimitedCharacterDataReader(const std::string &fn, std::string d="", size_t ns=0);
        
        const std::vector<std::string>&                 getNames(void);
        const std::vector<std::vector<std::string> >&   getData(void);
        
        
        
    protected:
        
        std::vector<std::string>                        names;
        std::vector<std::vector<std::string> >          data;
        
    };
    
}

#endif
