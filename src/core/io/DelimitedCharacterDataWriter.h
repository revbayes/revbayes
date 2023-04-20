#ifndef DelimitedCharacterDataWriter_H
#define DelimitedCharacterDataWriter_H

#include "AbstractHomologousDiscreteCharacterData.h"
#include "HomologousCharacterData.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

namespace RevBayesCore {
    
    /**
     * This class represents an object that writers character data objects into delimited text files.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Will Freyman)
     * @since 2013-04-15, version 1.0
     */
    class DelimitedCharacterDataWriter {
        
    public:
        DelimitedCharacterDataWriter();
        
        void                    writeData(const path& fn, const HomologousCharacterData &d, std::string del="\t");
//        void                    writeData(const path& fn, const Abstr &d, char del='\t');

        
    };
    
}


#endif
