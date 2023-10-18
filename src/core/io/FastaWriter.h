#ifndef FastaWriter_H
#define FastaWriter_H

#include "AbstractHomologousDiscreteCharacterData.h"
#include "AbstractNonHomologousDiscreteCharacterData.h"
#include "RbFileManager.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

namespace RevBayesCore {
    
    /**
     * This class represents the writer object of character data objects into files in Fasta format.
     *
     * This class currently has only one functionality,
     * to write character data objects into a file in Fasta format.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2013-04-15, version 1.0
     */
    class FastaWriter {
        
    public:
        FastaWriter();
        
        void                    writeData(const path& fn, const AbstractHomologousDiscreteCharacterData& d);
        void                    writeData(const path& fn, const AbstractNonHomologousDiscreteCharacterData& d);
        
        
    };
    
}


#endif
