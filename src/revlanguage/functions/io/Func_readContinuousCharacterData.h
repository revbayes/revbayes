#ifndef Func_readContinuousCharacterData_H
#define Func_readContinuousCharacterData_H

#include "Procedure.h"
#include "RbFileManager.h"

#include <map>
#include <string>
#include <vector>


namespace RevLanguage {
    
    /**
     * This is the interface for a class that reads in phylogenetic continuous character matrix data.
     * It uses the Nexus Class Library (NCL) to do the heavy lifting
     * for reading files in Fasta, Phylip, or Nexus formats. The function will take in a string (RbString)
     * as an argument that indicates the file location or directory location. If the path is to a file, then
     * that file is read (or an attempt to read it is made). If the path is to a directory, then the program
     * attempts to read _all_ of the files in that directory, and does so recursively, burrowing into other
     * directories in that path and reading the files there too. The function attempts to figure out the file
     * format type (Nexus, Fasta, or Phylip), whether the data is interleaved or not (if the file is in Phylip
     * format), and the data type (if the file is in Fasta or Phylip format). Currently, the program can read the
     * following:
     *
     *            A   B
     *    Nexus   Y   Y   
     *    Phylip  N   N   
     *    Fasta   N   N
     *
     *    A: Continuous data (not interleaved);
     *    B: Continuous data (interleaved);
     *
     * Note that some data formats do not support interleaving of data (such as Fasta) or reading of certain data types, so
     * the table is unlikely to complete with Yes's (Y) in all of the cells.
     */
    
    class Func_readContinuousCharacterData : public Procedure {
        
    public:
        // Basic utility functions
        Func_readContinuousCharacterData*       clone(void) const;                                          //!< Clone the object
        static const std::string&               getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                  getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                             getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                         getTypeSpec(void) const;                                    //!< Get language type of the object
        
        // Regular functions
        RevPtr<RevVariable>                     execute(void);                                              //!< Execute function
        const ArgumentRules&                    getArgumentRules(void) const;                               //!< Get argument rules
        const TypeSpec&                         getReturnType(void) const;                                  //!< Get type of return value
                
    };
    
}

#endif

