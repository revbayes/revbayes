#ifndef Func_readFossilCountsFile_hpp
#define Func_readFossilCountsFile_hpp


#include "Procedure.h"
#include "RbFileManager.h"

#include <map>
#include <string>
#include <vector>


namespace RevLanguage {
    
    /**
     * The Rev procedure to read fossil counts from delimited files.
     * 
     * 
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Bruno do Rosario Petrucci)
     * @since 2025-02-13, version 1.2.5
     *
     */
    class Func_readFossilCountsFile : public Procedure {
        
    public:
        // Basic utility functions
        Func_readFossilCountsFile*          clone(void) const;                                          //!< Clone the object
        static const std::string&           getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&              getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                         getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                     getTypeSpec(void) const;                                    //!< Get language type of the object
        
        // Regular functions
        RevPtr<RevVariable>                 execute(void);                                              //!< Execute function
        const ArgumentRules&                getArgumentRules(void) const;                               //!< Get argument rules
        const TypeSpec&                     getReturnType(void) const;                                  //!< Get type of return value
        
    private:
        std::string                         bitToState(const std::string &s);
        
    };
    
}

#endif /* Func_readFossilCountsFile_hpp */
