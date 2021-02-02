#ifndef RlAbstractCharacterData_H
#define RlAbstractCharacterData_H

#include <vector>
#include <iosfwd>

#include "MethodTable.h"
#include "RevPtr.h"

namespace RevBayesCore { class AbstractCharacterData; }


namespace RevLanguage {
class Argument;
class RevVariable;
    
    
    class AbstractCharacterData {
        
    public:
    
        AbstractCharacterData(void);                                                                                                 //!< Constructor requires character type
        
        MethodTable                                             getCharacterDataMethods(void) const;                                                                            //!< Get the common character data methods
        RevPtr<RevVariable>                                     executeCharacterDataMethod(const std::string& name, const std::vector<Argument>& args, bool &found,  RevBayesCore::AbstractCharacterData *charDataObject);            //!< Override to map member methods to internal functions
        
    };
    
}


#endif

