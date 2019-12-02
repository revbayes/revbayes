#ifndef RlAbstractCharacterData_H
#define RlAbstractCharacterData_H

#include "AbstractCharacterData.h"
#include "ModelObject.h"

#include <set>
#include <string>
#include <vector>


namespace RevLanguage {
    
    
    class AbstractCharacterData {
        
    public:
    
        AbstractCharacterData(void);                                                                                                 //!< Constructor requires character type
        
        MethodTable                                             getCharacterDataMethods(void) const;                                                                            //!< Get the common character data methods
        RevPtr<RevVariable>                                     executeCharacterDataMethod(const std::string& name, const std::vector<Argument>& args, bool &found,  RevBayesCore::AbstractCharacterData *charDataObject);            //!< Override to map member methods to internal functions
        
    };
    
}


#endif

