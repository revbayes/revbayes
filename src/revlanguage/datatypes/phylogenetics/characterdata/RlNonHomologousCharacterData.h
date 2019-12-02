#ifndef RlNonHomologousCharacterData_H
#define RlNonHomologousCharacterData_H

#include "NonHomologousCharacterData.h"
#include "RlAbstractCharacterData.h"

#include <set>
#include <string>
#include <vector>


namespace RevLanguage {
    
    
    class NonHomologousCharacterData : public AbstractCharacterData {
        
    public:
        
        NonHomologousCharacterData(void);                                                                                                 //!< Constructor requires character type
        
        MethodTable                                                 getCharacterDataMethods(void) const;                                                                            //!< Get the common character data methods
        RevPtr<RevVariable>                                         executeCharacterDataMethod(const std::string& name, const std::vector<Argument>& args, bool &found, RevBayesCore::AbstractCharacterData *charDataObject);            //!< Override to map member methods to internal functions
    
//        virtual RevBayesCore::NonHomologousCharacterData&           getValue(void) = 0;
//        virtual const RevBayesCore::NonHomologousCharacterData&     getValue(void) const  = 0;

        
    };
    
}


#endif

