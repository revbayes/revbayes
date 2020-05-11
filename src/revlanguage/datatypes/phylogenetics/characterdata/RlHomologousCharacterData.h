#ifndef RlHomologousCharacterData_H
#define RlHomologousCharacterData_H

#include <vector>
#include <iosfwd>

#include "RlAbstractCharacterData.h"
#include "MethodTable.h"
#include "RevPtr.h"

namespace RevBayesCore { class HomologousCharacterData; }


namespace RevLanguage {
class Argument;
class RevVariable;
    
    
    class HomologousCharacterData : public AbstractCharacterData {
        
    public:
        
        HomologousCharacterData(void);                                                                                                                                                                                  //!< Constructor requires character type
        
        MethodTable                                                 getCharacterDataMethods(void) const;                                                                                                                //!< Get the common character data methods
        RevPtr<RevVariable>                                         executeCharacterDataMethod(const std::string& name, const std::vector<Argument>& args, bool &found, RevBayesCore::HomologousCharacterData *obj);    //!< Override to map member methods to internal functions
        
//        virtual RevBayesCore::HomologousCharacterData&              getValue(void) = 0;
//        virtual const RevBayesCore::HomologousCharacterData&        getValue(void) const  = 0;
        
    };
    
}


#endif

