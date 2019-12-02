#ifndef RlHomologousCharacterData_H
#define RlHomologousCharacterData_H

#include "HomologousCharacterData.h"
#include "RlAbstractCharacterData.h"

#include <set>
#include <string>
#include <vector>


namespace RevLanguage {
    
    
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

