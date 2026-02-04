#ifndef RlCharacterHistory_H
#define RlCharacterHistory_H

#include <string>
#include <vector>
#include <iosfwd>

#include "CharacterHistory.h"
#include "CharacterHistoryDiscrete.h"
#include "ModelObject.h"


namespace RevLanguage {

    class Argument;
    class RevVariable;
    class TypeSpec;
    
    
    class CharacterHistory : public ModelObject<RevBayesCore::CharacterHistoryDiscrete> {
        
    public:
        CharacterHistory(void);                                                                                                                         //!< Constructor requires character type
        CharacterHistory(const RevBayesCore::CharacterHistoryDiscrete &d);                                                                              //!< Constructor requires character type
        CharacterHistory(RevBayesCore::CharacterHistoryDiscrete *d);                                                                                    //!< Constructor requires character type
        CharacterHistory(RevBayesCore::TypedDagNode<RevBayesCore::CharacterHistoryDiscrete>*d);                                                         //!< Constructor requires character type
        
        virtual ~CharacterHistory();
        
        
        // Basic utility functions
        virtual CharacterHistory*                               clone(void) const;                                                                                             //!< Clone object
        static const std::string&                               getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                                  getClassTypeSpec(void);                                                                 //!< Get class type spec
        virtual const TypeSpec&                                 getTypeSpec(void) const;                                                                //!< Get language type of the object
        
        // Member method functions
        virtual RevPtr<RevVariable>                             executeMethod(const std::string& name, const std::vector<Argument>& args, bool &f);     //!< Override to map member methods to internal functions
                
        
    private:
        
        void                                                    initMethods(void);
        
    };
    
}


#endif


