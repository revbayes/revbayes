#ifndef RlContinuousTaxonData_H
#define RlContinuousTaxonData_H

#include "ContinuousTaxonData.h"
#include "ModelObject.h"

#include <set>
#include <string>
#include <vector>


namespace RevLanguage {
    
    
    class ContinuousTaxonData : public ModelObject<RevBayesCore::ContinuousTaxonData > {
        
        public:
                                                ContinuousTaxonData(void);                                                                  //!< Constructor requires character type
                                                ContinuousTaxonData(RevBayesCore::ContinuousTaxonData *v);                                  //!< Constructor requires character type
            
            typedef RevBayesCore::ContinuousTaxonData valueType;
            
            // Basic utility functions
            ContinuousTaxonData*                clone(void) const;                                                                          //!< Clone object
            static const std::string&           getClassType(void);                                                                         //!< Get Rev type
            static const TypeSpec&              getClassTypeSpec(void);                                                                     //!< Get class type spec
            const TypeSpec&                     getTypeSpec(void) const;                                                                    //!< Get language type of the object
            
            // Member method inits
            virtual RevPtr<RevVariable>         executeMethod(const std::string& name, const std::vector<Argument>& args, bool &found);     //!< Override to map member methods to internal functions
            
            std::string                         getGuiName(void) { return ""; }
            std::string                         getGuiUnicodeSymbol(void) { return ""; }
            std::string                         getGuiInfo(void) { return ""; }
    };
}

#endif

