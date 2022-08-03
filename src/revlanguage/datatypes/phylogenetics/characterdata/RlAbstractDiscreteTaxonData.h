#ifndef RlAbstractDiscreteTaxonData_H
#define RlAbstractDiscreteTaxonData_H

#include <string>
#include <vector>
#include <iosfwd>

#include "AbstractDiscreteTaxonData.h"
#include "ModelObject.h"
#include "ConstantNode.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "RevPtr.h"
#include "RlConstantNode.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"


namespace RevLanguage {
class Argument;
class RevVariable;
class TypeSpec;
    
    
    class AbstractDiscreteTaxonData : public ModelObject<RevBayesCore::AbstractDiscreteTaxonData> {

    public:
        
                                                AbstractDiscreteTaxonData(void);                                                        //!< Constructor requires character type
                                                AbstractDiscreteTaxonData(const RevBayesCore::AbstractDiscreteTaxonData &d);            //!< Constructor requires character type
                                                AbstractDiscreteTaxonData(RevBayesCore::AbstractDiscreteTaxonData *d);                  //!< Constructor requires character type
                                                AbstractDiscreteTaxonData(RevBayesCore::TypedDagNode<RevBayesCore::AbstractDiscreteTaxonData>*d);   //!< Constructor requires character type
        
        
//        typedef RevBayesCore::AbstractDiscreteTaxonData valueType;
        
        // Basic utility functions
        AbstractDiscreteTaxonData*              clone(void) const;                                                                      //!< Clone object
        static const std::string&               getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                  getClassTypeSpec(void);                                                                 //!< Get class type spec
        const TypeSpec&                         getTypeSpec(void) const;                                                                //!< Get language type of the object
        
        // Member method inits
        virtual RevPtr<RevVariable>             executeMethod(const std::string& name, const std::vector<Argument>& args, bool &found); //!< Override to map member methods to internal functions
        
        std::string                             getGuiName(void) { return "Taxon Data"; }
        std::string                             getGuiUnicodeSymbol(void) { return "S"; }
        std::string                             getGuiInfo(void) { return ""; }
    };
    
}


#endif

