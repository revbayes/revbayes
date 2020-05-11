#ifndef RlAbstractNonHomologousDiscreteCharacterData_H
#define RlAbstractNonHomologousDiscreteCharacterData_H

#include <string>
#include <vector>
#include <iosfwd>

#include "AbstractNonHomologousDiscreteCharacterData.h"
#include "ModelObject.h"
#include "RlNonHomologousCharacterData.h"
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
    
    
    class AbstractNonHomologousDiscreteCharacterData : public ModelObject<RevBayesCore::AbstractNonHomologousDiscreteCharacterData>, public NonHomologousCharacterData {
        
    public:
        AbstractNonHomologousDiscreteCharacterData(void);                                                                                                  //!< Constructor requires character type
        AbstractNonHomologousDiscreteCharacterData(const RevBayesCore::AbstractNonHomologousDiscreteCharacterData &d);                                              //!< Constructor requires character type
        AbstractNonHomologousDiscreteCharacterData(RevBayesCore::AbstractNonHomologousDiscreteCharacterData *d);                                              //!< Constructor requires character type
        AbstractNonHomologousDiscreteCharacterData(RevBayesCore::TypedDagNode<RevBayesCore::AbstractNonHomologousDiscreteCharacterData>*d);                   //!< Constructor requires character type
        
        virtual ~AbstractNonHomologousDiscreteCharacterData();
        
        
        typedef RevBayesCore::AbstractNonHomologousDiscreteCharacterData valueType;

        // Basic utility functions
        virtual AbstractNonHomologousDiscreteCharacterData*     clone(void) const;                                                                      //!< Clone object
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

