#ifndef RlAbstractNonHomologousDiscreteCharacterData_H
#define RlAbstractNonHomologousDiscreteCharacterData_H

#include "AbstractNonHomologousDiscreteCharacterData.h"
#include "ModelObject.h"
#include "RlNonHomologousCharacterData.h"

#include <set>
#include <string>
#include <vector>


namespace RevLanguage {
    
    
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

