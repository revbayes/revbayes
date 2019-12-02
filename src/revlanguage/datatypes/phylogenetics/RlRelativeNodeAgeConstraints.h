#ifndef RlRelativeNodeAgeConstraints_H
#define RlRelativeNodeAgeConstraints_H

#include "ModelObject.h"
#include "RelativeNodeAgeConstraints.h"
#include "TypedDagNode.h"

#include <ostream>
#include <string>

namespace RevLanguage {
    
    class Real;
    
    class RlRelativeNodeAgeConstraints : public ModelObject<RevBayesCore::RelativeNodeAgeConstraints>  {
        
        public:
            
            RlRelativeNodeAgeConstraints(void);                                                                                                              //!< Default constructor
            RlRelativeNodeAgeConstraints(const RevBayesCore::RelativeNodeAgeConstraints &m);                                                                                        //!< Default constructor
            RlRelativeNodeAgeConstraints(RevBayesCore::RelativeNodeAgeConstraints *m);                                                                                        //!< Default constructor
            RlRelativeNodeAgeConstraints(RevBayesCore::TypedDagNode<RevBayesCore::RelativeNodeAgeConstraints> *d);                                                            //!< Default constructor
            
            
            // Basic utility functions
            RlRelativeNodeAgeConstraints*                   clone(void) const;                                                                      //!< Clone object
            static const std::string&           getClassType(void);                                                                     //!< Get Rev type
            static const TypeSpec&              getClassTypeSpec(void);                                                                 //!< Get class type spec
            const TypeSpec&                     getTypeSpec(void) const;                                                                //!< Get language type of the object
            
            
            // STL-like vector functions provided here
            virtual size_t                      getNumberOfConstraints(void) const;                                               //!< Get number of elements in container
            
           /* Might be desirable to have a function like that, but then we need to define a class for an individual RelativeNodeAgeConstraint
            const Real*                         getConstraint(size_t idx) const;
            */
            
            // Member method inits
            RevPtr<RevVariable>                 executeMethod(const std::string& name, const std::vector<Argument>& args, bool &f);     //!< Override to map member methods to internal functions

            std::string                         getGuiName(void) { return "Age Constraints"; }
            std::string                         getGuiUnicodeSymbol(void) { return "C"; }
            std::string                         getGuiInfo(void) { return ""; }
            
        private:
        
    };
    
    
    
}
#endif /* defined(RlRelativeNodeAgeConstraints_H) */
