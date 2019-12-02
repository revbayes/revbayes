#ifndef RlStoppingRule_H
#define RlStoppingRule_H

#include "StoppingRule.h"
#include "TypedDagNode.h"
#include "WorkspaceToCoreWrapperObject.h"

#include <ostream>
#include <string>

namespace RevLanguage {
    
    class StoppingRule : public WorkspaceToCoreWrapperObject<RevBayesCore::StoppingRule> {
        
    public:
        
        StoppingRule(void);                                                                                                                     //!< Default constructor
        
        // Basic utility functions
        virtual StoppingRule*                       clone(void) const = 0;                                                              //!< Clone object
        virtual void                                constructInternalObject(void) = 0;                                                  //!< We construct the a new internal object.
        static const std::string&                   getClassType(void);                                                                 //!< Get Rev type
        static const TypeSpec&                      getClassTypeSpec(void);                                                             //!< Get class type spec
        
    protected:
        
        
    };
    
}

#endif
