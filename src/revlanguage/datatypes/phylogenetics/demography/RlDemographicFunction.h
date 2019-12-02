#ifndef RlDemographicFunction_H
#define RlDemographicFunction_H

#include "DemographicFunction.h"
#include "RevVariable.h"
#include "TypedDagNode.h"
#include "WorkspaceToCoreWrapperObject.h"

#include <ostream>
#include <cstring>

namespace RevLanguage {
    
    class DemographicFunction : public WorkspaceToCoreWrapperObject<RevBayesCore::DemographicFunction> {
        
    public:
        
        DemographicFunction(void);                                                                                                          //!< Default constructor
        DemographicFunction(RevBayesCore::DemographicFunction *f);                                                                          //!< constructor with internal value
        
        // Basic utility functions
        virtual DemographicFunction*                clone(void) const = 0;                                                                  //!< Clone object
        static const std::string&                   getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                      getClassTypeSpec(void);                                                                 //!< Get class type spec
        virtual void                                printValue(std::ostream& o, bool user) const;                                           //!< Print value (for user)
        virtual void                                printValue(std::ostream& o) const;                                                      //!< Print value (for user)


        // Member method functions
        virtual RevPtr<RevVariable>                 executeMethod(const std::string& name, const std::vector<Argument>& args, bool &f);         //!< Map member methods to internal functions
        
        std::vector<std::string>                    getConstructorFunctionAliases(void) const;                                                  //!< Get the alternative names used for the constructor function in Rev.
        std::string                                 getConstructorFunctionName(void) const;
        virtual std::vector<std::string>            getDemographicFunctionAliases(void) const { return std::vector<std::string>(); }            //!< Get the alternative names used for the constructor function in Rev.
        virtual std::string                         getDemographicFunctionName(void) const = 0;
        
    protected:
        
    };
    
}

#endif
