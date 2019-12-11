#ifndef RlValidationAnalysis_H
#define RlValidationAnalysis_H

#include "ValidationAnalysis.h"
#include "TypedDagNode.h"
#include "WorkspaceToCoreWrapperObject.h"

#include <ostream>
#include <string>

namespace RevLanguage {
    
    
    /**
     * @brief RevLanguage wrapper class for the validation analysis object.
     *
     * @copydetails RevBayesCore::ValidationAnalysis
     * @see RevBayesCore::ValidationAnalysis for the internal object
     */
    class ValidationAnalysis : public WorkspaceToCoreWrapperObject<RevBayesCore::ValidationAnalysis> {
        
    public:
        
        ValidationAnalysis(void);
        
        // Basic utility functions
        virtual ValidationAnalysis*                 clone(void) const;                                                                      //!< Deep copy of the object
        void                                        constructInternalObject(void);                                                          //!< Construct a new internal ValidationAnalysis object.
        std::string                                 getConstructorFunctionName(void) const;                                                 //!< Get the name used for the constructor function in Rev.
        static const std::string&                   getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                      getClassTypeSpec(void);                                                                 //!< Get class type spec
        const MemberRules&                          getParameterRules(void) const;                                                          //!< Get member rules (const)
        virtual const TypeSpec&                     getTypeSpec(void) const;                                                                //!< Get language type of the object
        
        // Member method inits
        virtual RevPtr<RevVariable>                 executeMethod(const std::string& name, const std::vector<Argument>& args, bool &f);     //!< Map member methods to internal functions
        
    protected:
        
        void                                        initializeMethods(void);                                                                //!< Initialize the member methods
        virtual void                                printValue(std::ostream& o) const;                                                      //!< Print value (for user)
        void                                        setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        RevPtr<const RevVariable>                   simulations; //!< number of analyses to run
        RevPtr<const RevVariable>                   sampler;  //!< sampler object to run
        
    };
    
}

#endif
