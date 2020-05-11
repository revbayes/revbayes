/**
 * @file
 * This file contains the declaration of RevLanguage wrapper of the AncestralStateTrace class.
 *
 * @brief Declaration of RlAncestralStateTrace
 *
 * (c) Copyright 2014-
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 *
 * $Id: RlAncestralStateTrace.h $
 */


#ifndef RlAncestralStateTrace_H
#define RlAncestralStateTrace_H

#include <ostream>
#include <vector>

#include "Trace.h"
#include "WorkspaceToCoreWrapperObject.h"
#include "RevObject.h"
#include "RevPtr.h"

namespace RevLanguage {
class Argument;
class RevVariable;
class TypeSpec;
    
    
    class AncestralStateTrace : public WorkspaceToCoreWrapperObject<RevBayesCore::AncestralStateTrace> {
        
    public:
        
        AncestralStateTrace(void);                                                                                                          //!< Default constructor
        AncestralStateTrace(const RevBayesCore::AncestralStateTrace& x);													            	//!< Copy constructor                                                                                           //!< Copy constructor
        AncestralStateTrace(const AncestralStateTrace& x);
		
        // Basic utility functions
        virtual AncestralStateTrace*                clone(void) const;                                                          //!< Clone object
        void                                        constructInternalObject(void);                                                          //!< We construct the a new internal MCMC object.
        static const std::string&                   getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                      getClassTypeSpec(void);                                                                 //!< Get class type spec
        const MemberRules&                          getParameterRules(void) const;                                                             //!< Get member rules (const)
        virtual const TypeSpec&                     getTypeSpec(void) const;                                                                //!< Get language type of the object
        
        // Member method inits
        virtual RevPtr<RevVariable>                 executeMethod(const std::string& name, const std::vector<Argument>& args, bool &f);  //!< Override to map member methods to internal functions
        
    protected:
        
        virtual void                                printValue(std::ostream& o) const;                                                      //!< Print value (for user)
        void                                        setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);     //!< Set member variable
		
    };
    
	// Global functions using the class
    std::ostream&                       operator<<(std::ostream& o, const AncestralStateTrace& x);                                            //!< Overloaded output operator

}

#endif
