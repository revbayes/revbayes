/**
 * @file
 * This file contains the declaration of the RevLanguage wrapper of the ModelTrace class.
 *
 * @brief Declaration of RlTreeModelTrace
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date: 2012-08-06 20:14:22 +0200 (Mon, 06 Aug 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2012-09-06, version 1.0
 * @extends RbObject
 *
 * $Id: Real.h 1746 2012-08-06 18:14:22Z hoehna $
 */

#ifndef RlModelTrace_H
#define RlModelTrace_H

#include "Trace.h"
#include "TypedDagNode.h"
#include "WorkspaceToCoreWrapperObject.h"

#include <ostream>
#include <string>

namespace RevLanguage {
    
    class ModelTrace : public WorkspaceToCoreWrapperObject<RevBayesCore::ModelTrace> {
        
    public:
        
        ModelTrace(void);                                                                                                                        //!< Default constructor
        ModelTrace(const RevBayesCore::ModelTrace& x);                                                                                                //!< Copy constructor
        
        // Basic utility functions
        virtual ModelTrace*                         clone(void) const;                                                                      //!< Clone object
        void                                        constructInternalObject(void);                                                          //!< We construct the a new internal MCMC object.
        static const std::string&                   getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                      getClassTypeSpec(void);                                                                 //!< Get class type spec
        const MemberRules&                          getParameterRules(void) const;                                                          //!< Get member rules (const)
        virtual const TypeSpec&                     getTypeSpec(void) const;                                                                //!< Get language type of the object
        
        // Member method inits
        virtual RevPtr<RevVariable>                 executeMethod(const std::string& name, const std::vector<Argument>& args, bool &f);     //!< Override to map member methods to internal functions
        
    protected:
        
        virtual void                                printValue(std::ostream& o) const;                                                      //!< Print value (for user)
        void                                        setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);          //!< Set member variable
        
    };
    
    // Global functions using the class
    std::ostream&                       operator<<(std::ostream& o, const ModelTrace& x);                                //!< Overloaded output operator
    
    
}

#endif
