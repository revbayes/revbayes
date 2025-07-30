#ifndef RlMove_H
#define RlMove_H

#include "Move.h"
#include "RbHelpMove.h"
#include "TypedDagNode.h"
#include "WorkspaceToCoreWrapperObject.h"

#include <ostream>
#include <string>

namespace RevLanguage {
    
    class Move : public WorkspaceToCoreWrapperObject<RevBayesCore::Move> {
        
    public:
        
                                                    Move(void);                                                                         //!< Default constructor (0.0)
        
        // Basic utility functions
        virtual Move*                               clone(void) const = 0;                                                              //!< Clone object
        virtual void                                constructInternalObject(void) = 0;                                                  //!< We construct the a new internal object.
        static const std::string&                   getClassType(void);                                                                 //!< Get Rev type
        static const TypeSpec&                      getClassTypeSpec(void);                                                             //!< Get class type spec
        virtual const MemberRules&                  getParameterRules(void) const;                                                      //!< Get member rules (const)
        virtual void                                printValue(std::ostream& o, bool user) const;                                       //!< Print value (for user)
        virtual void                                printValue(std::ostream& o) const;                                       //!< Print value (for user)

        std::vector<std::string>                    getConstructorFunctionAliases(void) const;                                              //!< Get the alternative names used for the constructor function in Rev.
        std::string                                 getConstructorFunctionName(void) const;
        virtual std::vector<std::string>            getMoveAliases(void) const { return std::vector<std::string>(); }   //!< Get the alternative names used for the constructor function in Rev.
        virtual std::string                         getMoveName(void) const = 0;

    protected:
        
        virtual void                                setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);   //!< Set member variable
        
        RevPtr<const RevVariable>                   delay;
        RevPtr<const RevVariable>                   weight;
        RevPtr<const RevVariable>                   tune_target;

    };
    
}

#endif
