#ifndef Move_SubtreeSwapExtendingNonclock_H
#define Move_SubtreeSwapExtendingNonclock_H

#include "RlMove.h"
#include "TypedDagNode.h"

#include <ostream>
#include <string>

namespace RevLanguage {
    
    class Move_SubtreeSwapExtendingNonclock : public Move {
        
    public:
        
        Move_SubtreeSwapExtendingNonclock(void);                                                                                                                 //!< Default constructor
        
        // Basic utility functions
        virtual Move_SubtreeSwapExtendingNonclock*                       clone(void) const;                                                                      //!< Clone object
        void                                            constructInternalObject(void);                                                          //!< We construct the a new internal move.
        static const std::string&                       getClassName(void);                                                                     //!< Get class name
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                     getMoveName(void) const;                                                                //!< Get the name used for the constructor function in Rev.
        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)
        virtual const TypeSpec&                         getTypeSpec(void) const;                                                                //!< Get language type of the object
        virtual void                                    printValue(std::ostream& o) const;                                                      //!< Print value (for user)
        
    protected:
        
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        RevPtr<const RevVariable>                       tree;
        RevPtr<const RevVariable>                       extension_prob;
    };
    
}

#endif
