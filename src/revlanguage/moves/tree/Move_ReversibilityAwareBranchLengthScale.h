#ifndef Move_ReversibilityAwareBranchLengthScale_H
#define Move_ReversibilityAwareBranchLengthScale_H

#include "RlMove.h"
#include "TypedDagNode.h"

#include <ostream>
#include <string>

namespace RevLanguage {
    
    class Move_ReversibilityAwareBranchLengthScale : public Move {
        
    public:
        
        Move_ReversibilityAwareBranchLengthScale(void);                                                                                                               //!< Default constructor
        
        // Basic utility functions
        virtual Move_ReversibilityAwareBranchLengthScale*       clone(void) const;                                                                      //!< Clone object
        void                                                    constructInternalObject(void);                                                          //!< We construct the a new internal Move.
        static const std::string&                               getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                                  getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                             getMoveName(void) const;                                                                //!< Get the name used for the constructor function in Rev.
        const MemberRules&                                      getParameterRules(void) const;                                                          //!< Get member rules (const)
        virtual const TypeSpec&                                 getTypeSpec(void) const;                                                                //!< Get language type of the object
        virtual void                                            printValue(std::ostream& o) const;                                                      //!< Print value (for user)
                    
    protected:
                    
        void                                                    setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
                    
        /* No q member. This move no longer reads the rate matrix: it scales the two
           branches descending from the root together unconditionally, which is a
           valid proposal whether or not the rate matrix is currently time
           reversible. Holding the rate matrix would have forced RevBayes to touch
           it on every proposal, dirtying the whole CTMC and turning each
           branch-length change into a full-tree likelihood recomputation. */
        RevPtr<const RevVariable>                               tree;
        RevPtr<const RevVariable>                               delta;
        RevPtr<const RevVariable>                               tuning;
        
    };
    
}

#endif

