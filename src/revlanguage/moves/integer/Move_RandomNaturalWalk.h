#ifndef Move_RandomNaturalWalk_H
#define Move_RandomNaturalWalk_H

#include "RlMove.h"
#include "TypedDagNode.h"

#include <ostream>
#include <string>

namespace RevLanguage {
    
    
    /**
     * The RevLanguage wrapper of the random-integer-walk move.
     *
     * The RevLanguage wrapper of the random-integer-walk move simply
     * manages the interactions through the Rev with our core.
     * That is, the internal move object can be constructed and hooked up
     * in a DAG-nove (variable) that it works on.
     * See the RandomNaturalWalkMove.h for more details.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-01-28, version 1.0
     *
     */
    class Move_RandomNaturalWalk : public Move {
        
    public:
        
        Move_RandomNaturalWalk(void);                                                                                                           //!< Default constructor
        
        // Basic utility functions
        virtual Move_RandomNaturalWalk*             clone(void) const;                                                                          //!< Clone object
        void                                        constructInternalObject(void);                                                              //!< We construct the a new internal SlidingMove.
        static const std::string&                   getClassType(void);                                                                         //!< Get Rev type
        static const TypeSpec&                      getClassTypeSpec(void);                                                                     //!< Get class type spec
        std::vector<std::string>                    getMoveAliases( void ) const;                                                               //!< Get aliases for move name
        std::string                                 getMoveName(void) const;                                                                    //!< Get the name used for the constructor function in Rev.
        const MemberRules&                          getParameterRules(void) const;                                                              //!< Get member rules (const)
        virtual const TypeSpec&                     getTypeSpec(void) const;                                                                    //!< Get language type of the object
        virtual void                                printValue(std::ostream& o) const;                                                          //!< Print value (for user)
        
    protected:
        
        void                                        setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);           //!< Set member variable
        
        RevPtr<const RevVariable>                   x;                                                                                          //!< The variable on which the move works
        
    };
    
}

#endif
