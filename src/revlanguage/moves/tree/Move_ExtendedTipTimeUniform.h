#ifndef Move_ExtendedTipTimeUniform_H
#define Move_ExtendedTipTimeUniform_H

#include <ostream>
#include <string>
#include <vector>

#include "RlMove.h"
#include "RevPtr.h"
#include "RevVariable.h"

namespace RevLanguage {
    class TypeSpec;

    /**
     * @brief Rev wrapper class for the ExtendedTipTimeUniform move.
     *
     * Draws a new extinction time for a random extinct tip of an extended tree.
     */
    class Move_ExtendedTipTimeUniform : public Move {

    public:

        Move_ExtendedTipTimeUniform(void);                                                                                          //!< Default constructor

        virtual Move_ExtendedTipTimeUniform*        clone(void) const;                                                                      //!< Clone object
        void                                        constructInternalObject(void);                                                          //!< We construct the a new internal Move.
        static const std::string&                   getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                      getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                 getMoveName(void) const;                                                                //!< Get the name used for the constructor function in Rev.
        const MemberRules&                          getParameterRules(void) const;                                                          //!< Get member rules (const)
        virtual const TypeSpec&                     getTypeSpec(void) const;                                                                //!< Get language type of the object
        virtual void                                printValue(std::ostream& o) const;                                                      //!< Print value (for user)

    protected:

        void                                        setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable

        RevPtr<const RevVariable>                   tree;
    };

}

#endif
