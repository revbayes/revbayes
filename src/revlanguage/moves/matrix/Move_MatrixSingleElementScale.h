/* 
 * File:   Move_MatrixSingleElementScale.h
 * Author: nl
 *
 * Created on 13 juillet 2014, 18:13
 */

#ifndef Move_MatrixSingleElementScale_H
#define Move_MatrixSingleElementScale_H

#include "RlMove.h"
#include "TypedDagNode.h"

#include <ostream>
#include <string>

namespace RevLanguage {
    
    class Move_MatrixSingleElementScale : public Move {
        
    public:
        
        Move_MatrixSingleElementScale(void);                                                                                                //!< Default constructor
        
        // Basic utility functions
        virtual Move_MatrixSingleElementScale*      clone(void) const;                                                                      //!< Clone object
        void                                        constructInternalObject(void);                                                          //!< We construct the a new internal move.
        static const std::string&                   getClassType(void);                                                                     //!< Get class name
        static const TypeSpec&                      getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                 getMoveName(void) const;                                                                //!< Get the name used for the constructor function in Rev.
        const MemberRules&                          getParameterRules(void) const;                                                          //!< Get member rules (const)
        virtual const TypeSpec&                     getTypeSpec(void) const;                                                                //!< Get language type of the object
        virtual void                                printValue(std::ostream& o) const;                                                      //!< Print value (for user)
        
    protected:
        
        void                                        setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        RevPtr<const RevVariable>                   v;
        RevPtr<const RevVariable>                   lambda;
        RevPtr<const RevVariable>                   weight;
        RevPtr<const RevVariable>                   tune;
        
    };
    
}


#endif	/* Move_MatrixSingleElementScale_H */

