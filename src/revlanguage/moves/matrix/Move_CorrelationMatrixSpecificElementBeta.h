/* 
 * File:   Move_CorrelationMatrixSpecificElementBeta.h
 * Author: Michael R. May
 *
 * Created on 5 August 2017
 */

#ifndef Move_CorrelationMatrixSpecificElementBeta_H
#define Move_CorrelationMatrixSpecificElementBeta_H

#include "RlMove.h"
#include "TypedDagNode.h"

#include <ostream>
#include <string>

namespace RevLanguage {
    
    class Move_CorrelationMatrixSpecificElementBeta : public Move {
        
    public:
        
        Move_CorrelationMatrixSpecificElementBeta(void);                                                                                                //!< Default constructor
        
        // Basic utility functions
        virtual Move_CorrelationMatrixSpecificElementBeta* clone(void) const;                                                                      //!< Clone object
        void                                             constructInternalObject(void);                                                          //!< We construct the a new internal move.
        static const std::string&                        getClassType(void);                                                                     //!< Get class name
        static const TypeSpec&                           getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                      getMoveName(void) const;                                                                //!< Get the name used for the constructor function in Rev.
        const MemberRules&                               getParameterRules(void) const;                                                          //!< Get member rules (const)
        virtual const TypeSpec&                          getTypeSpec(void) const;                                                                //!< Get language type of the object
        virtual void                                     printValue(std::ostream& o) const;                                                      //!< Print value (for user)
        
    protected:
        
        void                                        setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        RevPtr<const RevVariable>                   v;
        RevPtr<const RevVariable>                   row;
        RevPtr<const RevVariable>                   col;
        RevPtr<const RevVariable>                   alpha;
        RevPtr<const RevVariable>                   weight;
        RevPtr<const RevVariable>                   tune;
        
    };
    
}


#endif	/* MOVE_CorrelationMatrixSpecificElementBeta_H */

