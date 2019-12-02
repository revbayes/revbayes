/* 
 * File:   Move_CorrelationMatrixSingleElementBeta.h
 * Author: Michael R. May
 *
 * Created on 5 August 2017
 */

#ifndef Move_CorrelationMatrixSingleElementBeta_H
#define Move_CorrelationMatrixSingleElementBeta_H

#include <ostream>
#include <vector>

#include "RlMove.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"

namespace RevBayesCore { class RbHelpReference; }

namespace RevLanguage {
class TypeSpec;
    
    class Move_CorrelationMatrixSingleElementBeta : public Move {
        
    public:
        
        Move_CorrelationMatrixSingleElementBeta(void);                                                                                                //!< Default constructor
        
        // Basic utility functions
        virtual Move_CorrelationMatrixSingleElementBeta* clone(void) const;                                                                      //!< Clone object
        void                                             constructInternalObject(void);                                                          //!< We construct the a new internal move.
        static const std::string&                        getClassType(void);                                                                     //!< Get class name
        static const TypeSpec&                           getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                      getMoveName(void) const;                                                                //!< Get the name used for the constructor function in Rev.
        const MemberRules&                               getParameterRules(void) const;                                                          //!< Get member rules (const)
        virtual const TypeSpec&                          getTypeSpec(void) const;                                                                //!< Get language type of the object
        virtual void                                     printValue(std::ostream& o) const;                                                      //!< Print value (for user)
        
    protected:
        
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        RevPtr<const RevVariable>                       v;
        RevPtr<const RevVariable>                       alpha;
        RevPtr<const RevVariable>                       weight;
        RevPtr<const RevVariable>                       tune;
        
    };
    
}


#endif	/* MOVE_CorrelationMatrixSingleElementBeta_H */

