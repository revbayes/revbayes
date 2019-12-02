/* 
 * File:   Move_CorrelationMatrixSingleElementBeta.h
 * Author: Michael R. May
 *
 * Created on 5 August 2017
 */

#ifndef Move_CorrelationMatrixSingleElementBeta_H
#define Move_CorrelationMatrixSingleElementBeta_H

#include "RlMove.h"
#include "TypedDagNode.h"

#include <ostream>
#include <string>

namespace RevLanguage {
    
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
        
        
        std::vector<std::string>                        getHelpAuthor(void) const;                                                              //!< Get the author(s) of this function
        std::string                                     getHelpDescription(void) const;                                                         //!< Get the description for this function
        std::string                                     getHelpDetails(void) const;                                                             //!< Get the more detailed description of the function
        std::string                                     getHelpExample(void) const;                                                             //!< Get an executable and instructive example
        std::vector<RevBayesCore::RbHelpReference>      getHelpReferences(void) const;                                                          //!< Get some references/citations for this function
        std::vector<std::string>                        getHelpSeeAlso(void) const;                                                             //!< Get suggested other functions
        std::string                                     getHelpTitle(void) const;                                                               //!< Get the title of this help entry

        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        RevPtr<const RevVariable>                       v;
        RevPtr<const RevVariable>                       alpha;
        RevPtr<const RevVariable>                       weight;
        RevPtr<const RevVariable>                       tune;
        
    };
    
}


#endif	/* MOVE_CorrelationMatrixSingleElementBeta_H */

