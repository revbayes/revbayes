#ifndef Move_LayeredScaleProposal_H
#define Move_LayeredScaleProposal_H

#include "RlMove.h"
#include "TypedDagNode.h"

#include <ostream>
#include <string>

namespace RevLanguage {


    /**
     * The RevLanguage wrapper of the LayeredScaleProposal move.
     *
     * The RevLanguage wrapper of the LayeredScaleProposal move simply
     * manages the interactions through the Rev with our core.
     * That is, the internal move object can be constructed and hooked up
     * in a DAG-nove (variable) that it works on.
     * See the LayeredScaleProposal.h for more details.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Bastien Boussau)
     * @since 2015-06-23, version 1.0
     *
     */
    class Move_LayeredScaleProposal : public Move {

    public:

        Move_LayeredScaleProposal(void);                                                                                                     //!< Default constructor

        // Basic utility functions
        virtual Move_LayeredScaleProposal*           clone(void) const;                                                                      //!< Clone object
        void                                        constructInternalObject(void);                                                          //!< We construct the a new internal move.
        static const std::string&                   getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                      getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                 getMoveName(void) const;                                                                //!< Get the name used for the constructor function in Rev.
        const MemberRules&                          getParameterRules(void) const;                                                          //!< Get member rules (const)
        virtual const TypeSpec&                     getTypeSpec(void) const;                                                                //!< Get language type of the object
        virtual void                                printValue(std::ostream& o) const;                                                      //!< Print value (for user)

    protected:

        void                                              setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable

        RevPtr<const RevVariable>                         tree;
        RevPtr<const RevVariable>                         lambda;                             //!< The tuning parameter
        RevPtr<const RevVariable>                         tune;                               //!< If autotuning should be used.

    };

}

#endif
