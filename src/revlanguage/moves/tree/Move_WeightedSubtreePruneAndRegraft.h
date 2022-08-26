/**
 * @file
 * This file contains the declaration of RevLanguage wrapper of WeightedSubtreePruneAndRegraft.
 *
 * @brief Declaration of RlSlidingMove
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date: 2012-08-06 20:14:22 +0200 (Mon, 06 Aug 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2009-11-20, version 1.0
 * @extends RbObject
 *
 * $Id: Real.h 1746 2012-08-06 18:14:22Z hoehna $
 */

#ifndef Move_WeightedSubtreePruneAndRegraft_H
#define Move_WeightedSubtreePruneAndRegraft_H

#include "RlMove.h"
#include "TypedDagNode.h"

#include <ostream>
#include <string>

namespace RevLanguage {
    
    class Move_WeightedSubtreePruneAndRegraft : public Move {
        
    public:
        
        Move_WeightedSubtreePruneAndRegraft(void);                                                                                                                 //!< Default constructor
        
        // Basic utility functions
        virtual Move_WeightedSubtreePruneAndRegraft*    clone(void) const;                                                                      //!< Clone object
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
        RevPtr<const RevVariable>                       alpha;
        RevPtr<const RevVariable>                       num_breaks;
        
    };
    
}

#endif
