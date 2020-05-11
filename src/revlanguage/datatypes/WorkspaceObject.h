/**
 * @file
 * This file contains the declaration of WorkspaceObject, which is the
 * the Rev abstract base class for all language objects wrapping core
 * datatypes that are used to control the analysis and manipulate model
 * objects, such as monitors, moves, models. All workspace objects have
 * in common that they cannot be included in model DAGs.
 *
 * @brief Declaration of WorkspaceObject
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date: 2012-08-03 17:51:49 +0200 (Fri, 03 Aug 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @since Version 1.0, 2012-08-06
 *
 * $Id: RevObject.h 1734 2012-08-03 15:51:49Z hoehna $
 */

#ifndef WorkspaceObject_H
#define WorkspaceObject_H

#include <iosfwd>

#include "RevObject.h"

namespace RevLanguage {
class TypeSpec;
    
    class WorkspaceObject : public RevObject {
        
    public:
        virtual                            ~WorkspaceObject(void);                                      //!< Virtual destructor
        
        // Basic utility functions you have to override (also getClassTypeSpec()!)
        virtual WorkspaceObject*            clone(void) const = 0;                                      //!< Clone object
        static const std::string&           getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&              getClassTypeSpec(void);                                     //!< Get class type spec
        virtual const TypeSpec&             getTypeSpec(void) const = 0;                                //!< Get the type spec of the instance
        
        virtual bool                        isWorkspaceObject(void) const;                                                            //!< Is this a Workspace object (compared with a workspace object)

    protected:

        WorkspaceObject(void);

    };
}

#endif


