/**
 * @file
 * This file contains the declaration of MethodTable, which is
 * used to hold member functions of complex objects.
 *
 * @brief Declaration of MethodTable
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @since Version 1.0, 2009-09-09
 * @extends RbObject
 *
 * $Id$
 */

#ifndef MethodTable_H
#define MethodTable_H

#include "FunctionTable.h"

#include <map>
#include <ostream>
#include <string>
#include <vector>

namespace RevLanguage {

    class Function;

    class MethodTable : public FunctionTable {

    public:
        MethodTable(MethodTable* parent = NULL);                                            //!< Constructor of empty table
        virtual         ~MethodTable() {}                                                   //!< Virtual destructor


        // Basic utility functions
        MethodTable*    clone(void) const { return new MethodTable(*this); }                //!< Clone object

        // MethodTable functions
        void            insertInheritedMethods(const MethodTable& inheritedMethods);        //!< Insert inherited methods
        
    };

}

#endif

