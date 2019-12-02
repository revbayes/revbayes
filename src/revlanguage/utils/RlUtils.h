/*!
 * @file
 * This file contains global constants of the RevLanguage.
 *
 * @brief Names of variable types
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-02-25 10:17:07 +0100 (Sat, 25 Feb 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2010-09-09, version 1.0
 *
 * $Id: RbUtil.h 1247 2012-02-25 09:17:07Z hoehna $
 */

#ifndef RlUtil_H
#define RlUtil_H

#include "RevNullObject.h"
#include "RevObject.h"

#include <string>

namespace RevLanguage {
    
    namespace RlUtils {

        // Empty return type spec
        static const TypeSpec& Void    = TypeSpec( "void", NULL );

//        const std::string Void                                   = "void";

        class RlTypeConverter {
            public:
            static RevObject*                toReal(double x);
            static RevObject*                toString(const std::string &x);
        };
    }

}

#endif
