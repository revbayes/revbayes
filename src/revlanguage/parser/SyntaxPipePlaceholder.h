/**
 * @file
 * This file contains the declaration of SyntaxUnaryExpr, which is
 * used to hold unary expressions in the syntax tree.
 *
 * @brief Declaration of SyntaxUnaryExpr
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 *
 * $Id$
 */

#ifndef SyntaxPipePlaceholder_H
#define SyntaxPipePlaceholder_H

#include "SyntaxElement.h"

#include <iostream>
#include <vector>


/**
 * This is the class used to hold binary expressions in the syntax tree.
 *
 * We store the operands and a flag signalling the type of operation to
 * be performed when getValue is called or to be represented when
 * getDAGNodeExpr is called.
 *
 */

namespace RevLanguage {

    /**
     * @brief Pipe placeholder
     *
     * A pipe placeholder can occur only at the top level function call
     * of a pipe expression.
     */
    class SyntaxPipePlaceholder : public SyntaxElement {

    public:
        // Basic utility functions
        SyntaxPipePlaceholder*      clone() const;                                                  //!< Clone object
        
        // Regular functions
        RevPtr<RevVariable>         evaluateContent(const std::shared_ptr<Environment>& env, bool dynamic=false);          //!< Get semantic value
    };

}

#endif

