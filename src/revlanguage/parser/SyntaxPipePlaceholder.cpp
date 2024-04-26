#include <cstddef>
#include <iostream>
#include <vector>
#include <set>
#include <string>

#include "Environment.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RbException.h"
#include "SyntaxPipePlaceholder.h"

using namespace RevLanguage;


/** Type-safe clone of syntax element */
SyntaxPipePlaceholder* SyntaxPipePlaceholder::clone () const
{
    return new SyntaxPipePlaceholder(*this);
}


RevPtr<RevVariable> SyntaxPipePlaceholder::evaluateContent( Environment& /* env */, bool /* dynamic */ )
{
    throw RbException()<<"Invalid use of pipe placeholder '_'";
}

