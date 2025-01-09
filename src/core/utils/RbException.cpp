#include "RbException.h"

#include <iostream>


/** Message constructor */
RbException::RbException(const std::string& msg) :
    message(msg)
{
}


/** General constructor */
RbException::RbException(ExceptionType type, const std::string& msg) :
    exception_type(type),
    message(msg)
{
}

/** Copy constructor: this is used when we throw the exception, but should not be used otherwise **/
RbException::RbException(const RbException& E)
    :message(E.message.str())
{
    // Copy formatting flags.
    message.flags(E.message.flags());
    message.precision(E.message.precision());
    message.width(E.message.width());
    // We haven't copies the locate.
}

RbException::ExceptionType RbException::getExceptionType(void) const
{

    return exception_type;
}


std::string RbException::getMessage(void) const
{

    return message.str();
}


void RbException::print(std::ostream &o) const
{
    
    std::string error_type;
    switch (exception_type)
    {
        case DEFAULT:
            error_type = "Error";
            break;
        case BUG:
            error_type = "Bug";
            break;
        case MATH_ERROR:
            error_type = "Mathematical Error";
            break;
        case MISSING_VARIABLE:
            error_type = "Missing Variable";
            break;
        case QUIT:
            error_type = "Quit";
            break;
            
        default:
            error_type = "Error";
    }
    
    o << error_type << ":\t" << message.str();
}

void RbException::setMessage(const std::string& msg)
{
    message.str(msg);
}


void RbException::prepend(const std::string& prefix)
{
    message.str(prefix + message.str());
}
