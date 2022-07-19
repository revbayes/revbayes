#include "RbException.h"

#include <iostream>


/** Default constructor */
RbException::RbException(void) :
    exception_type(DEFAULT),
    message()
{
}


/** Message constructor */
RbException::RbException(const std::string& msg) :
    exception_type(DEFAULT),
    message(msg)
{
}


/** General constructor */
RbException::RbException(ExceptionType type, const std::string& msg) :
    exception_type(type),
    message(msg)
{
}


RbException::ExceptionType RbException::getExceptionType(void) const
{

    return exception_type;
}


std::string RbException::getMessage(void) const
{

    return message;
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
    
    o << error_type << ":\t" << message;
}

void RbException::setMessage(const std::string& msg)
{
    message = msg;
}


void RbException::prepend(const std::string& prefix)
{
    message = prefix + message;
}
