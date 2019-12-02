/**
 * @file
 * This file contains the implementation of RbException, which
 * is used to handle eceptions in RevBayes.
 *
 * @brief Implementation of RbException
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 *
 * $Id$
 */

#include "RbException.h"

#include <string>
#include <iostream>


/** Static string with names of exception types for printing */
std::string RbException::exceptionName[] = { "Default", "Quit", "Missing Variable" };


/** Default constructor */
RbException::RbException(void) :
    exceptionType(DEFAULT),
    message()
{
}


/** Message constructor */
RbException::RbException(const std::string& msg) :
    exceptionType(DEFAULT),
    message(msg)
{
}


/** General constructor */
RbException::RbException(exceptionT type, const std::string& msg) :
    exceptionType(type),
    message(msg)
{
}


std::string RbException::getMessage(void) const
{

    return message;
}


void RbException::print(std::ostream &o) const
{
    
    std::string errorType;
    switch (exceptionType)
    {
        case DEFAULT:
            errorType = "Error";
            break;
        case QUIT:
            errorType = "Quit";
            break;
        case MISSING_VARIABLE:
            errorType = "Missing Variable";
            break;
            
        default:
            errorType = "Error";
    }
    
    o << errorType << ":\t" << message;
}

void RbException::setMessage(std::string msg)
{

    message = msg;
}


