#include "RlCommandLineOutputStream.h"

#include <iostream>

void CommandLineOutputStream::output(const std::string &o) const
{
    
    std::cout << o;
    
}


void CommandLineOutputStream::outputEndOfLine( void ) const
{
    
    std::cout << std::endl;
    
}
