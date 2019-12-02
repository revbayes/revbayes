#include "RlStringOutputStream.h"


void StringOutputStream::emptyOutputString(void)
{
    current_output = "";
}


const std::string& StringOutputStream::getOutputString( void ) const
{
    
    return current_output;
}

void StringOutputStream::output(const std::string &o) const
{
    
    current_output += o;
    
}


void StringOutputStream::outputEndOfLine( void ) const
{
    
    current_output += "\n";
    
}
