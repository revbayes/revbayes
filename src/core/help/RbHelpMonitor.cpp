#include "RbHelpMonitor.h"

using namespace RevBayesCore;

RbHelpMonitor::RbHelpMonitor( void ) : RbHelpType()
{
    
}



RbHelpMonitor* RbHelpMonitor::clone( void ) const
{
    
    // call the copy constructor
    return new RbHelpMonitor(*this);
}
