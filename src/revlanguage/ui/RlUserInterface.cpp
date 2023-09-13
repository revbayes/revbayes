#include <cstddef>
#include <iostream>
#include <string>

#include "RbSettings.h"
#include "RbUtil.h"
#include "StringUtilities.h"
#include "RlUserInterface.h"
#include "RlUserInterfaceOutputStream.h"

#if defined (RB_MPI)
#include <mpi.h>
#endif

using namespace RevLanguage;


UserInterface::UserInterface( void ) :
    process_id( 0 ),
    output_stream( NULL )
{
#if defined (RB_MPI)
    MPI_Comm_rank(MPI_COMM_WORLD, &process_id);
#endif
}

UserInterface::UserInterface( const UserInterface &u ) :
    process_id( u.process_id )
{
}

/** Ask user a yes/no question */
bool UserInterface::ask(std::string msg)
{

    std::string answer, dummy;
    std::cout << RevBayesCore::RbUtils::PAD << (msg + "? (yes/no) ");     // not using RBOUT or output because we do not want a newline
    std::cin >> answer;
    for (size_t i=0; i<answer.size(); i++)
    {
        answer[i] = char( tolower(answer[i])) ;
    }
    
    while (answer!="y" && answer!="yes" && answer!="n" && answer!="no") 
    {
        std::getline(std::cin, dummy);
        std::cout << std::endl;
		RBOUT("Please answer yes or no.");
        std::cout << RevBayesCore::RbUtils::PAD << (msg + "? (yes/no) "); // see above for choice of std::cout

        std::cin >> answer;
        for (size_t i=0; i<answer.size(); i++)
        {
            answer[i] = char( tolower(answer[i]) );
        }
        
    }
    std::getline(std::cin, dummy);

    if (answer[0] == 'y')
    {
        return true;
    }
    else
    {
        return false;
    }
    
}


/** Print a message and a newline */
void UserInterface::output(std::string msg)
{

    if ( process_id == 0 )
    {
        std::string tmp = StringUtilities::formatStringForScreen( msg, RevBayesCore::RbUtils::PAD, RevBayesCore::RbUtils::PAD, RbSettings::userSettings().getLineWidth() );
        output_stream->output(tmp);
    }
    
}


/** Print a message and a newline without the padding */
void UserInterface::output(std::string msg, const bool hasPadding)
{

    if ( process_id == 0 )
    {
        if (hasPadding == true)
        {
            output(msg);
        }
        else
        {
            output_stream->output(msg);
            output_stream->outputEndOfLine();
        }
        
    }
}


/** Convert to string and then call output to print message string */
void UserInterface::output(std::ostringstream msg)
{

    output( msg.str() );
}


void UserInterface::setOutputStream(UserInterfaceOutputStream *o)
{
    output_stream = o;
}
