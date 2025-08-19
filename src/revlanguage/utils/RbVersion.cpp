#include "RbVersion.h"

#include <string>

#include "GitVersion.h"
#include "StringUtilities.h"

#ifdef RB_MPI
#include <mpi.h>
#endif

RbVersion::RbVersion( void )
{
    
}

std::string RbVersion::getDate( void ) const
{
    std::string date = build_date;
    return date;
}

std::string RbVersion::getGitBranch( void ) const
{
    std::string git_branch = build_git_branch;
    return git_branch;
}

std::string RbVersion::getGitCommit( void ) const
{
    std::string git_commit = build_git_sha;
    return git_commit;
}

std::string RbVersion::getVersion( void ) const
{
    return "1.2.6-preview";
}


std::string RbVersion::getHeader( void ) const
{
    
    std::string header = "";
    header += "\n";
    header += "RevBayes version (" + getVersion() + ")\n";
    header += "Build from " + getGitBranch() + " (" + getGitCommit() + ") on " + getDate() + "\n";
#ifdef RB_MPI
    int num_processes = -1;
    MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
    header += "MPI enabled: using " + StringUtilities::to_string(num_processes) + " processes.\n";
#endif
    header += "\n";
    header += "Visit the website www.RevBayes.com for more information about RevBayes.\n";
    header += "\n";
    header += "RevBayes is free software released under the GPL license, version 3. Type 'license()' for details.\n";
    header += "\n";
    header += "To quit RevBayes type 'quit()' or 'q()'.\n";
    
    return header;
    
}



