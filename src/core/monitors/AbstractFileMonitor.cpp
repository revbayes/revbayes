#include "AbstractFileMonitor.h"

#include <string>

#include "RbFileManager.h"
#include "Cloneable.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;


AbstractFileMonitor::AbstractFileMonitor(DagNode *n, std::uint64_t g, const path &fname, bool ap, bool wv) : Monitor(g,n),
    out_stream(),
    filename( fname ),
    working_file_name( fname ),
    append(ap),
    flatten( true ),
    write_version( wv )
{}


AbstractFileMonitor::AbstractFileMonitor(const std::vector<DagNode *> &n, std::uint64_t g, const path &fname, bool ap, bool wv) : Monitor(g,n),
    out_stream(),
    filename( fname ),
    working_file_name( fname ),
    append(ap),
    flatten( true ),
    write_version( wv )
{}


AbstractFileMonitor::AbstractFileMonitor(const AbstractFileMonitor &f) : Monitor( f ),
    out_stream()
{    
    filename            = f.filename;
    working_file_name   = f.working_file_name;
    append              = f.append;
    flatten             = f.flatten;
    write_version       = f.write_version;
    
    if ( f.out_stream.is_open() == true )
    {
        openStream( true );
    }
    
}


AbstractFileMonitor::~AbstractFileMonitor(void)
{
    // we should always close the stream when the object is deleted
    if (out_stream.is_open())
    {
        closeStream();
    }   
}


void AbstractFileMonitor::addFileExtension(const std::string &s, bool dir)
{   
    if (dir)
        working_file_name = filename.parent_path() / s / filename.filename();
    else
        working_file_name = appendToStem(filename, s);
}


void AbstractFileMonitor::closeStream()
{
    out_stream.close();
}


bool AbstractFileMonitor::isFileMonitor( void ) const
{
    return true;
}


void AbstractFileMonitor::openStream( bool reopen )
{
    createDirectoryForFile( working_file_name );
            
    // open the stream to the file
    if ( append == true || reopen == true )
    {
        out_stream.open( working_file_name.string(), std::fstream::in | std::fstream::out | std::fstream::app);
    }
    else
    {
        out_stream.open( working_file_name.string(), std::fstream::out);
        out_stream.close();
        out_stream.open( working_file_name.string(), std::fstream::in | std::fstream::out);
    }
        
}


/**
 * Set whether to append to an existing file.
 *
 * \param[in]   tf   new flag value
 */
void AbstractFileMonitor::setAppend(bool tf)
{  
    append = tf;    
}


/**
 * Set whether to print the software version.
 *
 * \param[in]   tf   new flag value
 */
void AbstractFileMonitor::setPrintVersion(bool tf)
{
    
    write_version = tf;
    
}
