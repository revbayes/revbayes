#include "AbstractFileMonitor.h"

#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <string>

#include "nlohmann-json.h"

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


path AbstractFileMonitor::getMonitorFileName( void ) const
{
    return working_file_name;
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


void AbstractFileMonitor::truncateAfterGeneration(std::uint64_t lastGen)
{
    std::ifstream infile( working_file_name.string(), std::ios::binary );
    if ( !infile.good() ) return;

    std::string line;

    // detect format from the first data line
    bool is_json = false;
    std::streampos data_start;
    while ( true )
    {
        data_start = infile.tellg();
        if ( !safeGetline(infile, line) ) return;
        if ( line.empty() ) continue;
        if ( std::isdigit(static_cast<unsigned char>(line[0])) ) { is_json = false; break; }
        if ( line[0] == '{' )                                     { is_json = true;  break; }
        // otherwise a header/comment line — keep scanning
    }
    infile.seekg(data_start);

    if ( is_json )
    {
        while ( true )
        {
            std::streampos pos = infile.tellg();
            if ( !safeGetline(infile, line) ) break;
            if ( line.empty() ) continue;

            try
            {
                auto j = nlohmann::json::parse(line);
                if ( j.contains("Iteration") && j["Iteration"].is_number_unsigned() )
                {
                    std::uint64_t gen = j["Iteration"].get<std::uint64_t>();
                    if ( gen > lastGen )
                    {
                        infile.close();
                        std::filesystem::resize_file(working_file_name, pos);
                        return;
                    }
                }
            }
            catch ( const nlohmann::json::exception& ) { }
        }
    }
    else
    {
        while ( true )
        {
            std::streampos pos = infile.tellg();
            if ( !safeGetline(infile, line) ) break;
            if ( line.empty() ) continue;

            char *end_ptr = NULL;
            std::uint64_t gen = static_cast<std::uint64_t>( std::strtoull(line.c_str(), &end_ptr, 10) );
            if ( end_ptr != line.c_str() && gen > lastGen )
            {
                infile.close();
                std::filesystem::resize_file(working_file_name, pos);
                return;
            }
        }
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
