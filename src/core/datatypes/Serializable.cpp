#include "Serializable.h"
#include "RbFileManager.h"

using namespace RevBayesCore;

// Serialize (resurrect) the object from a file
void Serializable::initFromFile( const std::string &dir, const std::string &fn )
{
    RbFileManager fm = RbFileManager(dir, fn + ".out");
    fm.createDirectoryForFile();
    
    // open the stream to the file
    std::fstream inStream;
    inStream.open( fm.getFullFileName().c_str(), std::fstream::in);
    
    
    std::string s = "";
    while ( inStream.good() )
    {
        
        // Read a line
        std::string line;
        fm.safeGetline( inStream, line );
        
        // append
        s += line;
        
    }
    
    return initFromString( s );
}

// Write this object into a file in its default format.
void Serializable::writeToFile( const std::string &dir, const std::string &fn ) const
{
    
    RbFileManager fm = RbFileManager(dir, fn + ".out");
    fm.createDirectoryForFile();
    
    // open the stream to the file
    std::fstream outStream;
    outStream.open( fm.getFullFileName().c_str(), std::fstream::out);
    
    // write the value of the node
    outStream << this;
    outStream << std::endl;
    
    // close the stream
    outStream.close();
    
}
