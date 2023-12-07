#include "Serializable.h"

#include <ostream>
#include <string>

using namespace RevBayesCore;

// Serialize (resurrect) the object from a file
void Serializable::initFromFile( const path &dir, const std::string &fn )
{
    path filename = dir / (fn + ".out");
    
    // open the stream to the file
    std::ifstream inStream( filename.string() );
    
    std::string s = "";
    while ( inStream.good() )
    {
        
        // Read a line
        std::string line;
        safeGetline( inStream, line );
        
        // append
        s += line;
        
    }
    
    return initFromString( s );
}

// Write this object into a file in its default format.
void Serializable::writeToFile( const path &dir, const std::string &fn ) const
{
    path filename = dir / (fn + ".out");
    create_directories( dir );
    
    // open the stream to the file
    std::ofstream outStream( filename.string() );
    
    // write the value of the node
    outStream << this;
    outStream << std::endl;
    
    // close the stream
    outStream.close();
    
}
