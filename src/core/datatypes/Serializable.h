#ifndef Serializable_H
#define Serializable_H

#include <iosfwd>
#include "RbFileManager.h"

namespace RevBayesCore {
    
    class Serializable {
        
    public:
        virtual                         ~Serializable(void) {}
        
        virtual void                    initFromString( const std::string &s ) = 0;                      //!< Serialize (resurrect) the object from a string value
        virtual void                    initFromFile( const path &dir, const std::string &fn );          //!< Read and resurrect this object from a file in its default format.
        virtual void                    writeToFile( const path &dir, const std::string &fn ) const;     //!< Write this object into a file in its default format.
    };
    
}

#endif

