/**
 * @file Serializer.h
 *
 *
 * @brief Declaration of Serializer. This class clones objects by either calling clone, if the
 * object is derived from Cloneable, or calling the copy constructor. Hence, the Serializer can be
 * used in templates when it is unknown if either the copy constructor (basic classes)
 * or the clone function (abstract classes) exists.
 *
 * (c) Copyright 2009-
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 */

#ifndef Serializer_H
#define Serializer_H

#include <cstdlib>
#include <cstdint>
#include <string>

#include "RbException.h"
#include "RbFileManager.h"

namespace RevBayesCore {
    
    template <typename objType, int type>
    // general case: T is not derived from Cloneable
    // calls copy constructor
    class Serializer {
        
    public:
        //!< Create a clone of the given object.
        static void                     ressurectFromString( objType* /*obj*/, const std::string &s ) { throw RbException("Could not resurrect object from string value:\n" + s); }
        
        static void                     ressurectFromFile( objType* obj, const path &dir, const std::string &fn )
        {
            auto p = dir / (fn + ".txt");

            std::ifstream inStream( p.string() );
            
            std::string s = "";
            while ( inStream.good() )
            {
                
                // Read a line
                std::string line;
                safeGetline( inStream, line );
                
                // append
                s += line;
                
            }
            
            // simply delegate
            Serializer<objType, type>::ressurectFromString( obj, s );
        }

        static void                     writeToFile( const objType &obj, const path &dir, const std::string &fn )
        {
            auto p = dir / (fn + ".txt");

            createDirectoryForFile( p );
            
            // open the stream to the file
            std::ofstream outStream( p.string() );
            
            // write the value of the node
//            Printer<objType, IsDerivedFrom<objType, Printable>::Is >::printForComplexStoring( obj, outStream, "", -1, true );
            outStream << obj;
            outStream << std::endl;
            
            // close the stream
            outStream.close();
        }
        
    };
    
    template<>
    inline void Serializer<double,0>::ressurectFromString( double *obj, const std::string &s ) { *obj = atof( s.c_str()); }
    
    template<>
    inline void Serializer<int,0>::ressurectFromString( int *obj, const std::string &s ) { *obj = atoi( s.c_str()); }

    template<>
    inline void Serializer<std::int64_t,0>::ressurectFromString( std::int64_t *obj, const std::string &s ) { *obj = atoi( s.c_str()); }

    
    template <typename objType>
    // T is derived from Cloneable
    // calls clone
    class Serializer<objType,1> {
        
    public:
        //!< Create a clone of the given object.
        static void                     ressurectFromString( objType *obj, const std::string &s ) { obj->initFromString(s); }
        static void                     ressurectFromFile( objType *obj, const path &dir, const std::string &fn ) { obj->initFromFile(dir,fn); }
        static void                     writeToFile( const objType &obj, const path &dir, const std::string &fn ) { obj.writeToFile(dir,fn);  }
    };
    
    
}



#endif

