/**
 * @file Printer.h
 *
 *
 * @brief Declaration of Printer. This class clones objects by either calling clone, if the
 * object is derived from Cloneable, or calling the copy constructor. Hence, the Printer can be
 * used in templates when it is unknown if either the copy constructor (basic classes)
 * or the clone function (abstract classes) exists.
 *
 * (c) Copyright 2009-
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 */

#ifndef Printer_H
#define Printer_H

#include <string>
#include <cstdint>
#include <limits>

#include "StringUtilities.h"

namespace RevBayesCore {
    
    template <typename objType, int>
    // general case: objType is not derived from Printable
    // implements printing for a general object
    class Printer {
        
    public:
        //!< Printer the given object.
        static void                     printForUser( const objType &a, std::ostream &o, const std::string & /*sep*/, int /*l*/, bool /*left*/ ) {
            
            std::int64_t previousPrecision = o.precision();
            std::ios_base::fmtflags previousFlags = o.flags();
            
            std::fixed( o );
            o.precision( 3 );
            
            o << a;
            
            o.setf( previousFlags );
            o.precision( previousPrecision );
        }

        // print with rounding
        static void                     printForSimpleStoring( const objType &a, std::ostream &o, const std::string & /*sep*/, int l, bool left, bool flatten = true)
        {
            std::stringstream ss;
            ss << a;
            std::string s = ss.str();
            if ( l > 0 )
            {
                StringUtilities::fillWithSpaces(s, l, left);
            }
            o << s;
        }

        // print without rounding
        static void                     printForComplexStoring( const objType &a, std::ostream &o, const std::string & /*sep*/, int l, bool left, bool flatten = true )
        {
            std::stringstream ss;
            // set precision of stringstream to max
            ss.precision(std::numeric_limits<double>::digits10);
            ss << a;
            std::string s = ss.str();
            if ( l > 0 )
            {
                StringUtilities::fillWithSpaces(s, l, left);
            }
            o << s;
        }
        
    };
    
    template <typename objType>
    // T is derived from Printable
    // calls printFor... for the corresponding object type
    class Printer<objType,1> {
        
    public:
        //!< Printer the given object.
        static void                     printForUser( const objType &a, std::ostream &o, const std::string &sep, int l, bool left ) { a.printForUser(o, sep, l, left); }
        static void                     printForSimpleStoring( const objType &a, std::ostream &o, const std::string &sep, int l, bool left, bool flatten = true) { a.printForSimpleStoring(o, sep, l, left, flatten); }
        static void                     printForComplexStoring( const objType &a, std::ostream &o, const std::string &sep, int l, bool left, bool flatten = true ) { a.printForComplexStoring(o, sep, l, left, flatten); }

    };

}

#endif

