#include "TypedDagNode.h"

#include "NexusWriter.h"
#include "RbSettings.h"
#include <cstdint>
#include "RbUtil.h"
#include "Simplex.h"
#include "StringUtilities.h"
#include "TraceNumeric.h"
#include "TraceTree.h"

#include <ostream>
#include <string>
#include <limits>

///////////////////////
// createTraceObject //
///////////////////////

using namespace RevBayesCore;

template<>
AbstractTrace*  TypedDagNode<std::int64_t>::createTraceObject(void) const
{
    return new TraceNumericInteger();
}

template<>
AbstractTrace*  TypedDagNode<double>::createTraceObject(void) const
{
    return new TraceNumeric();
}

template<>
AbstractTrace*  TypedDagNode<RbVector<double> >::createTraceObject(void) const
{
    return new TraceNumericVector();
}
    
template<>
AbstractTrace*  TypedDagNode<Simplex>::createTraceObject(void) const
{
    return new TraceSimplex();
}
    
template<>
AbstractTrace*  TypedDagNode<Tree>::createTraceObject(void) const
{
    return new TraceTree( getValue().isRooted() );
}

    
/////////////////////
// isSimpleNumeric //
/////////////////////
template<>
bool TypedDagNode<std::int64_t>::isSimpleNumeric(void) const
{
    return true;
} 
    
template<>
bool TypedDagNode<double>::isSimpleNumeric(void) const
{
    return true;
}

template<>
bool TypedDagNode<RbVector<std::int64_t> >::isSimpleNumeric(void) const
{
    return true;
}
    
template<>
bool TypedDagNode<RbVector<double> >::isSimpleNumeric(void) const
{
    return true;
}
    
template<>
bool TypedDagNode<Simplex>::isSimpleNumeric(void) const
{
    return true;
}
    
////////////////
// printValue //
////////////////
template<>
void TypedDagNode<double>::printValue(std::ostream &o, const std::string & /*sep*/, int l, bool left, bool /*user*/, bool simple, bool flatten) const
{
    std::stringstream ss;

    // if simple == FALSE, print with maximum precision allowed
    if (!simple)
    {
        ss.precision(std::numeric_limits<double>::digits10);
    }

    // otherwise, use standard RB precision
    else
    {
        ss.precision(RbSettings::userSettings().getOutputPrecision());
    }
    ss << getValue();
    std::string s = ss.str();
    if ( l > 0 )
    {
        StringUtilities::fillWithSpaces(s, l, left);
    }
    o << s;
}

    
template<>
void TypedDagNode<std::int64_t>::printValue(std::ostream &o, const std::string & /*sep*/, int l, bool left, bool /*user*/, bool /*simple*/, bool /*flatten*/) const
{
        
    std::stringstream ss;
    ss << getValue();
    std::string s = ss.str();
    if ( l > 0 )
    {
        StringUtilities::fillWithSpaces(s, l, left);
    }
    o << s;
}
    
    
template<>
void TypedDagNode<unsigned int>::printValue(std::ostream &o, const std::string & /*sep*/, int l, bool left, bool /*user*/, bool /*simple*/, bool /*flatten*/) const
{
        
    std::stringstream ss;
    ss << getValue();
    std::string s = ss.str();
    if ( l > 0 )
    {
        StringUtilities::fillWithSpaces(s, l, left);
    }
    o << s;
}
    
    
template<>
void TypedDagNode<std::string>::printValue(std::ostream &o, const std::string & /*sep*/, int l, bool left, bool /*user*/, bool /*simple*/, bool /*flatten*/) const
{
        
    std::stringstream ss;
//        ss << "\"" << getValue() << "\"";
    ss << getValue();
    std::string s = ss.str();
    if ( l > 0 )
    {
        StringUtilities::fillWithSpaces(s, l, left);
    }
    o << s;
}
    
    
    
