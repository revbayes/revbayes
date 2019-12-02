#include "Simplex.h"

#include <iomanip>
#include <string>

#include "TypedDagNode.h"
#include "IsDerivedFrom.h"
#include "RbException.h"
#include "RbVectorImpl.h"
#include "Serializer.h"
#include "StringUtilities.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { class Serializable; }

using namespace RevBayesCore;


Simplex::Simplex( void ) : RbVector<double>()
{
    
}


Simplex::Simplex( size_t n ) : RbVector<double>( n )
{
    normalize();
}


Simplex::Simplex( size_t n, double v) : RbVector<double>(n,v)
{
    normalize();
}


Simplex::Simplex( const std::vector<double> &s ) : RbVector<double>(s)
{
    normalize();
}


Simplex::Simplex( const Simplex &m ) : RbVector<double>( m )
{
    
}


Simplex::~Simplex( void )
{
}


RbVector<double> Simplex::operator+(double a) const
{
    size_t n = size();
    RbVector<double> result(n, 0.0);
    
    for (size_t i = 0; i < n; ++i)
    {
        result[i] = this->operator[](i) + a;
    }
    
    return result;
}

RbVector<double> Simplex::operator-(double a) const
{
    size_t n = size();
    RbVector<double> result(n, 0.0);
    
    for (size_t i = 0; i < n; ++i)
    {
        result[i] = this->operator[](i) - a;
    }
    
    return result;
}

RbVector<double> Simplex::operator*(double a) const
{
    size_t n = size();
    RbVector<double> result(n, 0.0);
    
    for (size_t i = 0; i < n; ++i)
    {
        result[i] = this->operator[](i) * a;
    }
    
    return result;
}

RbVector<double> RevBayesCore::operator*(const double &a, const Simplex& b)
{
    size_t n = b.size();
    RbVector<double> result(n, 0.0);
    
    for (size_t i = 0; i < n; ++i)
    {
        result[i] = b[i] * a;
    }
    
    return result;
}


Simplex* Simplex::clone( void ) const
{
    return new Simplex(*this);
}


void Simplex::executeMethod(const std::string &n, const std::vector<const DagNode *> &args, double &rv) const
{
    
    if ( n == "[]" )
    {
        long index = static_cast<const TypedDagNode<long> *>( args[0] )->getValue()-1;
        rv = this->operator[](index);
    }
    else
    {
        throw RbException("A simplex object does not have a member method called '" + n + "'.");
    }
    
}

void Simplex::initFromString( const std::string &s )
{
    this->clear();
    std::string sub = s.substr( 1, s.size()-2);
    std::vector<std::string> elements;
    StringUtilities::stringSplit(sub,",", elements);
    for (size_t i=0; i<elements.size(); ++i)
    {
        double value;
        RevBayesCore::Serializer<double, IsDerivedFrom<double, Serializable>::Is >::ressurectFromString( &value, elements[i] );
        this->push_back( value );
    }
    normalize();
    
}


void Simplex::normalize( void )
{
    double sum = 0.0;
    for (size_t i=0; i<size(); ++i)
    {
        double v = this->operator[](i);
        if ( v < 0.0)
        {
            throw RbException("Cannot construct a simplex from negative numbers!");
        }
        sum += v;
    }
    for (size_t i=0; i<size(); ++i)
    {
         this->operator[](i) /= sum;
    }
}


std::ostream& RevBayesCore::operator<<(std::ostream& o, const Simplex& x)
{
    
    std::streamsize previousPrecision = o.precision();
    std::ios_base::fmtflags previousFlags = o.flags();
    
    o << "[ ";
    o << std::fixed;
    o << std::setprecision(4);
    
    // print the RbMatrix with each column of equal width and each column centered on the decimal
    for (size_t i=0; i < x.size(); i++)
    {
        if (i > 0)
            o << ",  ";
        
        o << x[i];
        
        if (i == x.size()-1)
            o << " ]";
        
    }
    
    o.setf(previousFlags);
    o.precision(previousPrecision);
    
    return o;

}
