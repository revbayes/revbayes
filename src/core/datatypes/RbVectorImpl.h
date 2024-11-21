#ifndef RbVectorImpl_H
#define RbVectorImpl_H

#include "Cloner.h"
#include "Cloneable.h"
#include "IsDerivedFrom.h"
#include "Printable.h"
#include "Printer.h"
#include "RbConstIterator.h"
#include "RbContainer.h"
#include "RbIterator.h"
#include "Serializable.h"
#include "Serializer.h"
#include "StringUtilities.h"

#include <vector>
#include <algorithm>
#include <iostream>
#include <string>

namespace RevBayesCore {
    
    template <class valueType, int indicator>
    // general case: T is not abstract
    // use actual objects
    class RbVectorImpl : public std::vector<valueType>, public Cloneable, public Serializable, public Printable, public Container {
        
    public:
        
        typedef typename std::vector<valueType>                     vectorType;
        typedef typename std::vector<valueType>::iterator           iterator;
        typedef typename std::vector<valueType>::const_iterator     const_iterator;

        // constructor(s)
        RbVectorImpl() {  }
        RbVectorImpl(size_t n) { for (size_t i = 0; i < n; ++i) this->push_back( valueType() ); }
        RbVectorImpl(size_t n, const valueType &v) { for (size_t i = 0; i < n; ++i) this->push_back( v ); }
        RbVectorImpl(const vectorType &v) { size_t n=v.size(); for (size_t i = 0; i < n; ++i) this->push_back( v[i] ); }
        RbVectorImpl(const RbVectorImpl<valueType,indicator> &v):std::vector<valueType>() { size_t n=v.size(); for (size_t i = 0; i < n; ++i) this->push_back( v[i] ); }
        RbVectorImpl(RbVectorImpl<valueType,indicator> &&v) = default;
        virtual                                            ~RbVectorImpl(void) { }
        
        
        // public member functions
        virtual RbVectorImpl<valueType, indicator>*         clone(void) const = 0;                                                                      //!< Create an independent clone

        // Serialize (resurrect) the object from a string
        virtual void                                        initFromString( const std::string &s )
        {
            this->clear();
            std::string sub = s.substr( 1, s.size()-2);
            std::vector<std::string> elements;
            StringUtilities::stringSplit(sub,",", elements);
            for (size_t i=0; i<elements.size(); ++i)
            {
                valueType value;
                RevBayesCore::Serializer<valueType, IsDerivedFrom<valueType, Serializable>::Is >::ressurectFromString( &value, elements[i] );
                this->push_back( value );
            }

        }
        
        // public (stl-like) vector functions
        RbVectorImpl<valueType,indicator>&                  operator=(const RbVectorImpl<valueType,indicator> &v) = default;
        RbVectorImpl<valueType,indicator>&                  operator=(RbVectorImpl<valueType,indicator> &&v) = default;

        RbIterator<valueType>                               begin(void) { return RbIterator<valueType>( this->std::vector<valueType>::begin() ); }
        RbConstIterator<valueType>                          begin(void) const { return RbConstIterator<valueType>( this->std::vector<valueType>::begin() ); }
        RbIterator<valueType>                               end(void) { return RbIterator<valueType>( this->std::vector<valueType>::end() ); }
        RbConstIterator<valueType>                          end(void) const { return RbConstIterator<valueType>( this->std::vector<valueType>::end() ); }
        RbIterator<valueType>                               erase(RbIterator<valueType> pos) { return this->std::vector<valueType>::erase( pos.getStlIterator() ); }
        RbConstIterator<valueType>                          erase(RbConstIterator<valueType> pos) { return this->std::vector<valueType>::erase( pos ); }

        RbIterator<valueType>                               find(const valueType &x) { return RbIterator<valueType>( std::find(this->std::vector<valueType>::begin(), this->std::vector<valueType>::end(), x) ); }
        RbConstIterator<valueType>                          find(const valueType &x) const { return RbConstIterator<valueType>( std::find(this->std::vector<valueType>::begin(), this->std::vector<valueType>::end(), x) ); }
        virtual size_t                                      size(void) const { return this->std::vector<valueType>::size(); }
        valueType&                                          operator[](size_t i)
        {
            if ( i >= std::vector<valueType>::size() )
            {
                throw RbException() << "Vector index out of range. You tried to access index '" << StringUtilities::to_string(i) << "' for a vector of size '" <<  StringUtilities::to_string(std::vector<valueType>::size()) << "'.";
            }
            return std::vector<valueType>::operator [](i);
        }
        const valueType&                                    operator[](size_t i) const
        {
            if ( i >= std::vector<valueType>::size() )
            {
                throw(RbException("Vector index out of range: "+StringUtilities::to_string(i)+" of "+StringUtilities::to_string(std::vector<valueType>::size())));
            }
            return std::vector<valueType>::operator [](i);
        }
        void                                                setElement(int i, const valueType& v) { (*this)[i] = v;}
        void                                                swap( valueType& a, valueType& b)
        {
            valueType temp = a;
            a = b;
            b = temp;
        }

	json                                                toJSON() const
	{
	    json j;
            for (size_t i=0; i<size(); ++i)
	    {
		if constexpr (std::is_base_of_v<Printable,valueType>)
		    j.push_back(this->operator[](i).toJSON());
		else if constexpr (std::is_convertible_v<valueType,json>)
		    j.push_back(this->operator[](i));
		else
		{
		    std::stringstream ss;
		    ss.precision(std::numeric_limits<double>::digits10);
		    ss << this->operator[](i);
		    j.push_back(ss.str());
		}
	    }
	    return j;
	}

        void                                                printForUser( std::ostream &o, const std::string &sep, int l, bool left ) const
        {
            o << "[";
            for (size_t i=0; i<size(); ++i)
            {
                if (i > 0)
                {
                    o << ",";
                }
                o << " ";
                Printer<valueType, IsDerivedFrom<valueType, Printable>::Is >::printForUser( this->operator[](i), o, sep, l, left );
            }
            o << " ]";
        }
        void                                                printForSimpleStoring( std::ostream &o, const std::string &sep, int l, bool left, bool flatten = true ) const
        {
            if (flatten) {
                for (size_t i = 0; i < size(); ++i) {
                    if (i > 0) {
                        o << sep;
                    }
                    Printer<valueType, IsDerivedFrom<valueType, Printable>::Is>::printForSimpleStoring(
                            this->operator[](i), o, sep, l, left);
                }
            }
            else {
                o << "[";
                for (size_t i=0; i<size(); ++i)
                {
                    if (i > 0)
                    {
                        o << ",";
                    }
                    Printer<valueType, IsDerivedFrom<valueType, Printable>::Is >::printForSimpleStoring( this->operator[](i), o, sep, l, left );
                }
                o << "]";
            }
        }
        void                                                printForComplexStoring( std::ostream &o, const std::string &sep, int l, bool left, bool flatten = true ) const
        {
            o.precision( std::numeric_limits<double>::digits10 );

            if (flatten) {
                for (size_t i=0; i<size(); ++i)
                {
                    if (i > 0)
                    {
                        o << sep;
                    }
                    Printer<valueType, IsDerivedFrom<valueType, Printable>::Is >::printForComplexStoring( this->operator[](i), o, sep, l, left );
                }
            }
            else {
                o << "[";
                for (size_t i=0; i<size(); ++i)
                {
                    if (i > 0)
                    {
                        o << ",";
                    }
                    Printer<valueType, IsDerivedFrom<valueType, Printable>::Is >::printForComplexStoring( this->operator[](i), o, sep, l, left );
                }
                o << "]";
            }

        }
        
    };
    
    template <typename valueType>
    // T is abstract
    // uses pointers
    class RbVectorImpl<valueType,1> : public Cloneable, public Serializable, public Printable, public Container {
        
    public:

        typedef typename std::vector<valueType*>                     vectorType;
        typedef typename std::vector<valueType*>::iterator           iterator;
        typedef typename std::vector<valueType*>::const_iterator     const_iterator;

        // constructor(s)
        RbVectorImpl() {  }
        RbVectorImpl(size_t n) { values = std::vector<valueType*>(n, NULL); }
        RbVectorImpl(size_t n, const valueType &v) { for (size_t i = 0; i < n; ++i) values.push_back( Cloner<valueType, IsDerivedFrom<valueType, Cloneable>::Is >::createClone( v ) ); }
        RbVectorImpl(const vectorType &v) { size_t n=v.size(); for (size_t i = 0; i < n; ++i) values.push_back( Cloner<valueType, IsDerivedFrom<valueType, Cloneable>::Is >::createClone( v[i] ) ); }
        RbVectorImpl(const RbVectorImpl<valueType,1> &v) : values() { size_t n=v.size(); for (size_t i = 0; i < n; ++i) values.push_back( Cloner<valueType, IsDerivedFrom<valueType, Cloneable>::Is >::createClone( v[i] ) ); }
        RbVectorImpl(RbVectorImpl<valueType,1> &&v) = default;
        virtual                                            ~RbVectorImpl(void) { clear(); }
        

        
        // public member functions
        virtual RbVectorImpl<valueType, 1>*                 clone(void) const = 0;                                                                      //!< Create an independent clone
        virtual void                                        initFromString( const std::string &s ) {}                                                 //!< Serialize (resurrect) the object from a string

        // public (stl-like) vector functions
        RbVectorImpl<valueType, 1>&                         operator=(const RbVectorImpl<valueType, 1> &v) {
            if ( this != &v )
            {
                clear();
                size_t n=v.size();
                for (size_t i = 0; i < n; ++i) values.push_back( Cloner<valueType, IsDerivedFrom<valueType, Cloneable>::Is >::createClone( v[i] ) );
            }
            return *this;
        }
        RbVectorImpl<valueType, 1>&                         operator=(RbVectorImpl<valueType, 1> &&v) = default;

        valueType&                                          operator[](size_t i)
        {
            if ( i >= values.size() )
            {
                throw RbException("Vector index out of range") ;
            }
            return *values[i];
        }
        const valueType&                                    operator[](size_t i) const
        {
            if ( i >= values.size() )
            {
                throw RbException("Vector index out of range");
            }
            return *values[i];
        }
        bool                                                operator==(const RbVectorImpl<valueType,1>& x) const { return values == x.values; }                              //!< Equals operator
        bool                                                operator!=(const RbVectorImpl<valueType,1>& x) const { return !operator==(x); }                              //!< Not-Equals operator
        bool                                                operator<(const RbVectorImpl<valueType,1>& x) const { return values < x.values; }
        bool                                                operator<=(const RbVectorImpl<valueType,1>& x) const { return operator<(x) || operator==(x); }

        void                                                clear(void) {
            size_t n = values.size();
            for (size_t i = 0; i < n; ++i)
            {
                valueType* v = values[i];
                delete v;
            }

            values.clear();
        }
        void                                                insert(size_t i, const valueType &v) { delete values[i]; values[i] = Cloner<valueType, IsDerivedFrom<valueType, Cloneable>::Is >::createClone( v ); }
        void                                                push_back(const valueType &v) { values.push_back( Cloner<valueType, IsDerivedFrom<valueType, Cloneable>::Is >::createClone( v ) ); }
        RbIterator<valueType>                               begin(void) { return RbIterator<valueType>( this->values.begin() ); }
        RbConstIterator<valueType>                          begin(void) const { return RbConstIterator<valueType>( this->values.begin() ); }
        RbIterator<valueType>                               end(void) { return RbIterator<valueType>( this->values.end() ); }
        RbConstIterator<valueType>                          end(void) const { return RbConstIterator<valueType>( this->values.end() ); }
        void                                                erase(size_t i) { valueType *tmp=values[i]; values.erase(values.begin()+i); delete tmp; }
        RbIterator<valueType>                               erase(RbIterator<valueType> pos) { valueType *tmp=&(*pos); delete tmp; return RbIterator<valueType>( this->values.erase( pos.getStlIterator() ) );}
        RbConstIterator<valueType>                          erase(RbConstIterator<valueType> pos) { valueType *tmp=*pos; delete tmp; return RbConstIterator<valueType>( this->values.erase( pos.getStlIterator() ) ); }
        RbIterator<valueType>                               find(const valueType &x) { return RbIterator<valueType>( std::find(this->values.begin(), this->values.end(), &x) ); }
        RbConstIterator<valueType>                          find(const valueType &x) const { return RbConstIterator<valueType>( std::find(this->values.begin(), this->values.end(), &x) ); }
        size_t                                              size(void) const { return this->values.size(); }
        void                                                setElement(int i, const valueType& v) { delete values[i]; values[i] = Cloner<valueType, IsDerivedFrom<valueType, Cloneable>::Is >::createClone(v);}

        void                                                swap( valueType& a, valueType& b)
        {
            valueType *temp = Cloner<valueType, IsDerivedFrom<valueType, Cloneable>::Is >::createClone( a );
            a = b;
            b = *temp;

            delete temp;
        }

	json                                                toJSON() const
	{
	    json j;
            for (size_t i=0; i<size(); ++i)
	    {
		if constexpr (std::is_base_of_v<Printable,valueType>)
		    j.push_back(this->operator[](i).toJSON());
		else if constexpr (std::is_convertible_v<valueType,json>)
		    j.push_back(this->operator[](i));
		else
		{
		    std::stringstream ss;
		    ss.precision(std::numeric_limits<double>::digits10);
		    ss << this->operator[](i);
		    j.push_back(ss.str());
		}
	    }
	    return j;
	}

        void                                                printForUser( std::ostream &o, const std::string &sep, int l, bool left ) const
        {
            o << "[";
            for (size_t i=0; i<size(); ++i)
            {
                if (i > 0)
                {
                    o << ",";
                }
                o << " ";
                Printer<valueType, IsDerivedFrom<valueType, Printable>::Is >::printForUser( this->operator[](i), o, sep, l, left );
            }
            o << " ]";
        }
        void                                                printForSimpleStoring( std::ostream &o, const std::string &sep, int l, bool left, bool flatten = true ) const
        {
            // if flatten == TRUE, save each element of vector separately
            if (flatten) {
                for (size_t i = 0; i < size(); ++i) {
                    if (i > 0) {
                        o << sep;
                    }
                    Printer<valueType, IsDerivedFrom<valueType, Printable>::Is>::printForSimpleStoring(
                            this->operator[](i), o, sep, l, left);
                }
            }

            // otherwise, save full vector
            else {
                o << "[";
                for (size_t i=0; i<size(); ++i)
                {
                    if (i > 0)
                    {
                        o << ",";
                    }
                    Printer<valueType, IsDerivedFrom<valueType, Printable>::Is >::printForSimpleStoring( this->operator[](i), o, sep, l, left );
                }
                o << "]";
            }
        }
        void                                                printForComplexStoring( std::ostream &o, const std::string &sep, int l, bool left, bool flatten = true ) const
        {
            // set precision to maximum
            o.precision( std::numeric_limits<double>::digits10 );

            // if flatten == TRUE, save each element of vector separately
            if (flatten) {
                for (size_t i=0; i<size(); ++i)
                {
                    if (i > 0)
                    {
                        o << sep;
                    }
                    Printer<valueType, IsDerivedFrom<valueType, Printable>::Is >::printForComplexStoring( this->operator[](i), o, sep, l, left );
                }
            }

            // otherwise, save full vector
            else {
                o << "[";
                for (size_t i=0; i<size(); ++i)
                {
                    if (i > 0)
                    {
                        o << ",";
                    }
                    Printer<valueType, IsDerivedFrom<valueType, Printable>::Is >::printForComplexStoring( this->operator[](i), o, sep, l, left );
                }
                o << "]";
            }

        }

    protected:

        // private members
        std::vector<valueType*>                             values;
    };


}


#endif

