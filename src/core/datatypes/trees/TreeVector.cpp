#include "TreeVector.h"

#include <ostream>
#include <string>

#include "RbException.h"

using namespace RevBayesCore;

TreeVector::TreeVector(void)
{

}


bool TreeVector::operator==(const TreeVector &mve) const
{

    if ( num_trees != mve.num_trees )
    {
        return false;
    }

    if ( values.size() != mve.values.size() )
    {
        return false;
    }

    for (size_t i=0; i<values.size(); ++i)
    {
        if ( values[i] != mve.values[i] )
        {
            return false;
        }
    }

    return true;
}



bool TreeVector::operator!=(const TreeVector &mve) const
{
    return (operator==(mve) == false);
}



void TreeVector::addValues(const Tree *t, long n)
{
    for (size_t i=0; i<(size_t)n; i++)
    {
        Tree *tree_copy = t->clone();
        values.push_back( *tree_copy );
    }
}



void TreeVector::clear(void)
{
    values.clear();
}



TreeVector* TreeVector::clone(void) const
{
    return new TreeVector( *this );
}



size_t TreeVector::getNumberOfValues(void) const
{
    return values.size();
}



std::vector<Tree>& TreeVector::getValues(void)
{
    return values;
}



const std::vector<Tree>& TreeVector::getValues(void) const
{
    return values;
}



void TreeVector::setNumberOfTrees(long n)
{
    num_trees = n;
}



std::ostream& RevBayesCore::operator<<(std::ostream& o, const RevBayesCore::TreeVector& x)
{

    const std::vector<Tree>& trees = x.getValues();

    for ( size_t j=0; j<x.getNumberOfValues(); ++j )
    {
        o << trees[j].getNewickRepresentation();
    }

    return o;
}
