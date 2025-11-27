#include "BetaBrokenStickFunction.h"

#include "DistributionBeta.h"
#include "RbMathFunctions.h"

/**
 * Wrapper for dealing with BetaBrokenStickFunction used the TypedDagNode classes of type doubles
 * @param a double value
 * @param b double value
 * @param nc number of categories to usen
 */

RevBayesCore::BetaBrokenStickFunction::BetaBrokenStickFunction(const TypedDagNode<double> *a, const TypedDagNode<double> *b, const TypedDagNode<std::int64_t> *nc) : TypedFunction< Simplex >( new Simplex() ),
alpha( a ),
beta( b ),
numCats(nc)
{

    addParameter( alpha );
    addParameter( beta );
    addParameter( numCats );

}



RevBayesCore::BetaBrokenStickFunction* RevBayesCore::BetaBrokenStickFunction::clone( void ) const
{

    return new BetaBrokenStickFunction(*this);
}


void RevBayesCore::BetaBrokenStickFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP) {

    if (oldP == alpha)
    {
        alpha = static_cast<const TypedDagNode<double>* >( newP );
    }

    if (oldP == beta)
    {
        beta = static_cast<const TypedDagNode<double>* >( newP );
    }

    if (oldP == numCats)
    {
        numCats = static_cast<const TypedDagNode<std::int64_t>* >( newP );
    }

}

void RevBayesCore::BetaBrokenStickFunction::update( void ) {

    double a = alpha->getValue();
    double b = beta->getValue();

    bool flipped = false;

    // if (b > a)
    // {
    //   flipped = true;
    //   double tmp = a;
    //   a = b;
    //   b = tmp;
    // }

    int nCats = numCats->getValue();
    double interval = 1.0 / (nCats);

    double p_low = 0.0;
    double p_high = 0.0;

    // Empty current simplex (value is a std::vector<double>*)
    value->clear();
    (*value) = std::vector<double>(nCats,0.0);

    for (int i=0; i<nCats-1; i++)
    {
      p_high = RbStatistics::Beta::quantile(a, b, (i/double (nCats)) + interval);
      (*value)[i] = p_high - p_low;
      p_low = p_high;
    }
    p_high = 1.0;
    (*value)[nCats-1] = p_high - p_low;
    

    if ( flipped )
    {
      std::reverse((*value).begin(),(*value).end());
    }

}
