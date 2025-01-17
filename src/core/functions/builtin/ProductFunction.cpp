//
//  ProductFunction.cpp
//  revbayes_traitfig
//
//  Created by Michael Landis on 11/15/24.
//

#include "ProductFunction.h"
#include "RbConstants.h"
#include "RbConstIterator.h"
#include "RbConstIteratorImpl.h"
#include "RbVector.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;
/**
 * ProductFunction of a RbVector Constructor.
 * @param v the vector of values of type double
 */
ProductFunction::ProductFunction(const TypedDagNode<RbVector<double> > *v)
    : TypedFunction<double>( new double(0.0) ),
    vals( v ),
    use_safe_product(true)
{
    // add the parameters as parents
    this->addParameter( vals );
    
    update();
}


ProductFunction::~ProductFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



ProductFunction* ProductFunction::clone( void ) const
{
    return new ProductFunction( *this );
}


void ProductFunction::update( void )
{
    
    double m = 1.0;
    const RbVector<double> &v = vals->getValue();
    if (use_safe_product) {
  
        size_t n = v.size();
  
        m = 0.0;
        int sign = 1;
        bool containsZero = false;
        for ( RbConstIterator<double> it = v.begin(); it != v.end(); ++it) {
            // abort if val == 0
            if (*it == 0.0) {
                containsZero = true;
                break;
            }
            // track sign of result
            else if (*it < 0) {
                sign *= -1;
            }
            // sum-log
            m += std::log(std::abs(*it));
        }

        if (containsZero) {
            m = 0.0;
        } else {
            m = sign * std::exp(m);
        }
        
    } else {
        for ( RbConstIterator<double> it = v.begin(); it != v.end(); ++it)
        {
            if (*it == 0.0) {
                // product is zero if any element is zero
                m = 0.0;
                break;
            } else {
                // otherwise, multiply element into product
                m *= *it;
            }
        }
    }
        
    *this->value = m ;
    
}



void ProductFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if ( oldP == vals )
    {
        vals = static_cast<const TypedDagNode<RbVector<double> >* >( newP );
    }
    
}
