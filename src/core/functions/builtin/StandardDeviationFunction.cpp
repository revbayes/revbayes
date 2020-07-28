#include "StandardDeviationFunction.h"

#include <cmath>
#include <vector>

#include "RbVector.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/**
 * StandardDeviationFunction of a RbVector Constructor.
 * @param v the vector of values of type double
 */
StandardDeviationFunction::StandardDeviationFunction(const TypedDagNode< RbVector<double> > *v) : TypedFunction<double>( new double(0.0) ),
vals( v )
{
    // add the parameters as parents
    this->addParameter( vals );
    
    update();
}


/**
 * Empty destructor.
 */
StandardDeviationFunction::~StandardDeviationFunction( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



/**
 * Clone function for deep copies.
 */
StandardDeviationFunction* StandardDeviationFunction::clone( void ) const {
    return new StandardDeviationFunction( *this );
}


/**
 * Update the current value based on the current parameter values.
 */
void StandardDeviationFunction::update( void ) {
    
    double m = 0;
    const std::vector<double> &v = vals->getValue();
    for ( std::vector<double>::const_iterator it = v.begin(); it != v.end(); ++it) {
        m += *it;
    }
    m /= v.size();
    
    double var = 0.0;
    for ( std::vector<double>::const_iterator it = v.begin(); it != v.end(); ++it) {
        var += (*it-m)*(*it-m);
    }
    
    *this->value = sqrt( var / (int(v.size())-1) );
    
}



/**
 * Swap the internal parameters. This happens when the parameters are re-assigned or the entire model graph is cloned.
 * Here we only need to store the new pointer to the vector of real values.
 */
void StandardDeviationFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if ( oldP == vals )
    {
        vals = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    
}

