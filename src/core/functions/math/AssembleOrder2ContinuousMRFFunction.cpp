#include "AssembleOrder2ContinuousMRFFunction.h"

#include <cmath>

#include "Cloner.h"
#include "RbVectorImpl.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

AssembleOrder2ContinuousMRFFunction::AssembleOrder2ContinuousMRFFunction(const TypedDagNode< double > *theta1, const TypedDagNode< RbVector<double> > *increments, bool log_theta1) : TypedFunction< RbVector<double> >( new RbVector<double>(increments->getValue().size()+1,0.0) ),
    increments( increments ),
    theta1( theta1 ),
    theta_1_is_log( log_theta1 ),
    field_size( increments->getValue().size()+1 )
{
    // add the lambda parameter as a parent
    addParameter( increments );
    addParameter( theta1 );

    update();
}

AssembleOrder2ContinuousMRFFunction::~AssembleOrder2ContinuousMRFFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



AssembleOrder2ContinuousMRFFunction* AssembleOrder2ContinuousMRFFunction::clone( void ) const
{
    return new AssembleOrder2ContinuousMRFFunction( *this );
}


void AssembleOrder2ContinuousMRFFunction::update( void )
{

  RbVector<double> &log_theta = *value;

  RbVector<double> delta = increments->getValue();

  // Handle first value being on log-scale or not
  if (theta_1_is_log)
  {
    log_theta[0] = theta1->getValue();
  }
  else
  {
    log_theta[0] = log(theta1->getValue());
  }

  // Assemble the field
  log_theta[1] = log_theta[0] + delta[0];

  for (size_t i=0; i<(field_size-2); ++i)
  {
    log_theta[i+2] = -log_theta[i] + 2*log_theta[i+1] + delta[i+1];
  }

  // Transform off log-scale
  for (size_t i=0; i<field_size; ++i)
  {
    log_theta[i] = exp(log_theta[i]);
  }

}

void AssembleOrder2ContinuousMRFFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{

    if (oldP == increments)
    {
        increments = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == theta1)
    {
        theta1 = static_cast<const TypedDagNode< double >* >( newP );
    }

}
