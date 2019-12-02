#include "DagNode.h"
#include "VariableMonitor.h"
#include "MaximumLikelihoodEstimation.h"
#include "MoveSchedule.h"
#include "RandomMoveSchedule.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbException.h"
#include "RbIterator.h"
#include "RbMathLogic.h"
#include "RbOptions.h"
#include "RlUserInterface.h"
#include "SingleRandomMoveSchedule.h"
#include "RandomMoveSchedule.h"
#include "ExtendedNewickTreeMonitor.h"

#include <cmath>
#include <iomanip>
#include <sstream>
#include <typeinfo>
#include "SequentialMoveSchedule.h"

using namespace RevBayesCore;


/**
 * Constructor. We create an independent copy of the model and thus of all DAG nodes.
 * Someone might have wanted to run another MaximumLikelihoodEstimation with different settings on the same model.
 * Thus we also create our own copies of the monitors and moves.
 *
 * \param[in]    m    The model containing all DAG nodes.
 * \param[in]    mvs  The vector of moves.
 * \param[in]    mons The vector of monitors.
 */
MaximumLikelihoodEstimation::MaximumLikelihoodEstimation(void) : Parallelizable(),
generation(0)
{
    
}


/**
 * Destructor. Nothing to do here.
 */
MaximumLikelihoodEstimation::~MaximumLikelihoodEstimation(void)
{
    
}


/**
 * Get the current generation number.
 */
size_t MaximumLikelihoodEstimation::getCurrentGeneration( void ) const
{
    return generation;
}


std::ostream& RevBayesCore::operator<<(std::ostream& o, const MaximumLikelihoodEstimation& x)
{
    o << "MaximumLikelihoodEstimation";
    
    return o;
}



