#include <cstddef>
#include <sstream>

#include "MaximumLikelihoodEstimation.h"
#include "Parallelizable.h"

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



