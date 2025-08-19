#include "AbstractGibbsMove.h"

#include <cstddef>
#include <cmath>
#include <iomanip>
#include <vector>

#include "DagNode.h"
#include "RbException.h"
#include "Cloneable.h"
#include "RbConstants.h" // IWYU pragma: keep


using namespace RevBayesCore;


/**
 * Constructor
 *
 * Here we simply allocate and initialize the move object.
 *
 * \param[in]    w   The weight how often the proposal will be used (per iteration).
 * \param[in]    t   If auto tuning should be used.
 */
AbstractGibbsMove::AbstractGibbsMove( double w  ) : AbstractMove( w, false )
{

}



/**
 * Basic destructor doing nothing.
 */
AbstractGibbsMove::~AbstractGibbsMove( void )
{

}


double AbstractGibbsMove::getMoveTuningParameter( void ) const
{
    // Gibbs move has no tuning parameter
    return RbConstants::Double::nan;
}

/**
 * Check if move is allowable with given heating scheme.
 * Gibbs moves are not generally compatible with heats, so here we return false with any heating.
 * We allow individual moves to override this on a case-by-case basis.
 */
bool AbstractGibbsMove::heatsAreAllowable(double prHeat, double lHeat, double pHeat)
{
    if ( prHeat != 1.0 || lHeat != 1.0 || pHeat != 1.0 )
    {
      return false;
    }
    else
    {
      return true;
    }
}


/**
 * Perform the move.
 * Here we store some info and delegate to performMove.
 */
void AbstractGibbsMove::performMcmcMove( double prHeat, double lHeat, double pHeat )
{
    // check heating values
    if ( !heatsAreAllowable(prHeat, lHeat, pHeat) )
    {
        throw RbException("Cannot apply Gibbs sampler when the probability is heated.");
    }

    // delegate to derived class
    performGibbsMove();

}


/**
 * Print the summary of the move.
 *
 * The summary just contains the current value of the tuning parameter.
 * It is printed to the stream that it passed in.
 *
 * \param[in]     o     The stream to which we print the summary.
 */
void AbstractGibbsMove::printSummary(std::ostream &o, bool current_period) const
{
    std::streamsize previousPrecision = o.precision();
    std::ios_base::fmtflags previousFlags = o.flags();

    o << std::fixed;
    o << std::setprecision(4);

    // print the name
    const std::string &n = getMoveName();
    size_t spaces = 40 - (n.length() > 40 ? 40 : n.length());
    o << n;
    for (size_t i = 0; i < spaces; ++i)
    {
        o << " ";
    }
    o << " ";

    // print the DagNode name
    const std::vector<DagNode*> nodes = getDagNodes();
    std::string dn_name = "???";
    if ( nodes.size() > 0 )
    {
        dn_name = nodes[0]->getName();
    }
    spaces = 20 - (dn_name.length() > 20 ? 20 : dn_name.length());
    o << dn_name;
    for (size_t i = 0; i < spaces; ++i)
    {
        o << " ";
    }
    o << " ";

    // print the weight
    int w_length = 4;
    if (weight > 0) w_length -= (int)log10(weight);
    for (int i = 0; i < w_length; ++i)
    {
        o << " ";
    }
    o << weight;
    o << " ";

    size_t num_tried = num_tried_total;
    if (current_period == true)
    {
        num_tried = num_tried_current_period;
    }

    // print the number of tries
    int t_length = 9;
    if (num_tried > 0) t_length -= (int)log10(num_tried);
    for (int i = 0; i < t_length; ++i)
    {
        o << " ";
    }
    o << num_tried;
    o << " ";

    // print the number of accepted
    int a_length = 9;
    if (num_tried > 0) a_length -= (int)log10(num_tried);

    for (int i = 0; i < a_length; ++i)
    {
        o << " ";
    }
    o << num_tried;
    o << " ";

    // print the acceptance ratio
    double ratio = 1.0;
    if (num_tried == 0) ratio = 0;
    int r_length = 5;

    for (int i = 0; i < r_length; ++i)
    {
        o << " ";
    }
    o << ratio;
    o << " ";

    o << std::endl;

    o.setf(previousFlags);
    o.precision(previousPrecision);


}


void AbstractGibbsMove::setMoveTuningParameter(double tp)
{
    // Gibbs move has no tuning parameter: nothing to do
}


/**
 * Tune the move.
 * This is a dummy implementation because Gibbs move cannot be tuned.
 */
void AbstractGibbsMove::tune( void )
{

}
