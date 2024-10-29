#include <cstdlib>
#include <cmath>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <utility>
#include <vector>
#include <tuple>

#include <boost/optional.hpp>

#include "DagNode.h"
#include "SliceSamplingMove.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "AbstractMove.h"
#include "RbOrderedSet.h"
#include "StochasticNode.h"
#include "RbSettings.h" // for debugMCMC setting

using boost::optional;

using namespace RevBayesCore;

const double log_0 = RbConstants::Double::neginf;

// How can we make a run-time adjustable log level?
int log_level = 0;

/** 
 * Constructor
 *
 * Here we simply allocate and initialize the move object.
 *
 * \param[in]    w   The weight how often the proposal will be used (per iteration).
 * \param[in]    t   If auto tuning should be used.
 */
SliceSamplingMove::SliceSamplingMove( StochasticNode<double> *n, optional<double> lb, optional<double> ub, double window_, double weight_, BoundarySearchMethod s, bool t )
    : AbstractMove( std::vector<DagNode*>(), weight_ ,t),
      lower_bound(lb),
      upper_bound(ub),
      variable( n ),
      window( window_ ),
      total_movement( 0.0 ),
      search_method ( s )
{
    assert( not variable->isClamped() );

    addNode( n );
}


/**
 * Basic destructor doing nothing.
 */
SliceSamplingMove::~SliceSamplingMove( void )
{
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the SliceSamplingMove. 
 */
SliceSamplingMove* SliceSamplingMove::clone( void ) const 
{
    return new SliceSamplingMove( *this );
}


/**
 * Get moves' name of object 
 *
 * \return The moves' name.
 */
const std::string& SliceSamplingMove::getMoveName( void ) const 
{
    static std::string name = "SliceSampling";

    return name;
}


double SliceSamplingMove::getMoveTuningParameter( void ) const
{
    return window;
}


double uniform()
{
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    return rng->uniform01();
}

struct interval
{
    optional<double> lower_bound;

    optional<double> upper_bound;

    /// check if a value is below the lower bound on the range
    bool below_lower_bound(double x) const { return (lower_bound and x < *lower_bound); }
    /// check if a value is above the upper bound on the range
    bool above_upper_bound(double x) const { return (upper_bound and x > *upper_bound); }
    /// check if a value is in the range or not
    bool in_range(double x) const  { return (not below_lower_bound(x) and not above_upper_bound(x));}

    interval() {};
    interval(optional<double> lb, optional<double> ub): lower_bound(lb), upper_bound(ub) {}
};

std::map<const DagNode*, double> getNodePrs(DagNode* node, const RbOrderedSet<DagNode*>& affectedNodes)
{
    std::map<const DagNode*, double> Prs;
    Prs.insert({node, node->getLnProbability()});
    for(auto affectedNode: affectedNodes)
	Prs.insert({affectedNode, affectedNode->getLnProbability()});
    return Prs;
}

namespace  {

/// This object allow computing the probability of the current point, and also store the variable's range
    class slice_function: public interval
    {
        StochasticNode<double>* variable;
        double lHeat;
        double pHeat;
        double prHeat;
        RbOrderedSet<DagNode*> affectedNodes;
        int num_evals;

    public:

        int get_num_evals() const {return num_evals;}

        double operator()()
            {
                double lnPrior = 0.0;
                double lnLikelihood = 0.0;

                // 1. compute the probability of the current value for each node
                lnPrior += variable->getLnProbability();

                // 2. then we recompute the probability for all the affected nodes
                for (auto& node: affectedNodes)
                {
                    if ( node->isClamped() )
                        lnLikelihood += node->getLnProbability();
                    else
                        lnPrior += node->getLnProbability();
                }

                // 3. exponentiate with the chain heat
                double lnPosterior = pHeat * (lHeat * lnLikelihood + prHeat * lnPrior);

                return lnPosterior;
            }

        double operator()(double x)
            {
                assert( not variable->isClamped() );

                num_evals++;

                // Don't set the variable to out-of-bounds values.
                if (below_lower_bound(x) or above_upper_bound(x))
                    return RbConstants::Double::neginf;

		// NOTE: We could call variable->setValue(new double(x)),
		//       which calls distribution->setValue(x,true) and
		//       touches all the nodes.
		//       But this does less allocation.
                variable->getValue() = x;
                // touch all the nodes to set flags for recomputation
                variable->touch();

                double Pr_ = (*this)();

                // call accept for each node  --  automatically includes affected nodes
                variable->keep();

                return Pr_;
            }

        double current_value() const
            {
                return variable->getValue();
            }

        slice_function(StochasticNode<double> *n, optional<double> lb, optional<double> ub, double pr, double l, double p)
            :interval(lb,ub),
             variable(n),
             lHeat(l),
             pHeat(p),
             prHeat(pr),
             num_evals(0)
            {
                variable->initiateGetAffectedNodes( affectedNodes );

                if (RbSettings::userSettings().getDebugMCMC() > 0)
                {
                    std::map<const DagNode*, double> NodePrs;

                    double x = current_value();
                    double Pr = operator()();
                    auto Prs = getNodePrs(variable,affectedNodes);
                    double Prx = operator()(x);
                    auto Prsx = getNodePrs(variable,affectedNodes);
                    if (std::abs(Pr - Prx) > 1.0e-9)
                    {
                        std::cerr<<std::setprecision(10)<<std::endl;
                        std::cerr<<"mvSlice for "<<variable->getName()<<": probability is "<<Pr<<" but should be "<<Prx<<":  delta = "<<Pr - Prx<<"\n";
                        for(auto& [n,pr1]: Prs)
                        {
                            double pr2 = Prsx.at(n);
                            if (std::abs(pr1-pr2) > 1.0e-6)
                                std::cerr<<"         cause: probability for "<<n->getName()<<" is "<<pr1<<" but should be "<<pr2<<":  delta = "<<pr1-pr2<<"\n";
                        }
                        std::abort();
                    }
                }
            }
    };

}

std::pair<double,double> 
find_slice_boundaries_stepping_out(double x0,slice_function& g,double logy, double w,int m)
{
    assert(g.in_range(x0));

    double u = uniform()*w;
    double L = x0 - u;
    double R = x0 + (w-u);

    // Expand the interval until its ends are outside the slice, or until
    // the limit on steps is reached.

    if (m>1) {
        int J = uniform()*m;
        int K = (m-1)-J;

        while (J>0 and (not g.below_lower_bound(L)) and g(L)>logy) {
            L -= w;
            J--;
        }

        while (K>0 and (not g.above_upper_bound(R)) and g(R)>logy) {
            R += w;
            K--;
        }
    }
    else {
        while ((not g.below_lower_bound(L)) and g(L)>logy)
            L -= w;

        while ((not g.above_upper_bound(R)) and g(R)>logy)
            R += w;
    }

    // Shrink interval to lower and upper bounds.

    if (g.below_lower_bound(L)) L = *g.lower_bound;
    if (g.above_upper_bound(R)) R = *g.upper_bound;

    assert(L < R);

    return {L,R};
}

std::tuple<double,double,optional<double>,optional<double>>
find_slice_boundaries_doubling(double x0,slice_function& g,double logy, double w, int K)
{
    assert(x0 + w > x0);
    assert(g.in_range(x0));

    double u = uniform()*w;
    double L = x0 - u;
    double R = L + w;
    assert(L < x0);
    assert(x0 < R);

    optional<double> gL_cached;
    auto gL = [&]() {
        if (not gL_cached)
            gL_cached = g(L);
        return *gL_cached;
    };

    optional<double> gR_cached;
    auto gR = [&]() {
        if (not gR_cached)
            gR_cached = g(R);
        return *gR_cached;
    };

    auto too_large = [](double L, double R, double w)
        {
            double M = (L+R)/2;
            double W = R-L;
            assert(W > 0);
            bool ok = (L < M) and (M < R) and (W+w > W) and (L-w < L) and (R+w>R);
            return not ok;
        };

    while ( K > 0 and (gL() > logy or gR() > logy))
    {
        if (log_level >= 4)
            std::cerr<<"!!    L0 = "<<L<<" (g(L) = "<<gL()<<")  x0 = "<<x0<<"   R0 = "<<R<<" (g(R) = "<<gR()<<")\n";

        double W2 = (R-L);
        if (uniform() < 0.5)
        {
            double L2 = L - W2;
            if (too_large(L2, R, w))
                break;
            L = L2;
            gL_cached = {};
        }
        else
        {
            double R2 = R + W2;
            if (too_large(L, R2, w))
                break;
            R = R2;
            gR_cached = {};
        }

        K --;
    }

    assert(L < R);
    assert( L < (L+R)/2 and (L+R)/2 < R);

    //  std::cerr<<"[]    L0 = "<<L<<"   x0 = "<<x0<<"   R0 = "<<R<<"\n";

    // FIXME: GCC 5 complains if we don't write out the tuple type. GCC 7 does not need it.  How about GCC 6?
    return {L,R,gL_cached,gR_cached};
}

double search_interval(double x0,double& L, double& R, slice_function& g,double logy)
{
    //  assert(g(x0) > g(L) and g(x0) > g(R));
    assert(g(x0) >= logy);
    assert(L < R);
    assert(L <= x0 and x0 <= R);

    //double L0 = L, R0 = R;

    for (int i=0;i<200;i++)
    {
        double x1 = L + uniform()*(R-L);
        double gx1 = g(x1);

        if (gx1 >= logy) return x1;

        if (x1 > x0) 
            R = x1;
        else
            L = x1;
    }

    std::abort();

    return x0;
}


bool pre_slice_sampling_check_OK(double x0, slice_function& g)
{
    // If x is not in the range then this could be a range that is reduced to avoid loss of precision.
    if (not g.in_range(x0))
    {
        if (log_level >= 4) std::cerr<<x0<<" not in range!";
        return false;
    }
    else
    {
        assert(g.in_range(x0));
        return true;
    }
}

bool can_propose_same_interval_doubling(double x0, double x1, double w, double L, double R, optional<double> gL_cached, optional<double> gR_cached, slice_function& g, double log_y)
{
    bool D = false;

    auto gL = [&]() {
        if (not gL_cached)
            gL_cached = g(L);
        return *gL_cached;
    };

    auto gR = [&]() {
        if (not gR_cached)
            gR_cached = g(R);
        return *gR_cached;
    };

    bool ok = true;
    while (ok and R-L > 1.1*w)
    {
        double M = (R+L)/2;
        assert( L < M and M < R);

        // Check if x0 and x1 are in different halves of the interval.
        if ((x0 < M and x1 >= M) or (x0 >= M and x1 < M))
            D = true;

        if (x1 < M)
        {
            R = M;
            gR_cached = {};
        }
        else
        {
            L = M;
            gL_cached = {};
        }

        if (D and log_y >= gL() and log_y >= gR())
            ok = false;
    }

    // FIXME - this is clunky.  Do we really want to set x by evaluate g( )?
    if (D)
    {
        // We may have set x to L or R, so reset it to the right values.
        if (ok)
            g(x1);
        else
            g(x0);
    }

    return ok;
}

double slice_sample_stepping_out(double x0, slice_function& g,double w, int m)
{
    assert(g.in_range(x0));

    double gx0 = g();
#ifndef NDEBUG
    volatile double diff = gx0 - g(x0);
    assert(std::abs(diff) < 1.0e-9);
#endif

    // Determine the slice level, in log terms.

    double logy = gx0 + log(uniform()); // - exponential(1.0);

    // Find the initial interval to sample from.
    auto [L, R] = find_slice_boundaries_stepping_out(x0,g,logy,w,m);

    // Sample from the interval, shrinking it on each rejection

    return search_interval(x0,L,R,g,logy);
}


// We need to SET the value INSIDE this routine.
// Are we assuming that calling g sets the value?
double slice_sample_doubling(double x0, slice_function& g, double w, int m)
{
    // 0. Check that the values are OK
    if (not pre_slice_sampling_check_OK(x0, g))
        return x0;

    // 1. Determine the slice level, in log terms.
    double logy = g() + log(uniform()); // - exponential(1);

    // 2. Find the initial interval to sample from.
    auto [L, R, gL_cached, gR_cached]  = find_slice_boundaries_doubling(x0,g,logy,w,m);

    // 3. Sample from the interval, shrinking it on each rejection
    double x1 = search_interval(x0,L,R,g,logy);

    // 4. Check that we can propose the same interval from x2
    // We need to SET the value INSIDE this routine if we recompute g().
    if (can_propose_same_interval_doubling(x0, x1, w, L, R, gL_cached, gR_cached, g, logy))
        return x1;
    else
        return x0;
}

void SliceSamplingMove::performMcmcMove( double prHeat, double lHeat, double pHeat )
{
    slice_function g(variable, lower_bound, upper_bound, prHeat, lHeat, pHeat);

    double x1 = g.current_value();

    double x2;
    if (search_method == search_doubling)
        x2 = slice_sample_doubling(x1, g, window, 40);
    else
        x2 = slice_sample_stepping_out(x1, g, window, 40);

    total_movement += std::abs(x2 - x1);

    numPr += g.get_num_evals();

    if (auto_tuning and num_tried_total > 3 and search_method == search_stepping_out)
    {
        double predicted_window = 4.0*total_movement/num_tried_total;
        window = 0.95*window + 0.05*predicted_window;
    }
    
}


/**
 * Print the summary of the move.
 *
 * The summary just contains the current value of the tuning parameter.
 * It is printed to the stream that it passed in.
 *
 * \param[in]     o     The stream to which we print the summary.
 */
void SliceSamplingMove::printSummary(std::ostream &o, bool current_period) const
{
    std::streamsize previousPrecision = o.precision();
    std::ios_base::fmtflags previousFlags = o.flags();
    
    o << std::fixed;
    o << std::setprecision(4);
    
    // print the name
    const std::string &n = getMoveName();
    size_t spaces = 40 - (n.length() > 40 ? 40 : n.length());
    o << n;
    for (size_t i = 0; i < spaces; ++i) {
        o << " ";
    }
    o << " ";
    
    // print the DagNode name
    const std::string &dn_name = (*nodes.begin())->getName();
    spaces = 20 - (dn_name.length() > 20 ? 20 : dn_name.length());
    o << dn_name;
    for (size_t i = 0; i < spaces; ++i) {
        o << " ";
    }
    o << " ";
    
    // print the weight
    int w_length = 4;
    if (weight > 0) w_length -= (int)log10(weight);
    for (int i = 0; i < w_length; ++i) {
        o << " ";
    }
    o << weight;
    o << " ";
    
    // print the number of tries
    int t_length = 9;
    if (num_tried_total > 0) t_length -= (int)log10(num_tried_total);
    for (int i = 0; i < t_length; ++i) {
        o << " ";
    }
    o << num_tried_total;
    o << " ";
    
    // print the average distance moved
    o<<"\n";
    if (num_tried_current_period > 0)
    {
        o<<"  Ave. |x2-x1| = "<<total_movement/num_tried_current_period<<std::endl;
    }

    // print the average distance moved
    if (num_tried_current_period > 0)
    {
        o<<"  Ave. # of Pr evals = "<<double(numPr)/num_tried_current_period<<std::endl;
    }

    //    proposal->printParameterSummary( o );
    o<<"  window = "<<window<<std::endl;
    
    o << std::endl;
    
    o.setf(previousFlags);
    o.precision(previousPrecision);
    
    
}

/**
 * Reset the move counters. Here we only reset the counter for the number of accepted moves.
 *
 */
void SliceSamplingMove::resetMoveCounters( void )
{
    total_movement = 0.0;
    num_tried_current_period = 0;
    numPr = 0;
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void SliceSamplingMove::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    variable = static_cast<StochasticNode<double>* >(newN) ;
}


void SliceSamplingMove::setMoveTuningParameter(double tp)
{
    window = tp;
}


/**
 * Tune the move to accept the desired acceptance ratio.
 * We only compute the acceptance ratio here and delegate the call to the proposal.
 */
void SliceSamplingMove::tune( void ) 
{
    
    if ( num_tried_current_period > 2 )
    {
        double predicted_window = 4.0*total_movement/num_tried_current_period;
        
        double p = exp(-double(num_tried_current_period)*0.5);
        window = p*window + (1.0-p)*predicted_window;
    }
    
}

