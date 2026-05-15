#include <cstdlib>
#include <cmath>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <utility>
#include <vector>
#include <tuple>

#include "DagNode.h"
#include "DebugMove.h"
#include "SliceSamplingMove.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "AbstractMove.h"
#include "RbOrderedSet.h"
#include "StochasticNode.h"
#include "RbSettings.h" // for debugMCMC setting
#include <range/v3/all.hpp> // for ranges::views

namespace views = ranges::views;

using std::optional;

using namespace RevBayesCore;


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

        optional<LogDensity> operator()()
            {
                LogDensity lnPrior = 0.0;
                LogDensity lnLikelihood = 0.0;

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
                LogDensity lnPosterior = pHeat * (lHeat * lnLikelihood + prHeat * lnPrior);

                return lnPosterior;
            }

        optional<LogDensity> operator()(double x)
            {
                assert( not variable->isClamped() );

                num_evals++;

                // Don't set the variable to out-of-bounds values.
                if (below_lower_bound(x) or above_upper_bound(x))
                    return {};

		// NOTE: We could call variable->setValue(new double(x)),
		//       which calls distribution->setValue(x,true) and
		//       touches all the nodes.
		//       But this does less allocation.
                variable->getValue() = x;
                // touch all the nodes to set flags for recomputation
                variable->touch();

                auto Pr_ = (*this)();

                // call accept for each node  --  automatically includes affected nodes
                variable->keep();

                return Pr_;
            }

        NodePrMap getNodePrs()
        {
            return ::getNodePrs({variable}, affectedNodes);
        }

        void showNodePrs(const std::string& stage, const NodePrMap& pdfs) const
        {
            std::vector<const RevBayesCore::DagNode*> nodes;
            nodes.push_back(variable);
            ::showNodePrs(stage, pdfs, nodes);
        }

        void checkPrs()
        {
            double x = current_value();
            auto Prs = getNodePrs();
            auto Pr = operator()();
            auto Prx = operator()(x);
            auto Prsx = getNodePrs();
            compareNodePrs("SliceSampling(" + variable->getName() + ")", Prs, Prsx, "initial check");
        }

        double current_value() const
            {
                return variable->getValue();
            }

        std::string name() const {return variable->getName();}

        slice_function(StochasticNode<double> *n, optional<double> lb, optional<double> ub, double pr, double l, double p)
            :interval(lb,ub),
             variable(n),
             lHeat(l),
             pHeat(p),
             prHeat(pr),
             num_evals(0)
            {
                variable->initiateGetAffectedNodes( affectedNodes );
            }
    };

}

std::pair<double,double> 
find_slice_boundaries_stepping_out(double x0,slice_function& g, const LogDensity& logy, double w,int m)
{
    int logMCMC = RbSettings::userSettings().getLogMCMC();

    assert(x0 + w > x0);
    assert(g.in_range(x0));

    double u = uniform()*w;
    double L = x0 - u;
    double R = x0 + (w-u);
    assert(L < x0);
    assert(x0 < R);

    if (logMCMC >= 4) std::cerr<<"stepping:     L = "<<L<<"  x0 = "<<x0<<"   R0 = "<<R<<"\n";

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

std::ostream& operator<<(std::ostream& o, optional<LogDensity> l)
{
    if (not l)
        o<<"OOB";
    else
        o<<l.value();
    return o;
}

std::tuple<double,double,optional<optional<LogDensity>>,optional<optional<LogDensity>>>
find_slice_boundaries_doubling(double x0,slice_function& g, const LogDensity& logy, double w, int K)
{
    int logMCMC = RbSettings::userSettings().getLogMCMC();
    assert(x0 + w > x0);
    assert(g.in_range(x0));

    double u = uniform()*w;
    double L = x0 - u;
    double R = L + w;
    assert(L < x0);
    assert(x0 < R);

    optional<optional<LogDensity>> gL_cached;
    auto gL = [&]() {
        if (not gL_cached)
            gL_cached = g(L);
        return *gL_cached;
    };

    optional<optional<LogDensity>> gR_cached;
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
        if (logMCMC >= 4)
            std::cerr<<"doubling:    L0 = "<<L<<" (g(L) = "<<gL()<<")  x0 = "<<x0<<"   R0 = "<<R<<" (g(R) = "<<gR()<<")\n";

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

    if (logMCMC >= 1)
        std::cerr<<"doubling:    L0 = "<<L<<" (g(L) = "<<gL()<<")  x0 = "<<x0<<"   R0 = "<<R<<" (g(R) = "<<gR()<<")\n";

    return {L,R,gL_cached,gR_cached};
}

double search_interval(double x0,double& L, double& R, slice_function& g,const LogDensity& logy)
{
    int logMCMC = RbSettings::userSettings().getLogMCMC();
    int debugMCMC = RbSettings::userSettings().getDebugMCMC();

    if (debugMCMC >= 1)
	g.checkPrs();

//    Don't compute g(L) or g(R) because they might not be in the domain.
//    Arguably they should return Pr=0 and that would be fine.
//    But right now we have some asserts that trigger.
//    assert(g(x0) > g(L) and g(x0) > g(R));
    assert(g(x0) >= logy);
    assert(L < R);
    assert(L <= x0 and x0 <= R);

    //double L0 = L, R0 = R;

    for (int i=0;i<200;i++)
    {
        double x1 = L + uniform()*(R-L);

        auto gx1 = g(x1);

        if (logMCMC >= 4)
	    std::cerr<<"search:    L0 = "<<L<<"  x1 = "<<x1<<" g(x1) = "<<gx1<<"  R0 = "<<R<<"\n";

        if (gx1 and *gx1 >= logy)
        {
            assert(g.in_range(x1));
            return x1;
        }
                

        if (x1 > x0) 
            R = x1;
        else
            L = x1;
    }

    std::abort();

    assert(g.in_range(x0));
    return x0;
}


bool pre_slice_sampling_check_OK(double x0, slice_function& g)
{
    // If x is not in the range then this could be a range that is reduced to avoid loss of precision.
    if (not g.in_range(x0))
    {
        int logMCMC = RbSettings::userSettings().getLogMCMC();
        if (logMCMC >= 4) std::cerr<<"   "<<x0<<" not in range!";
        return false;
    }
    else
    {
        assert(g.in_range(x0));
        return true;
    }
}

bool can_propose_same_interval_doubling(double x0, double x1, double w, double L, double R, optional<optional<LogDensity>> gL_cached, optional<optional<LogDensity>> gR_cached, slice_function& g, const LogDensity& log_y)
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

    int logMCMC = RbSettings::userSettings().getLogMCMC();
    int debugMCMC = RbSettings::userSettings().getDebugMCMC();

    NodePrMap initialPdfs;
    if (logMCMC >= 4)
    	initialPdfs = g.getNodePrs();
    if (logMCMC >= 4)
        g.showNodePrs("BEFORE", initialPdfs);

    // 1. Determine the slice level, in log terms.
    auto gx0 = g();

    if (not gx0)
        throw RbException()<<"slice_sampling_stepping_out: x0 is out of bounds!";

    // Determine the slice level, in log terms.
    auto logy = *gx0 + log(uniform()); // - exponential(1.0);

    if (logMCMC >= 1)
	std::cerr<<"mvSlice("<<g.name()<<", stepping_out):  x0 = "<<x0<<"  g(x0) = "<<gx0<<"   logy = "<<logy<<"\n";

    // Check that g() and g(x) don't give different PDFs!
    if (debugMCMC >= 1 or debug_build)
        g.checkPrs();

    assert(*g(x0) > logy);

    // 2. Find the initial interval to sample from.
    auto [L, R] = find_slice_boundaries_stepping_out(x0,g,logy,w,m);

    // 3. Sample from the interval, shrinking it on each rejection
    double x = search_interval(x0,L,R,g,logy);

    NodePrMap finalPdfs;
    if (logMCMC >= 4)
	finalPdfs = g.getNodePrs();
    if (logMCMC >= 4)
        g.showNodePrs("FINAL", finalPdfs);

    return x;
}


// We need to SET the value INSIDE this routine.
// Are we assuming that calling g sets the value?
double slice_sample_doubling(double x0, slice_function& g, double w, int m)
{
    // 0. Check that the values are OK
    if (not pre_slice_sampling_check_OK(x0, g))
        return x0;

    int logMCMC = RbSettings::userSettings().getLogMCMC();
    int debugMCMC = RbSettings::userSettings().getDebugMCMC();

    NodePrMap initialPdfs;
    if (logMCMC >= 3 or debugMCMC >= 1)
    	initialPdfs = g.getNodePrs();
    if (logMCMC >= 3)
        g.showNodePrs("BEFORE", initialPdfs);

    // 1. Determine the slice level, in log terms.
    auto gx0 = g();

    if (not gx0)
        throw RbException()<<"slice_sampling_stepping_out: x0 is out of bounds!";

    // Determine the slice level, in log terms.
    LogDensity logy = *gx0 + log(uniform()); // - exponential(1);

    if (logMCMC >= 1) std::cerr<<"mvSlice("<<g.name()<<", doubling):  x0 = "<<x0<<"  g(x0) = "<<gx0<<"   logy = "<<logy<<"\n";

    // Check that g() and g(x) don't give different PDFs!
    if (debugMCMC >= 1 or debug_build)
        g.checkPrs();
    assert(*g(x0) > logy);

    // 2. Find the initial interval to sample from.
    auto [L, R, gL_cached, gR_cached]  = find_slice_boundaries_doubling(x0,g,logy,w,m);

    // 3. Sample from the interval, shrinking it on each rejection
    double x1 = search_interval(x0,L,R,g,logy);

    // 4. Check that we can propose the same interval from x2
    // We need to SET the value INSIDE this routine if we recompute g().
    if (not can_propose_same_interval_doubling(x0, x1, w, L, R, gL_cached, gR_cached, g, logy))
        x1 = x0;

    NodePrMap finalPdfs;
    if (logMCMC >=4)
	finalPdfs = g.getNodePrs();
    if (logMCMC >= 4)
        g.showNodePrs("FINAL", finalPdfs);

    return x1;
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

