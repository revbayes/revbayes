#include <stdlib.h>
#include <cmath>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <utility>
#include <vector>

#include "DagNode.h"
#include "SliceSamplingMove.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "AbstractMove.h"
#include "RbOrderedSet.h"
#include "StochasticNode.h"

using namespace RevBayesCore;


/** 
 * Constructor
 *
 * Here we simply allocate and initialize the move object.
 *
 * \param[in]    w   The weight how often the proposal will be used (per iteration).
 * \param[in]    t   If auto tuning should be used.
 */
SliceSamplingMove::SliceSamplingMove( StochasticNode<double> *n, double window_, double weight_, bool t ) : AbstractMove( std::vector<DagNode*>(), weight_ ,t),
    variable( n ),
    window( window_ ),
    total_movement( 0.0 )
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
  bool has_lower_bound;
  double lower_bound;

  bool has_upper_bound;
  double upper_bound;

  /// check if a value is below the lower bound on the range
  bool below_lower_bound(double x) const { return (has_lower_bound and x<lower_bound); }
  /// check if a value is above the upper bound on the range
  bool above_upper_bound(double x) const { return (has_upper_bound and x>upper_bound); }
  /// check if a value is in the range or not
  bool in_range(double x) const  { return (not below_lower_bound(x) and not above_upper_bound(x));}

  interval():has_lower_bound(false),lower_bound(0),has_upper_bound(false),upper_bound(0.0) {}
  interval(double l, double u):has_lower_bound(true),lower_bound(l),has_upper_bound(true),upper_bound(u) {}
  interval(bool hl,double l, bool hu, double u):has_lower_bound(hl),lower_bound(l),has_upper_bound(hu),upper_bound(u) {}
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
        num_evals++;

        variable->getValue() = x;

        // first we touch all the nodes
        // that will set the flags for recomputation
        variable->touch();

        assert( not variable->isClamped() );

        double Pr_ = (*this)();

        // call accept for each node  --  automatically includes affected nodes
        variable->keep();

        return Pr_;
    }

    double current_value() const
    {
        return variable->getValue();
    }

    slice_function(StochasticNode<double> *n, double pr, double l, double p, bool pos_only=false)
        :variable(n),
          lHeat(l),
          pHeat(p),
          prHeat(pr),
          num_evals(0)
    {
        variable->initiateGetAffectedNodes( affectedNodes );

        has_lower_bound = pos_only;
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

  if (g.below_lower_bound(L)) L = g.lower_bound;
  if (g.above_upper_bound(R)) R = g.upper_bound;

  assert(L < R);

  return std::pair<double,double>(L,R);
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

double slice_sample(double x0, slice_function& g,double w, int m)
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

  std::pair<double,double> interval = find_slice_boundaries_stepping_out(x0,g,logy,w,m);
  double L = interval.first;
  double R = interval.second;

  // Sample from the interval, shrinking it on each rejection

  return search_interval(x0,L,R,g,logy);
}


void SliceSamplingMove::performMcmcMove( double prHeat, double lHeat, double pHeat )
{
  slice_function g(variable, prHeat, lHeat, pHeat);

  double x1 = g.current_value();

  double x2 = slice_sample(x1, g, window, 100);

  total_movement += std::abs(x2 - x1);

  numPr += g.get_num_evals();

  if (auto_tuning and num_tried_total > 3)
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

