#include <float.h>
#include <stddef.h>
#include <algorithm>
#include <cmath>
#include <iosfwd>
#include <string>
#include <vector>

#include "DistributionExponential.h"
#include "OccurrenceBirthDeathProcess.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbMathCombinatorialFunctions.h"
#include "RbMathLogic.h"
#include "AbstractBirthDeathProcess.h"
#include "DagNode.h"
#include "RbException.h"
#include "RbVector.h"
#include "TopologyNode.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "StartingTreeSimulator.h"


#include "ComputeLikelihoodsLtMt.h"

namespace RevBayesCore { class Taxon; }

using namespace RevBayesCore;

/**
 * Constructor.
 * We delegate most parameters to the base class and initialize the members.
 * \param[in]    sa                        Start age of the process.
 * \param[in]    inspeciation              Speciation/birth rate(s).
 * \param[in]    inextinction              Extinction/death rate(s).
 * \param[in]    inserialsampling          Serial sampling rate(s).
 * \param[in]    intreatment               Probabilit(y|ies) of death upon sampling (treatment).
 * \param[in]    inoccurrence              Occurrence sampling rate(s).
 * \param[in]    ineventsampling           Sampling probability at present time.
 * \param[in]    ht                        Rate interval change times of the piecewise constant process.
 * \param[in]    cdt                       Condition of the process (survival/survival2).
 * \param[in]    tn                        Taxa.
 * \param[in]    uo                        If true the start age is the origin time otherwise the root age of the process.
 * \param[in]    t                         Starting tree if we want to avoid simulating trees.
 * \param[in]    n                         Maximum number of hidden lineages (algorithm accuracy).
 * \param[in]    O                         Vector of occurrence ages.

 */
OccurrenceBirthDeathProcess::OccurrenceBirthDeathProcess(                                  const TypedDagNode<double> *sa,
                                                                                           const DagNode *inspeciation,
                                                                                           const DagNode *inextinction,
                                                                                           const DagNode *inserialsampling,
                                                                                           const DagNode *intreatment,
                                                                                           const DagNode *inoccurrence,
                                                                                           const DagNode *ineventsampling,
                                                                                           const TypedDagNode< RbVector<double> > *ht,
                                                                                           const std::string &cdt,
                                                                                           const std::vector<Taxon> &tn,
                                                                                           bool uo,
                                                                                           Tree *t,
                                                                                           const TypedDagNode<long> *n,
                                                                                           const TypedDagNode<RbVector<double> > *O,
                                                                                           bool mt,
                                                                                           bool vb) : AbstractBirthDeathProcess( sa, cdt, tn, uo, t ),
    interval_times(ht),
    offset( 0.0 ),
    start_age(sa),
    maxHiddenLin(n),
    occurrence_ages(O),
    useMt ( mt ),
    cond (cdt),
    verbose ( vb )

{
    // initialize all the pointers to NULL
    homogeneous_lambda   = NULL;
    heterogeneous_lambda = NULL;

    homogeneous_mu       = NULL;
    heterogeneous_mu     = NULL;

    homogeneous_psi      = NULL;
    heterogeneous_psi    = NULL;

    homogeneous_r        = NULL;
    heterogeneous_r      = NULL;

    homogeneous_o        = NULL;
    heterogeneous_o      = NULL;

    homogeneous_rho      = NULL;


    std::vector<double> times = timeline;
    std::vector<double> times_sorted_ascending = times;

    sort(times_sorted_ascending.begin(), times_sorted_ascending.end() );

    if ( times != times_sorted_ascending )
    {
        throw(RbException("Rate change times must be provided in ascending order."));
    }

    addParameter( interval_times );

    heterogeneous_lambda = dynamic_cast<const TypedDagNode<RbVector<double> >*>(inspeciation);
    homogeneous_lambda = dynamic_cast<const TypedDagNode<double >*>(inspeciation);

    addParameter( homogeneous_lambda );
    addParameter( heterogeneous_lambda );

    heterogeneous_mu = dynamic_cast<const TypedDagNode<RbVector<double> >*>(inextinction);
    homogeneous_mu = dynamic_cast<const TypedDagNode<double >*>(inextinction);

    addParameter( homogeneous_mu );
    addParameter( heterogeneous_mu );

    heterogeneous_psi = dynamic_cast<const TypedDagNode<RbVector<double> >*>(inserialsampling);
    homogeneous_psi = dynamic_cast<const TypedDagNode<double >*>(inserialsampling);

    addParameter( homogeneous_psi );
    addParameter( heterogeneous_psi );

    heterogeneous_r = dynamic_cast<const TypedDagNode<RbVector<double> >*>(intreatment);
    homogeneous_r = dynamic_cast<const TypedDagNode<double >*>(intreatment);

    addParameter( homogeneous_r );
    addParameter( heterogeneous_r );

    heterogeneous_o = dynamic_cast<const TypedDagNode<RbVector<double> >*>(inoccurrence);
    homogeneous_o = dynamic_cast<const TypedDagNode<double >*>(inoccurrence);

    addParameter( homogeneous_o );
    addParameter( heterogeneous_o );

    homogeneous_rho = dynamic_cast<const TypedDagNode<double >*>(ineventsampling);
    addParameter( homogeneous_rho );


    //check that lengths of vector arguments are sane
    if ( heterogeneous_lambda != NULL && !(interval_times->getValue().size() == heterogeneous_lambda->getValue().size() - 1) )
    {
      throw(RbException("If provided as a vector, argument lambda must have one more element than timeline."));
    }

    if ( heterogeneous_mu != NULL && !(interval_times->getValue().size() == heterogeneous_mu->getValue().size() - 1) )
    {
      throw(RbException("If provided as a vector, argument mu must have one more element than timeline."));
    }

    if ( heterogeneous_psi != NULL && !(interval_times->getValue().size() == heterogeneous_psi->getValue().size() - 1) )
    {
      throw(RbException("If provided as a vector, argument psi must have one more element than timeline."));
    }

    if ( heterogeneous_r != NULL && !(interval_times->getValue().size() == heterogeneous_r->getValue().size() - 1) )
    {
      throw(RbException("If provided as a vector, argument r must have one more element than timeline."));
    }

    if ( heterogeneous_o != NULL && !(interval_times->getValue().size() == heterogeneous_o->getValue().size() - 1) )
    {
      throw(RbException("If provided as a vector, argument o must have one more element than timeline."));
    }

    updateVectorParameters();



    if (t == NULL) {

        // with piecewise constant parameters, redrawValue() sometimes fails to initialise, we use instead the following:
        if(timeline.size() > 1) {

            RbVector<Clade> constr;
            StartingTreeSimulator simulator;
            Tree *my_tree = simulator.simulateTree( taxa, constr );

            // store the new value
            delete value;
            value = my_tree;
        }

        else {
            redrawValue();
        }
        
      }
}





/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
OccurrenceBirthDeathProcess* OccurrenceBirthDeathProcess::clone( void ) const
{
    return new OccurrenceBirthDeathProcess( *this );
}





/**
 * Compute the log-transformed probability of the current value under the current parameter values.
 */
 double OccurrenceBirthDeathProcess::computeLnProbabilityDivergenceTimes( void )
 {

     // variable declarations and initialization
    double lnProbTimes = computeLnProbabilityTimes();

    return lnProbTimes;
 }





/**
 * Compute the log probability of the current value under the current parameter values.
 */
double OccurrenceBirthDeathProcess::computeLnProbabilityTimes( void ) const
{
    updateVectorParameters();

    // compute the log-likelihood : use ComputeLikelihoodsBackwardsLt (backward traversal of the tree) or ComputeLikelihoodsForwardsMt (forward traversal of the tree)
    const Tree& tree = *value;

    std::vector<double> occAges;
    if (occurrence_ages != NULL)
    {
        occAges = occurrence_ages->getValue();
    }
    else
    {
        occAges = std::vector<double>();
    }

    double logLikelihood = ComputeLnLikelihoodOBDP(start_age->getValue(), timeline, lambda, mu, psi, omega, homogeneous_rho, r, maxHiddenLin, cond, useMt, verbose, occAges, tree);

    return logLikelihood;
}





/**
 * return the index i so that t_{i-1} <= t < t_i
 * where t_i is the instantaneous sampling time (i = 0,...,l)
 * t_0 = 0.0
 * t_l = Inf
 */
size_t OccurrenceBirthDeathProcess::findIndex(double t) const
{
    // Linear search for interval because std::lower_bound is not cooperating
    if (timeline.size() == 2) // If timeline.size() were 1, we would have 0 break points and would be in constant-rate version
    {
      return(t < timeline[1] ? 0 : 1);
    }
    else
    {
      for (size_t i=0; i < timeline.size()-1; ++i)
      {
          if (t >= timeline[i] && t < timeline[i+1])
          {
              return(i);
          }
      }

      return(timeline.size() - 1);
    }
}





// calculate offset so we can set t_0 to time of most recent tip
void OccurrenceBirthDeathProcess::getOffset(void) const
{
    offset = RbConstants::Double::max;
    for (size_t i = 0; i < value->getNumberOfNodes(); i++)
    {
        const TopologyNode& n = value->getNode( i );

        if ( n.getAge() < offset )
        {
            offset = n.getAge();
        }
    }
}





/**
 * Here we ensure that we have a length-l vector of rates and events for all parameters and the timeline.
 * In the case of homogeneous/constant-rate parameters, we fill the vector with this rate.
 * In the case of homogenous/single events, we make all but the first event rates 0.
 *
 */
void OccurrenceBirthDeathProcess::updateVectorParameters( void ) const
{
    // clean and get timeline
    timeline.clear();
    timeline = interval_times->getValue();

    // Add t_0
    getOffset();
    timeline.insert(timeline.begin(),offset);

    // clean all the sets
    lambda.clear();
    mu.clear();
    psi.clear();
    r.clear();
    omega.clear();
    psi_event.clear();

    // Get vector of birth rates
    if ( heterogeneous_lambda != NULL )
    {
      lambda = heterogeneous_lambda->getValue();
    }
    else
    {
      lambda = std::vector<double>(timeline.size(),homogeneous_lambda->getValue());
    }

    // Get vector of death rates
    if ( heterogeneous_mu != NULL )
    {
      mu = heterogeneous_mu->getValue();
    }
    else
    {
      mu = std::vector<double>(timeline.size(),homogeneous_mu->getValue());
    }

    // Get vector of serial sampling rates
    if ( heterogeneous_psi != NULL )
    {
      psi = heterogeneous_psi->getValue();
    }
    else
    {
      psi = std::vector<double>(timeline.size(),homogeneous_psi->getValue());
    }

    // Get vector of conditional death upon sampling probabilities
    if ( heterogeneous_r != NULL )
    {
      r = heterogeneous_r->getValue();
    }
    else
    {
      r = std::vector<double>(timeline.size(),homogeneous_r->getValue());
    }

    // Get vector of occcurrence rates.
    if ( heterogeneous_o != NULL )
    {
      omega = heterogeneous_o->getValue();
    }
    else
    {
      omega = std::vector<double>(timeline.size(),homogeneous_o->getValue());
    }

    psi_event = std::vector<double>(timeline.size(),0.0);
    if ( homogeneous_rho != NULL )
    {
        // User specified the sampling fraction at the present
        psi_event[0] = homogeneous_rho->getValue();
    }
    else
    {
        // set the final sampling to one (for sampling at the present)
        psi_event[0] = 1.0;
    }
}





double OccurrenceBirthDeathProcess::pSurvival(double start, double end) const
{
    if(end != 0.0) {
        //To do: this only makes sense if the rate vectors match the timeline.
        //the vector parameters and timeline have to be adjusted.
        double time = start - end;
        std::vector<double> birth;
        std::vector<double> death;
        std::vector<double> ps;
        std::vector<double> om;
        std::vector<double> rp;
        std::vector<double> d;

    for(int i = 0; i < d.size(); i++) {

        if ((timeline[i] - end) > 0)  {
            d.push_back(timeline[i]);
            birth.push_back(lambda[i]);
            death.push_back(mu[i]);
            ps.push_back(psi[i]);
            om.push_back(omega[i]);
            rp.push_back(r[i]);
        }
    }
        std::vector<double> res = GetFunctionUandP(time, d, birth, death, ps, om, homogeneous_rho, rp);
        return 1 - res[0] ;
  }
    else {
        std::vector<double> res = GetFunctionUandP(start, timeline, lambda, mu, psi, omega, homogeneous_rho, r);
        return 1 - res[0] ;
    }
}





/**
 * Simulate new speciation times.
 */
double OccurrenceBirthDeathProcess::simulateDivergenceTime(double origin, double present) const
{
    // incorrect placeholder for constant SSBDP


    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;

    size_t i = findIndex(present);

    // get the parameters
    double age = origin - present;
    double b = lambda[i];
    double d = mu[i];
    double p_e = psi_event[i];


    // get a random draw
    double u = rng->uniform01();

    // compute the time for this draw
    // see Hartmann et al. 2010 and Stadler 2011
    double t = 0.0;
    if ( b > d )
    {
        if( p_e > 0.0 )
        {
            t = ( log( ( (b-d) / (1 - (u)*(1-((b-d)*exp((d-b)*age))/(p_e*b+(b*(1-p_e)-d)*exp((d-b)*age) ) ) ) - (b*(1-p_e)-d) ) / (p_e * b) ) )  /  (b-d);
        }
        else
        {
            t = log( 1 - u * (exp(age*(d-b)) - 1) / exp(age*(d-b)) ) / (b-d);
        }
    }
    else
    {
        if( p_e > 0.0 )
        {
            t = ( log( ( (b-d) / (1 - (u)*(1-(b-d)/(p_e*b*exp((b-d)*age)+(b*(1-p_e)-d) ) ) ) - (b*(1-p_e)-d) ) / (p_e * b) ) )  /  (b-d);
        }
        else
        {
            t = log( 1 - u * (1 - exp(age*(b-d)))  ) / (b-d);
        }
    }

    return present + t;
}





/**
 * Swap the parameters held by this distribution.
 *
 * \param[in]    oldP      Pointer to the old parameter.
 * \param[in]    newP      Pointer to the new parameter.
 */
void OccurrenceBirthDeathProcess::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    // Rate parameters
    if (oldP == heterogeneous_lambda)
    {
        heterogeneous_lambda = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == heterogeneous_mu)
    {
        heterogeneous_mu     = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == heterogeneous_psi)
    {
        heterogeneous_psi    = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == homogeneous_lambda)
    {
        homogeneous_lambda   = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == homogeneous_mu)
    {
        homogeneous_mu       = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == homogeneous_psi)
    {
        homogeneous_psi      = static_cast<const TypedDagNode<double>* >( newP );
    }
    // Treatment
    else if (oldP == heterogeneous_r)
    {
        heterogeneous_r      = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == homogeneous_r)
    {
        homogeneous_r        = static_cast<const TypedDagNode<double>* >( newP );
    }
    //occurrence
    else if (oldP == heterogeneous_o)
    {
        heterogeneous_o      = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == homogeneous_o)
    {
        homogeneous_o        = static_cast<const TypedDagNode<double>* >( newP );
    }
    // Event probability parameters
    else if (oldP == homogeneous_rho)
    {
        homogeneous_rho      = static_cast<const TypedDagNode<double>* >( newP );
    }
    else
    {
        // delegate the super-class
        AbstractBirthDeathProcess::swapParameterInternal(oldP, newP);
    }
}
