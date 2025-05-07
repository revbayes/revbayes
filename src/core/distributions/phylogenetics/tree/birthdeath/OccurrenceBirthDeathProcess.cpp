#include <cfloat>
#include <cstddef>
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
                                                                                           const TypedDagNode< RbVector<double> > *speciation_timeline,
                                                                                           const TypedDagNode< RbVector<double> > *extinction_timeline,
                                                                                           const TypedDagNode< RbVector<double> > *sampling_timeline,
                                                                                           const TypedDagNode< RbVector<double> > *treatment_timeline,
                                                                                           const TypedDagNode< RbVector<double> > *occurrence_timeline,
                                                                                           const std::string &cdt,
                                                                                           const std::vector<Taxon> &tn,
                                                                                           bool uo,
                                                                                           Tree *t,
                                                                                           const TypedDagNode<long> *n,
                                                                                           const TypedDagNode<RbVector<double> > *O,
                                                                                           bool mt,
                                                                                           bool vb) : AbstractBirthDeathProcess( sa, cdt, tn, uo, t ),
    interval_times_global(ht),
    offset( 0.0 ),
    start_age(sa),
    maxHiddenLin(n),
    occurrence_ages(O),
    useMt ( mt ),
    cond (cdt),
    verbose ( vb ),
    interval_times_speciation(speciation_timeline),
    interval_times_extinction(extinction_timeline),
    interval_times_sampling(sampling_timeline),
    interval_times_treatment(treatment_timeline),
    interval_times_occurrence(occurrence_timeline)

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

    homogeneous_omega    = NULL;
    heterogeneous_omega  = NULL;

    homogeneous_rho      = NULL;


    // We use a global timeline if
    //    1) the user provides one
    //    2) the user provides NO timeline arguments, meaning this is a constant-rate model
    if ( interval_times_global != NULL )
    {
      using_global_timeline = true;

      std::vector<double> times = interval_times_global->getValue();
      std::vector<double> times_sorted_ascending = times;

      sort(times_sorted_ascending.begin(), times_sorted_ascending.end() );

      if ( times != times_sorted_ascending )
      {
          throw RbException("Rate change times must be provided in ascending order.");
      }

    }

    addParameter( interval_times_global );

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

    heterogeneous_omega = dynamic_cast<const TypedDagNode<RbVector<double> >*>(inoccurrence);
    homogeneous_omega = dynamic_cast<const TypedDagNode<double >*>(inoccurrence);

    addParameter( homogeneous_omega );
    addParameter( heterogeneous_omega );

    homogeneous_rho = dynamic_cast<const TypedDagNode<double >*>(ineventsampling);
    addParameter( homogeneous_rho );

    prepareTimeline();

    if (t == NULL) {
          redrawValue();
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

    double lnProbTimes = computeLnProbabilityTimes();
    return lnProbTimes;
 }


/**
 * Compute the log probability of the current value under the current parameter values.
 */
double OccurrenceBirthDeathProcess::computeLnProbabilityTimes( void ) const
{
    prepareTimeline();

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

    double logLikelihood = ComputeLnLikelihoodOBDP(start_age->getValue(), global_timeline, lambda, mu, psi, omega, homogeneous_rho, r, maxHiddenLin, cond, useMt, verbose, occAges, tree);

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
    // @TODO @efficiency: this would be much faster if we can get std::lower_bound to work consistently
    // Linear search for interval because std::lower_bound is not cooperating
    if (global_timeline.size() == 1)
    {
        // If global_timeline.size() is 1, we have 0 break points and are in constant-rate version
        return 0;
    }
    else if (global_timeline.size() == 2)
    {
        return (t < (global_timeline[1]-1E-5) ? 0 : 1);
    }
    else
    {
        for (size_t i=0; i < global_timeline.size()-1; ++i)
        {
            if (t >= (global_timeline[i]-1E-5) && t < (global_timeline[i+1]-1E-5))
            {
                return i;
            }
        }

        return global_timeline.size() - 1;
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


/*
 * Here wepopulate all parameter vectors with their final values.
 * This requires that we:
 *    1) Clear out old values of all parameter vectors
 *    2) Refill and sort vector-valued parameters (leaving scalar parameters alone) to go from present to past
 *    3) Sort (assemble first if needed) the global timeline, attach the first time (the offset)
 * Then we can fill in our final vector for each parameter, which will be a vector of the same size as the global timeline
 */
void OccurrenceBirthDeathProcess::prepareTimeline( void ) const
{
    // clean all the sets
    lambda.clear();
    mu.clear();
    psi.clear();
    r.clear();
    omega.clear();
    psi_event.clear();

    lambda_times.clear();
    mu_times.clear();
    psi_times.clear();
    r_times.clear();
    omega_times.clear();
    global_timeline.clear();

    // put in current values for vector parameters so we can re-order them as needed
    if (heterogeneous_lambda != NULL)
    {
        lambda = heterogeneous_lambda->getValue();
    }
    if (heterogeneous_mu != NULL)
    {
        mu = heterogeneous_mu->getValue();
    }
    if (heterogeneous_psi != NULL)
    {
        psi = heterogeneous_psi->getValue();
    }
    if (heterogeneous_r != NULL)
    {
        r = heterogeneous_r->getValue();
    }
    if (heterogeneous_omega != NULL)
    {
        omega = heterogeneous_omega->getValue();
    }

    //@TODO we need to check that we have either a scalar or a vector for ALL of lambda/mu/phi/r/Phi (Lambda and Mu are allowed to be NULL), this should probably be done here
    // put in current values for vector parameters so we can re-order them as needed
    if (interval_times_global != NULL)
    {
        global_timeline = interval_times_global->getValue();
    }
    if (interval_times_speciation != NULL)
    {
        lambda_times = interval_times_speciation->getValue();
    }
    if (interval_times_extinction != NULL)
    {
        mu_times = interval_times_extinction->getValue();
    }
    if (interval_times_sampling != NULL)
    {
        psi_times = interval_times_sampling->getValue();
    }
    if (interval_times_treatment != NULL)
    {
        r_times = interval_times_treatment->getValue();
    }
    if (interval_times_occurrence != NULL)
    {
        omega_times = interval_times_occurrence->getValue();
    }


    // If it's a constant-rate process, make sure we only have scalars
    bool using_constant_rate_process = isConstantRate();
    if ( using_constant_rate_process )
    {
        global_timeline = std::vector<double>(0,0.0);
    }
    // If we have a real
    else if ( using_global_timeline )
    {
        // std::cout << "using global timeline" << std::endl;
        if ( interval_times_speciation != NULL ||
             interval_times_extinction != NULL ||
             interval_times_sampling != NULL ||
             interval_times_treatment != NULL ||
             interval_times_occurrence != NULL)
        {
            throw RbException("Both heterogeneous and homogeneous rate change times provided");
        }

        // check that the number of provided parameters matches the global timeline
        // Right now, the global timeline is only the interval times, i.e. breaks between pieces/episodes/windows
        checkVectorSizes(heterogeneous_lambda,interval_times_global,1,spn,true);
        checkVectorSizes(heterogeneous_mu,interval_times_global,1,exn,true);
        checkVectorSizes(heterogeneous_psi,interval_times_global,1,smp,true);
        checkVectorSizes(heterogeneous_r,interval_times_global,1,trt,false);
        checkVectorSizes(heterogeneous_omega,interval_times_global,1,occ,false);

        sortGlobalTimesAndVectorParameter();

        // we are done with setting up the timeline (i.e., using the provided global timeline) and checking all dimension of parameters
    }
    // We only need to assemble a global timeline if
    else if ( interval_times_speciation != NULL ||
              interval_times_extinction != NULL ||
              interval_times_sampling != NULL ||
              interval_times_treatment != NULL ||
              interval_times_occurrence != NULL  )
    {
        // check if correct number of speciation rates were provided
        // if provided as a vector, sort to the correct timescale
        if ( heterogeneous_lambda == NULL && homogeneous_lambda == NULL)
        {
            throw RbException("Speciation rate must be of type RealPos or RealPos[]");
        }
        else if ( heterogeneous_lambda != NULL )
        {
            if ( interval_times_speciation == NULL ) throw RbException("No time intervals provided for piecewise constant speciation rates");
            checkVectorSizes(heterogeneous_lambda,interval_times_speciation,1,spn,true);
            sortNonGlobalTimesAndVectorParameter(lambda_times,lambda);
        }

        // check if correct number of extinction rates were provided
        // if provided as a vector, sort to the correct timescale
        if ( heterogeneous_mu == NULL && homogeneous_mu == NULL)
        {
            throw RbException("Extinction rate must be of type RealPos or RealPos[]");
        }
        else if ( heterogeneous_mu != NULL )
        {
            if ( interval_times_extinction == NULL ) throw RbException("No time intervals provided for piecewise constant extinction rates");
            checkVectorSizes(heterogeneous_mu,interval_times_extinction,1,exn,true);
            sortNonGlobalTimesAndVectorParameter(mu_times,mu);
        }

        // check if correct number of sampling rates were provided
        // if provided as a vector, sort to the correct timescale
        if ( heterogeneous_psi == NULL && homogeneous_psi == NULL)
        {
            throw RbException("Sampling rate must be of type RealPos or RealPos[]");
        }
        else if ( heterogeneous_psi != NULL )
        {
            if ( interval_times_sampling == NULL ) throw RbException("No time intervals provided for piecewise constant sampling rates");
            checkVectorSizes(heterogeneous_psi,interval_times_sampling,1,smp,true);
            sortNonGlobalTimesAndVectorParameter(psi_times,psi);
        }

        // check if correct number of serial treatment probabilities were provided
        // if provided as a vector, sort to the correct timescale
        if ( heterogeneous_r == NULL && homogeneous_r == NULL)
        {
            throw RbException("Treatment probability for serial sampling rate must be of type Probability or Probability[]");
        }
        else if ( heterogeneous_r != NULL )
        {
            if ( interval_times_treatment == NULL ) throw RbException("No time intervals provided for piecewise constant sampling rates");
            checkVectorSizes(heterogeneous_r,interval_times_treatment,1,trt,false);
            sortNonGlobalTimesAndVectorParameter(r_times,r);
        }

        // check if correct number of occurrence sampling rates were provided
        // if provided as a vector, sort to the correct timescale
        if ( heterogeneous_omega == NULL && homogeneous_omega == NULL)
        {
            throw RbException("Occurrence sampling rate must be of type Probability or Probability[]");
        }
        else if ( heterogeneous_omega != NULL )
        {
            if ( interval_times_occurrence == NULL ) throw RbException("No time intervals provided for piecewise constant occurrence sampling rates");
            checkVectorSizes(heterogeneous_omega,interval_times_occurrence,1,occ,false);
            sortNonGlobalTimesAndVectorParameter(omega_times,omega);
        }

        // check that extant sampling probability is provided
        if ( homogeneous_rho == NULL)
        {
          throw RbException("Event-sampling probabilities must be of type Probability");
        }


        // now we start assembling the global timeline by finding the union of unique intervals for all parameters
        std::set<double> event_times;
        addTimesToGlobalTimeline(event_times,interval_times_speciation);
        addTimesToGlobalTimeline(event_times,interval_times_extinction);
        addTimesToGlobalTimeline(event_times,interval_times_sampling);
        addTimesToGlobalTimeline(event_times,interval_times_treatment);
        addTimesToGlobalTimeline(event_times,interval_times_occurrence);

        // we are done with setting up the timeline (i.e., using the all the provided timeline) and checking all dimension of parameters

    }

    // @TODO: @ANDY: Double check the offset works
    // Add s_0
    getOffset();
    global_timeline.insert(global_timeline.begin(),offset);

    // For each parameter vector, we now make sure that its size matches the size of the global vector
    // For a RATE parameter, there are three cases
    //     1) It is a vector and it matches the size of the global timeline, in which case it is already sorted and we can use it
    //     2) It is a vector and it DOES NOT match the size of the global timeline, in which case we must expand it to match
    //     3) It is a scalar, in which case we simply populate a vector of the correct size with the value

    // Get vector of birth rates
    // @TODO: @SEBASTIAN: would it be better here to check if interval_times_parameter == NULL instead of checking the size? They should be equivalent
    if ( heterogeneous_lambda != NULL )
    {
        if ( lambda.size() != global_timeline.size() )
        {
            expandNonGlobalRateParameterVector(lambda,lambda_times);
        } // else it matches in size and is already sorted and is thus ready to be used
    }
    else
    {
        lambda = std::vector<double>(global_timeline.size(),homogeneous_lambda->getValue());
    }

    // Get vector of death rates
    if ( heterogeneous_mu != NULL )
    {
        if ( mu.size() != global_timeline.size() )
        {
            expandNonGlobalRateParameterVector(mu,mu_times);
        } // else it matches in size and is already sorted and is thus ready to be used
    }
    else
    {
        mu = std::vector<double>(global_timeline.size(),homogeneous_mu->getValue());
    }

    // Get vector of sampling rates
    if ( heterogeneous_psi != NULL )
    {
        if ( psi.size() != global_timeline.size() )
        {
            expandNonGlobalRateParameterVector(psi,psi_times);
        } // else it matches in size and is already sorted and is thus ready to be used
    }
    else
    {
        psi = std::vector<double>(global_timeline.size(),homogeneous_psi->getValue());
    }

    // Get vector of treatment probabilities
    if ( heterogeneous_r != NULL )
    {
        if ( r.size() != global_timeline.size() )
        {
            // r is not a rate parameter, but it behaves like them for this function, as it is defined in intervals
            expandNonGlobalRateParameterVector(r,r_times);
        } // else it matches in size and is already sorted and is thus ready to be used
    }
    else
    {
        r = std::vector<double>(global_timeline.size(),homogeneous_r->getValue());
    }

    // Get vector of treatment probabilities
    if ( heterogeneous_omega != NULL )
    {
        if ( omega.size() != global_timeline.size() )
        {
            // r is not a rate parameter, but it behaves like them for this function, as it is defined in intervals
            expandNonGlobalRateParameterVector(omega,omega_times);
        } // else it matches in size and is already sorted and is thus ready to be used
    }
    else
    {
        omega = std::vector<double>(global_timeline.size(),homogeneous_omega->getValue());
    }
    psi_event = std::vector<double>(global_timeline.size(),0.0);
    if ( homogeneous_rho != NULL )
    {
        // User specified the sampling fraction at the present
        psi_event[0] = homogeneous_rho->getValue();
    }




}

/**
 * Adds parameter-specific timeline to the set
 */
void OccurrenceBirthDeathProcess::addTimesToGlobalTimeline(std::set<double> &event_times, const TypedDagNode<RbVector<double> > *par_times) const
{
  if ( par_times != NULL )
  {
      const std::vector<double>& times = par_times->getValue();
      for (std::vector<double>::const_iterator it = times.begin(); it != times.end(); ++it)
      {
          event_times.insert( *it );
      }
  }

}

bool OccurrenceBirthDeathProcess::isConstantRate(void) const
{
  bool has_no_interval_times = false;
  // For there to be no intervals, every timeline must either be NULL or have size 0
  if ( (interval_times_global == NULL           || interval_times_global->getValue().size() == 0 ) &&
       (interval_times_speciation == NULL       || interval_times_speciation->getValue().size() == 0 ) &&
       (interval_times_extinction == NULL       || interval_times_extinction->getValue().size() == 0 ) &&
       (interval_times_sampling == NULL         || interval_times_sampling->getValue().size() == 0 ) &&
       (interval_times_treatment == NULL        || interval_times_treatment->getValue().size() == 0 ) &&
       (interval_times_occurrence == NULL       || interval_times_occurrence->getValue().size() == 0 ))
  {
    has_no_interval_times = true;
  }

  bool all_parameters_are_scalars = false;
  // For all parameters to be scalars,
  // 1) rate parameters must either be homogenous or they must have size <= 1 (1 for scalar, 0 if it's null)
  if ( (heterogeneous_lambda == NULL || lambda.size() <= 1) &&
       (heterogeneous_mu == NULL     || mu.size() <= 1) &&
       (heterogeneous_psi == NULL    || psi.size() <= 1) &&
       (heterogeneous_r == NULL      || r.size() <= 1) &&
       (heterogeneous_omega == NULL  || omega.size() <= 1))
       {
         all_parameters_are_scalars = true;
       }

  // std::cout << "has_no_interval_times == " << has_no_interval_times << "; all_parameters_are_scalars == " << all_parameters_are_scalars << std::endl;

  if (has_no_interval_times && !all_parameters_are_scalars)
  {
    throw RbException("No timeline(s) was (were) provided but there are non-scalar parameters.");
  }

  return has_no_interval_times && all_parameters_are_scalars;
}


/**
 * return the index i so that x_{i-1} <= t < x_i
 * where x is one of the input vector timelines
 */
size_t OccurrenceBirthDeathProcess::findIndex(double t, const std::vector<double>& timeline) const
{

    // Linear search for interval because std::lower_bound is not cooperating
    if (timeline.size() == 1)
    {
        return 0;
    }
    else if (timeline.size() == 2)
    {
        return(t < timeline[1] ? 0 : 1);
    }
    else
    {
        for (size_t i=0; i < timeline.size()-1; ++i)
        {
            if (t >= timeline[i] && t < timeline[i+1])
            {
                return i;
            }
        }

        return timeline.size() - 1;
    }
}

/**
 * Checks that v1 is the correct size compared to reference vector ref, given the expected size difference.
 * If the first vector parameter is not a vector (is null), does nothing.
 * If the sizes are wrong, throws an exception.
 * Uses param_name and is_rate to make a sensible error message
 */
void OccurrenceBirthDeathProcess::checkVectorSizes(const TypedDagNode<RbVector<double> >* v, const TypedDagNode<RbVector<double> >* ref, int v1_minus_ref, const std::string& param_name, bool is_rate) const
{
  if ( v != NULL )
  {
    if ( v->getValue().size() - ref->getValue().size() != v1_minus_ref )
    {
      std::string vec_type = is_rate ? "rates" : "probabilities";
      std::stringstream ss;
      ss << "Number of " << param_name << " " << vec_type << " (" << v->getValue().size() << ") does not match number of time intervals (" << ref->getValue().size() << ")";
      throw RbException(ss.str());
    }
  }
}

/**
 * Sorts times to run from present to past (0->inf) and orders par to match this.
 */
void OccurrenceBirthDeathProcess::sortNonGlobalTimesAndVectorParameter(std::vector<double> &times, std::vector<double>& par) const
{
  std::vector<double> times_sorted_ascending = times;
  std::vector<double> times_sorted_descending = times;

  sort(times_sorted_ascending.begin(), times_sorted_ascending.end() );
  sort(times_sorted_descending.rbegin(), times_sorted_descending.rend() );

  // We want times in ascending order, so if they already are we're done here
  if ( times != times_sorted_ascending )
  {
      // If times are sorted in descending order, we just flip the parameter and time vectors
      if ( times == times_sorted_ascending )
      {
        std::reverse(times.begin(),times.end());
        std::reverse(par.begin(),par.end());
      }
      else
      {
        // Pair up the times and the parameter values so we can sort them together
        std::vector<std::pair<double,double> > times_par;
        for (size_t i=0; i<times.size(); ++i)
        {
          times_par.push_back(std::make_pair(times[i],par[i]));
        }

        std::sort(times_par.begin(),times_par.end());

        // Replace times with sorted times
        for (size_t i=0; i<times.size(); ++i)
        {
          times[i] = times_par[i].first;
          par[i] = times_par[i].second;
        }
      }
  }

  if ( times[0] < DBL_EPSILON )
  {
    throw RbException("User-specified interval times cannot include time = 0");
  }

}



/**
 * Takes a par.size() < global_timeline.size() vector and makes it the correct size to work with our global timeline.
 * The parameter has its own reference timeline, which we use to find the rate in the global intervals.
 * This works only for parameters (lambda,mu,phi,r), where the global timeline is simply a finer grid than the variable-specific timelines.
 */
void OccurrenceBirthDeathProcess::expandNonGlobalRateParameterVector(std::vector<double> &par, const std::vector<double> &par_times) const
{
    // Store the original values so we can overwrite the vector
    std::vector<double> old_par = par;

    // For each time in the global timeline, find the rate according to this variable's own timeline
    for (size_t i=0; i<global_timeline.size(); ++i)
    {
      // Where is this global time interval in the variable's timeline?
      size_t idx = findIndex(global_timeline[i],par_times);
      par[i] = old_par[idx];
    }

}

/**
 * Sorts global times to run from present to past (0->inf) and orders ALL vector parameters to match this.
 * These can only be sorted after the local copies have values in them.
 */
void OccurrenceBirthDeathProcess::sortGlobalTimesAndVectorParameter( void ) const
{
    std::vector<double> times_sorted_ascending  = global_timeline;
    std::vector<double> times_sorted_descending = global_timeline;

    sort(times_sorted_ascending.begin(), times_sorted_ascending.end() );
    sort(times_sorted_descending.rbegin(), times_sorted_descending.rend() );

    // We want times in ascending order, so if they already are we're done here
    if ( global_timeline != times_sorted_ascending )
    {
        // If times are sorted in descending order, we just flip the parameter and time vectors
        if ( global_timeline == times_sorted_descending )
        {
            // Reverse timeline
            std::reverse(global_timeline.begin(),global_timeline.end());

            // @TODO: @ANDY: These checks for NULL might be superfluous because the std vectors are initialized to empty vectors by default in c++
            // Reverse all vector parameters
            if (heterogeneous_lambda != NULL)
            {
                sort(lambda.rbegin(),lambda.rend());
            }
            if (heterogeneous_mu != NULL)
            {
                sort(mu.rbegin(),mu.rend());
            }
            if (heterogeneous_psi != NULL)
            {
                sort(psi.rbegin(),psi.rend());
            }
            if (heterogeneous_r != NULL)
            {
                sort(r.rbegin(),r.rend());
            }
            if (heterogeneous_omega != NULL)
            {
                sort(omega.rbegin(),omega.rend());
            }

        }
        else
        {
            // Find ordering of times vector
            std::vector<size_t> ordering;
            for (size_t i=0; i<global_timeline.size(); ++i)
            {
                for (size_t j=0; j<global_timeline.size(); ++j)
                {
                    if ( times_sorted_ascending[i] == global_timeline[j] )
                    {
                        ordering.push_back(j);
                        break;
                    }
                }
            }

            // Replace times with sorted times
            global_timeline = times_sorted_ascending;

            // Sort all vector parameters
            if (heterogeneous_lambda != NULL)
            {
                std::vector<double> old_lambda = lambda;
                for (size_t i=0; i<global_timeline.size(); ++i)
                {
                    lambda[i] = old_lambda[ordering[i]];
                }
            }
            if (heterogeneous_mu != NULL)
            {
                std::vector<double> old_mu = mu;
                for (size_t i=0; i<global_timeline.size(); ++i)
                {
                    mu[i] = old_mu[ordering[i]];
                }
            }
            if (heterogeneous_psi != NULL)
            {
                std::vector<double> old_psi = psi;
                for (size_t i=0; i<global_timeline.size(); ++i)
                {
                    psi[i] = old_psi[ordering[i]];
                }
            }
            if (heterogeneous_r != NULL)
            {
                std::vector<double> old_r = r;
                for (size_t i=0; i<global_timeline.size(); ++i)
                {
                    r[i] = old_r[ordering[i]];
                }
            }
            if (heterogeneous_omega != NULL)
            {
                std::vector<double> old_omega = omega;
                for (size_t i=0; i<global_timeline.size(); ++i)
                {
                    omega[i] = old_omega[ordering[i]];
                }
            }

        }
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
    //We use a coalescent simulator

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
    else if (oldP == heterogeneous_omega)
    {
        heterogeneous_omega      = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == homogeneous_omega)
    {
        homogeneous_omega        = static_cast<const TypedDagNode<double>* >( newP );
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
