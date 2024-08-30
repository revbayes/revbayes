#include <cfloat>
#include <cstddef>
#include <algorithm>
#include <cmath>
#include <iosfwd>
#include <string>
#include <vector>
#include <cassert>

#include "AbstractBirthDeathProcess.h"
#include "BirthDeathForwardSimulator.h"
#include "BirthDeathSamplingTreatmentProcess.h"
#include "Clade.h"
#include "DagNode.h"
#include "DistributionExponential.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbMathCombinatorialFunctions.h"
#include "RbMathLogic.h"
#include "RbException.h"
#include "RbSettings.h"
#include "RbVector.h"
#include "StartingTreeSimulator.h"
#include "TopologyNode.h"
#include "Tree.h"
#include "TreeUtilities.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class Taxon; }

using namespace RevBayesCore;

/**
 * Constructor.
 * We delegate most parameters to the base class and initialize the members.
 *
 * \param[in]    s              Speciation rates.
 * \param[in]    e              Extinction rates.
 * \param[in]    p              Serial sampling rates.
 * \param[in]    r              Instantaneous sampling probabilities.
 * \param[in]    t              Rate change times.
 * \param[in]    cdt            Condition of the process (none/survival/#Taxa).
 * \param[in]    tn             Taxa.
 * \param[in]    c              Clades conditioned to be present.
 */
BirthDeathSamplingTreatmentProcess::BirthDeathSamplingTreatmentProcess(const TypedDagNode<double> *ra,
                                                                        const DagNode *in_speciation,
                                                                        const DagNode *in_extinction,
                                                                        const DagNode *in_sampling,
                                                                        const DagNode *in_treatment,
                                                                        const DagNode *in_event_speciation,
                                                                        const DagNode *in_event_extinction,
                                                                        const DagNode *in_event_sampling,
                                                                        const DagNode *in_event_treatment,
                                                                        const TypedDagNode< RbVector<double> > *timeline,
                                                                        const TypedDagNode< RbVector<double> > *speciation_timeline,
                                                                        const TypedDagNode< RbVector<double> > *extinction_timeline,
                                                                        const TypedDagNode< RbVector<double> > *sampling_timeline,
                                                                        const TypedDagNode< RbVector<double> > *treatment_timeline,
                                                                        const TypedDagNode< RbVector<double> > *event_speciation_timeline,
                                                                        const TypedDagNode< RbVector<double> > *event_extinction_timeline,
                                                                        const TypedDagNode< RbVector<double> > *event_sampling_timeline,
                                                                        const std::string &cdt,
                                                                        const std::vector<Taxon> &tn,
                                                                        bool uo,
                                                                        Tree *t,
                                                                        long age_check_precision) : AbstractBirthDeathProcess( ra, cdt, tn, uo, t ),
    interval_times_global(timeline),
    interval_times_speciation(speciation_timeline),
    interval_times_extinction(extinction_timeline),
    interval_times_sampling(sampling_timeline),
    interval_times_treatment(treatment_timeline),
    interval_times_event_speciation(event_sampling_timeline),
    interval_times_event_extinction(event_extinction_timeline),
    interval_times_event_sampling(event_sampling_timeline),
    offset( 0.0 ),
    taxa(tn),
    age_check_precision(age_check_precision)
{
    // initialize all the pointers to NULL
    homogeneous_lambda   = NULL;
    homogeneous_mu       = NULL;
    homogeneous_phi      = NULL;
    homogeneous_r        = NULL;
    // homogeneous_Lambda   = NULL;
    // homogeneous_Mu       = NULL;
    homogeneous_Phi      = NULL;
    heterogeneous_lambda = NULL;
    heterogeneous_mu     = NULL;
    heterogeneous_phi    = NULL;
    heterogeneous_r      = NULL;
    heterogeneous_Lambda = NULL;
    heterogeneous_Mu     = NULL;
    heterogeneous_Phi    = NULL;
    heterogeneous_R      = NULL;

    //@TODO @SEBASTIAN: sometime we might want to allow "homogeneous" aka scalar Mu/Lambda

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

    heterogeneous_lambda = dynamic_cast<const TypedDagNode<RbVector<double> >*>(in_speciation);
    homogeneous_lambda   = dynamic_cast<const TypedDagNode<double >*>(in_speciation);

    addParameter( homogeneous_lambda );
    addParameter( heterogeneous_lambda );

    heterogeneous_mu = dynamic_cast<const TypedDagNode<RbVector<double> >*>(in_extinction);
    homogeneous_mu   = dynamic_cast<const TypedDagNode<double >*>(in_extinction);

    addParameter( homogeneous_mu );
    addParameter( heterogeneous_mu );

    heterogeneous_phi = dynamic_cast<const TypedDagNode<RbVector<double> >*>(in_sampling);
    homogeneous_phi   = dynamic_cast<const TypedDagNode<double >*>(in_sampling);

    addParameter( homogeneous_phi );
    addParameter( heterogeneous_phi );

    heterogeneous_r = dynamic_cast<const TypedDagNode<RbVector<double> >*>(in_treatment);
    homogeneous_r   = dynamic_cast<const TypedDagNode<double >*>(in_treatment);

    addParameter( homogeneous_r );
    addParameter( heterogeneous_r );

    heterogeneous_Lambda = dynamic_cast<const TypedDagNode<RbVector<double> >*>(in_event_speciation);

    addParameter( heterogeneous_Lambda );

    heterogeneous_Mu = dynamic_cast<const TypedDagNode<RbVector<double> >*>(in_event_extinction);

    addParameter( heterogeneous_Mu );

    heterogeneous_Phi = dynamic_cast<const TypedDagNode<RbVector<double> >*>(in_event_sampling);
    homogeneous_Phi   = dynamic_cast<const TypedDagNode<double >*>(in_event_sampling);  

    addParameter( homogeneous_Phi );
    addParameter( heterogeneous_Phi );

    heterogeneous_R = dynamic_cast<const TypedDagNode<RbVector<double> >*>(in_event_treatment);

    addParameter( heterogeneous_R );

    //TODO: should check that the first interval time is not less than the first tip in offset computation

    //TODO: returning neginf and nan are not currently consistently used for different issues with invalid values.

    // updateVectorParameters();
    prepareTimeline();
    prepareProbComputation();
    
    delete value;

    if ( starting_tree == NULL )
    {
        RbVector<Clade> constr;
        // We employ a coalescent simulator to guarantee that the starting tree matches all time constraints
        StartingTreeSimulator simulator;
        RevBayesCore::Tree *my_tree = simulator.simulateTree( taxa, constr );
        // store the new value
        value = my_tree;
    }
    else
    {
        try
        {
            RevBayesCore::Tree *my_tree = TreeUtilities::startingTreeInitializer( *t, taxa, age_check_precision );
            value = my_tree->clone();
        }
        catch (RbException &e)
        {
            value = nullptr;
            // The line above is to prevent a segfault when ~AbstractRootedTreeDistribution() tries to delete
            // a nonexistent starting_tree
            throw RbException( e.getMessage() );
        }
    }

    countAllNodes();

}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
BirthDeathSamplingTreatmentProcess* BirthDeathSamplingTreatmentProcess::clone( void ) const
{
    return new BirthDeathSamplingTreatmentProcess( *this );
}

/**
 * Adds parameter-specific timeline to the set
 */
void BirthDeathSamplingTreatmentProcess::addTimesToGlobalTimeline(std::set<double> &event_times, const TypedDagNode<RbVector<double> > *par_times) const
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

/**
 * Checks that v1 is the correct size compared to reference vector ref, given the expected size difference.
 * If the first vector parameter is not a vector (is null), does nothing.
 * If the sizes are wrong, throws an exception.
 * Uses param_name and is_rate to make a sensible error message
 */
void BirthDeathSamplingTreatmentProcess::checkVectorSizes(const TypedDagNode<RbVector<double> >* v, const TypedDagNode<RbVector<double> >* ref, int v1_minus_ref, const std::string& param_name, bool is_rate) const
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
 * Compute the log-transformed probability of the current value under the current parameter values.
 *
 */
double BirthDeathSamplingTreatmentProcess::computeLnProbabilityDivergenceTimes( void ) const
{
    
    // update parameter vectors
    prepareTimeline();

    // Assign nodes to sets
    if ( countAllNodes() )
    {
        return RbConstants::Double::neginf;
    }

    if ( offset > DBL_EPSILON && phi_event[0] > DBL_EPSILON )
    {
        throw RbException("Event sampling fraction at the present is non-zero but there are no tips at the present.");
    }

    // precompute A_i, B_i, C_i, E_i(t_i)
    prepareProbComputation();

    // variable declarations and initialization
    double lnProbTimes = computeLnProbabilityTimes();
    
    return lnProbTimes;
}


/**
 * Compute the log probability of the current value under the current parameter values.
 */
double BirthDeathSamplingTreatmentProcess::computeLnProbabilityTimes( void ) const
{
    
    // variable declarations and initialization
    double lnProbTimes = 0.0;
    size_t num_initial_lineages = 0;
    TopologyNode* root = &value->getRoot();

    // Make sure there is only one kind of event
    // Only check if we have any events of any kind, this means we can avoid checking in purely serial models
    if ( heterogeneous_Lambda != NULL || heterogeneous_Mu != NULL || heterogeneous_Phi != NULL)
    {
        for (size_t i=0; i<global_timeline.size(); ++i)
        {
            if ( (phi_event[i] > DBL_EPSILON && (mu_event[i] > DBL_EPSILON || lambda_event[i] > DBL_EPSILON)) || (mu_event[i] > DBL_EPSILON && lambda_event[i] > DBL_EPSILON) )
            {
                return RbConstants::Double::neginf;
            }
        }
    }

    if ( use_origin == true )
    {
        // If we are conditioning on survival from the origin,
        // then we must divide by 2 the log survival probability computed by AbstractBirthDeathProcess
        // TODO: Generalize AbstractBirthDeathProcess to allow conditioning on the origin
        num_initial_lineages = 1;

        double t = getOriginAge();
        size_t index = findIndex(t);
        lnProbTimes += lnD(index,t);

        t = value->getRoot().getAge();
        index = findIndex(t);
        lnProbTimes -= lnD(index,t);
    }
    // if conditioning on root, root node must be a "true" bifurcation event
    else
    {
        if ( root->isSampledAncestorTipOrParent() )
        {
            return RbConstants::Double::neginf;
        }

        num_initial_lineages = 2;
    }

    // get node/time variables
    size_t num_nodes = value->getNumberOfNodes();

    // add the event-sampling terms (iia)
    for (size_t i = 0; i < global_timeline.size(); ++i)
    {
        // Only compute extinction probability when there is an extinction event
        if (mu_event[i] > DBL_EPSILON)
        {
            if ( RbMath::isFinite(lnProbTimes) == false )
            {
                return RbConstants::Double::nan;
            }

            // Calculate probability of the survivors
            int active_lineages_at_t = survivors(global_timeline[i]);
                    
            lnProbTimes += active_lineages_at_t * log(1 - mu_event[i]);
        }
        
    }
    // add the event-sampling terms (iib)
    for (size_t i = 0; i < global_timeline.size(); ++i)
    {
        // Only compute sampling probability when there is a sampling event
        if (phi_event[i] > DBL_EPSILON)
        {
            if ( RbMath::isFinite(lnProbTimes) == false )
            {
                return RbConstants::Double::nan;
            }

            // Calculate probability of the samples
            double ln_sampling_event_prob = 0.0;
            int S_i = int(event_sampled_ancestor_ages[i].size());
            int T_i = int(event_tip_ages[i].size());
            int I_i = S_i + T_i;
            int L_i = survivors(global_timeline[i]); //A(t_{\rho_i})
            
            // Make sure that we aren't claiming to have sampled all lineages without having sampled all lineages
            if (phi_event[i] >= (1.0 - DBL_EPSILON) && (L_i != I_i) )
            {

                return RbConstants::Double::neginf;
            }
            else
            {
                ln_sampling_event_prob += I_i * log(phi_event[i]);

                // Instead of adding the sampling probability to ln_D we add it here.
                if ( i > 0 && (L_i - I_i) > 0 )
                {
                    ln_sampling_event_prob += (L_i - I_i) * log(1 - phi_event[i]);
                }

            }

            // Calculate probability of the sampled ancestors
            if ( r_event[i] > (1.0 - DBL_EPSILON) && S_i > 0 )
            {
                // Cannot have sampled ancestors if r(t) == 1
                return RbConstants::Double::neginf;
            }
            if ( global_timeline[i] > DBL_EPSILON )
            {
                // only add these terms for sampling that is not at the present
                if ( S_i > 0 )
                {
                    ln_sampling_event_prob += S_i * log(1 - r_event[i]);
                }
                if ( T_i > 0 )
                { 
                    ln_sampling_event_prob += T_i * log(r_event[i] + (1 - r_event[i])*E_previous[i]);
                }
                
            }
            lnProbTimes += ln_sampling_event_prob;
        }
    }

    // add the serial tip age terms (iii)
    for (size_t i = 0; i < serial_tip_ages.size(); ++i)
    {
        if ( RbMath::isFinite(lnProbTimes) == false )
        {
            return RbConstants::Double::nan;
        }

        double t = serial_tip_ages[i];
        size_t index = findIndex(t);

        // add the log probability for the serial sampling events
        if ( phi[index] == 0.0 )
        {
            return RbConstants::Double::neginf;
        }
        else
        {
            double this_prob = r[index] + (1 - r[index]) * E(index,t);
            this_prob *= phi[index];
            lnProbTimes += log( this_prob );
        }
    }

    // add the serial sampled ancestor terms (iv)
    for (size_t i=0; i < serial_sampled_ancestor_ages.size(); ++i)
    {
        if ( RbMath::isFinite(lnProbTimes) == false )
        {
            return RbConstants::Double::nan;
        }

        double t = serial_sampled_ancestor_ages[i];
        size_t index = findIndex(t);

        if ( r[index] > 1.0 - DBL_EPSILON )
        {
            return RbConstants::Double::neginf;
        }

        lnProbTimes += log(phi[index]) + log(1 - r[index]);

    }

    // add the burst bifurcation age terms (v)
    for (size_t i = 0; i < global_timeline.size(); ++i)
    {
        
        if ( lambda_event[i] > RbSettings::userSettings().getTolerance() || event_bifurcation_times[i].size() > 0 )
        {

            lnProbTimes += event_bifurcation_times[i].size() * log(lambda_event[i]);

            // Instead of adding the burst probability to ln_D we add it here.
            int active_lineages_at_t = survivors(global_timeline[i]); //A(t_{\rho_i})
            int A_minus_K = active_lineages_at_t - int(event_bifurcation_times[i].size());
            lnProbTimes += A_minus_K * log(2*lambda_event[i]*E_previous[i]+(1.0 - lambda_event[i]));

        }

    }

    // add the non-burst bifurcation age terms (vi)
    for (size_t i = 0; i < serial_bifurcation_times.size(); ++i)
    {
        if ( RbMath::isFinite(lnProbTimes) == false )
        {
            return RbConstants::Double::nan;
        }

        double t = serial_bifurcation_times[i];
        size_t index = findIndex(t);
        lnProbTimes += log( lambda[index] );
    }

    // Compute probabilities of branch segments on all branches
    for (size_t i=0; i<num_nodes; ++i)
    {
        const TopologyNode& n = value->getNode( i );

        double t = n.getAge();
        
        size_t index = findIndex(t);
        // @TODO: We may need to check more carefully for the boundary of the epochs!
        double diff = t - global_timeline[index];
        if ( index > 0 && fabs(diff) < 1E-4 )
        {
            --index;
        }
        double this_ln_D = lnD(index,t);

        if ( n.isTip() )
        {
            lnProbTimes -= this_ln_D;
        }
        else
        {
            lnProbTimes += this_ln_D;
        }

    }
    lnProbTimes += lnD(findIndex(value->getRoot().getAge()),value->getRoot().getAge());

    // condition on survival
    if ( condition == "survival" )
    {
        double age = use_origin ? getOriginAge() : value->getRoot().getAge();
        // conditioning on survival depends on if we are using the origin or root age
        // origin: we condition on a single lineage surviving to the present and being sampled
        // root: we condition on the above plus the other root child leaving a sampled descendant
        
        lnProbTimes -= num_initial_lineages * log( pSurvival(age,0.0) );
    }
    else if ( condition == "sampling" )
    {
        // conditioning on sampling depends on if we are using the origin or root age
        // origin: the conditioning suggested by Stadler 2011 and used by Gavryuskina (2014), sampling at least one lineage
        // root age: sampling at least one descendent from each child of the root
        double age = use_origin ? getOriginAge() : value->getRoot().getAge();
        lnProbTimes -= num_initial_lineages * log( pSampling(age) );
    }

    if ( RbMath::isFinite(lnProbTimes) == false )
    {
        return RbConstants::Double::nan;
    }


    return lnProbTimes;
}

/*
 * Counts all nodes in a number of sets
 * 1) bifurcation (birth) events NOT at burst times
 * 2) bifurcation (birth) events at burst times
 * 3) sampled ancestor nodes NOT at sampling events
 * 4) sampled nodes from extinct lineages NOT at sampling events
 * 5) sampled ancestor nodes at sampling events
 * 6) sampled nodes from extinct lineages at sampling events
 *
 * Non-burst trackers (1,3,4) are vectors of times of the samples.
 * All burst trackers (2,5,6) are vectors of vectors of samples, each vector corresponding to an event
 */
bool BirthDeathSamplingTreatmentProcess::countAllNodes(void) const
{
  // get node/time variables
  size_t num_nodes = value->getNumberOfNodes();

  size_t num_extant_taxa = 0;

  // clear current
  serial_tip_ages.clear();
  serial_sampled_ancestor_ages.clear();
  event_tip_ages.clear();
  event_sampled_ancestor_ages.clear();
  serial_bifurcation_times.clear();
  event_bifurcation_times.clear();

  // initialize vectors of vectors
  event_sampled_ancestor_ages = std::vector<std::vector<double> >(global_timeline.size(),std::vector<double>(0,1.0));
  event_tip_ages = std::vector<std::vector<double> >(global_timeline.size(),std::vector<double>(0,1.0));
  event_bifurcation_times = std::vector<std::vector<double> >(global_timeline.size(),std::vector<double>(0,1.0));

  // Assign all node times (bifurcations, sampled ancestors, and tips) to their sets
  for (size_t i = 0; i < num_nodes; i++)
  {
      const TopologyNode& n = value->getNode( i );

      double t = n.getAge();

      if ( n.isTip() && n.isFossil() && n.isSampledAncestorTip() )
      {
        // node is sampled ancestor
          int at_event = whichIntervalTime(t);

          // If this tip is not at an event time (and specifically at an event time with Phi[i] > 0), it's a serial tip
          if (at_event == -1 || phi_event[at_event] < DBL_EPSILON)
          {
              serial_sampled_ancestor_ages.push_back(t);
          }
          else
          {
              event_sampled_ancestor_ages[at_event].push_back(t);
          }
      }
      else if ( n.isTip() && n.isFossil() && !n.isSampledAncestorTip() )
      {
          // node is serial leaf
          int at_event = whichIntervalTime(t);

          // If this tip is not at an event time (and specifically at an event time with Phi[i] > 0), it's a serial tip
          if (at_event == -1 || phi_event[at_event] < DBL_EPSILON)
          {
              serial_tip_ages.push_back(t);
          }
          else
          {
              event_tip_ages[at_event].push_back(t);
          }
      }
      else if ( n.isTip() && !n.isFossil() )
      {
          // Node is at present, this can happen even if Phi[0] = 0, so we check if there is really a sampling event at the present
          if (phi_event[0] >= DBL_EPSILON)
          {
              // node is extant leaf
              num_extant_taxa++;
              event_tip_ages[0].push_back(0.0);
          }
          else
          {
              serial_tip_ages.push_back(0.0);
          }
      }
      else if ( n.isInternal() && !n.getChild(0).isSampledAncestorTip() && !n.getChild(1).isSampledAncestorTip() )
      {
          if ( n.isRoot() == false || use_origin == true )
          {
              // node is bifurcation event (a "true" node)
              int at_event = whichIntervalTime(t);

              // If this bifurcation is not at an event time (and specifically at an event time with Lambda[i] > 0), it's a serial bifurcation
              if ( at_event == -1 || lambda_event[at_event] < DBL_EPSILON)
              {
                  serial_bifurcation_times.push_back(t);
              }
              else
              {
                  event_bifurcation_times[at_event].push_back(t);
              }
          }
      }
      else if ( n.isInternal() && n.getChild(0).isSampledAncestorTip() && n.getChild(1).isSampledAncestorTip() )
      {
          return true;
      }
  }

  return false;
}


/**
 * Compute ln(D_i(t))
 */
double BirthDeathSamplingTreatmentProcess::lnD(size_t i, double t) const
{
    // D(0) = 1
    if ( t < DBL_EPSILON )
    {
        return 0.0;
    }
    else
    {
        double s = global_timeline[i];
        double this_lnD_i = 0.0;
        if (i > 0)
        {
            this_lnD_i = lnD_previous[i];
        }
        this_lnD_i += 2*RbConstants::LN2 + (-A_i[i] * (t - s));
        this_lnD_i -= 2 * log(1 + B_i[i] + exp(-A_i[i] * (t - s)) * (1 - B_i[i]));

        return this_lnD_i;
    }
}

/**
 * Compute E_i(t)
 */
double BirthDeathSamplingTreatmentProcess::E(size_t i, double t, bool computeSurvival) const
{
    double E_i;
    double s = global_timeline[i];

    // Are we computing E(t) for survival conditioning?
    if (computeSurvival == true)
    {
        E_i = lambda[i] + mu[i] + phi[i]*r[i];
        E_i -= A_survival_i[i] * (1 + B_survival_i[i] - exp(-A_survival_i[i] * (t - s)) * (1 - B_survival_i[i])) / (1 + B_survival_i[i] + exp(-A_survival_i[i] * (t - s)) * (1 - B_survival_i[i]));
        E_i /= (2 * lambda[i]);
    }
    else
    {
        E_i = lambda[i] + mu[i] + phi[i];
        E_i -= A_i[i] * (1 + B_i[i] - exp(-A_i[i] * (t - s)) * (1 - B_i[i])) / (1 + B_i[i] + exp(-A_i[i] * (t - s)) * (1 - B_i[i]));
        E_i /= (2 * lambda[i]);
    }

    return E_i;

}

/**
 * return the index i so that s_{i-1} <= t < s_i
 * where s_i is the global timeline of events
 * s_0 = 0.0
 * s_l = Inf
 */
size_t BirthDeathSamplingTreatmentProcess::findIndex(double t) const
{
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

/**
 * return the index i so that x_{i-1} <= t < x_i
 * where x is one of the input vector timelines
 */
size_t BirthDeathSamplingTreatmentProcess::findIndex(double t, const std::vector<double> &timeline) const
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

// calculate offset so we can set t_0 to time of most recent tip
void BirthDeathSamplingTreatmentProcess::getOffset(void) const
{
    // On first pass, there is no tree, so we can't loop over nodes
    // Get taxon ages directly from taxa instead
    if ( value->getNumberOfNodes() == 0 )
    {
        offset = RbConstants::Double::max;
        for (size_t i = 0; i < taxa.size(); i++)
        {
            const Taxon& n = taxa[i];

            if ( n.getAge() < offset )
            {
                offset = n.getAge();
            }
        }
    }
    // On later passes we have the tree, to avoid any issues with tree and taxon age mismatch, get ages from tree
    else
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

}

bool BirthDeathSamplingTreatmentProcess::isConstantRate(void) const
{
  bool has_no_interval_times = false;
  // For there to be no intervals, every timeline must either be NULL or have size 0
  if ( (interval_times_global == NULL           || interval_times_global->getValue().size() == 0 ) &&
       (interval_times_speciation == NULL       || interval_times_speciation->getValue().size() == 0 ) &&
       (interval_times_extinction == NULL       || interval_times_extinction->getValue().size() == 0 ) &&
       (interval_times_sampling == NULL         || interval_times_sampling->getValue().size() == 0 ) &&
       (interval_times_treatment == NULL        || interval_times_treatment->getValue().size() == 0 ) &&
       (interval_times_event_speciation == NULL || interval_times_event_speciation->getValue().size() == 0 ) &&
       (interval_times_event_extinction == NULL || interval_times_event_extinction->getValue().size() == 0 ) &&
       (interval_times_event_sampling == NULL   || interval_times_event_sampling->getValue().size() == 0 ) )
  {
    has_no_interval_times = true;
  }

  bool all_parameters_are_scalars = false;
  // For all parameters to be scalars,
  // 1) rate parameters must either be homogenous or they must have size <= 1 (1 for scalar, 0 if it's null)
  // 2) Lambda/Mu must be of size 0 or NULL
  // 3) Phi must be of size 1 or a scalar
  if ( (heterogeneous_lambda == NULL || lambda.size() <= 1) &&
       (heterogeneous_mu == NULL     || mu.size() <= 1) &&
       (heterogeneous_phi == NULL    || phi.size() <= 1) &&
       (heterogeneous_r == NULL      || r.size() <= 1) &&
       (heterogeneous_Lambda == NULL || lambda_event.size() == 0) &&
       (heterogeneous_Mu     == NULL || mu_event.size() == 0) &&
       (heterogeneous_Phi == NULL    || phi_event.size() <= 1) &&
       (heterogeneous_R      == NULL || r_event.size() == 0) )
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

// double BirthDeathSamplingTreatmentProcess::lnProbNumTaxa(size_t n, double start, double end, bool MRCA) const
// {
//   throw(RbException("Cannot compute lnProbNumTaxa."));
//   return RbConstants::Double::neginf;
// }

double BirthDeathSamplingTreatmentProcess::lnProbTreeShape(void) const
{
    // the birth death divergence times density is derived for a (ranked) unlabeled oriented tree
    // so we convert to a (ranked) labeled non-oriented tree probability by multiplying by 2^{n+m-1} / (n+m)!
    // where n is the number of extant tips, m is the number of extinct tips

    int num_taxa = (int)value->getNumberOfTips();
    int num_sa = (int)value->getNumberOfSampledAncestors();

    //Gavryushkina (2014) uses the following
    return (num_taxa - num_sa - 1) * RbConstants::LN2 - RbMath::lnFactorial(num_taxa);
}

/**
 * Takes a par.size() < global_timeline.size() vector and makes it the correct size to work with our global timeline.
 * The parameter has its own reference timeline, which we use to find the rate in the global intervals.
 * This works only for parameters Lambda,Mu,Phi,R where missing values are 0.0
 */
void BirthDeathSamplingTreatmentProcess::expandNonGlobalProbabilityParameterVector(std::vector<double> &par, const std::vector<double> &par_times) const
{
    // @TODO @efficiency: this works but it would be faster to auto-advance indices rather than have an internal loop
    // Store the original values so we can overwrite the vector
    std::vector<double> old_par = par;

    // For each time in the global timeline, find the rate according to this variable's own timeline
    for (size_t i=0; i<global_timeline.size(); ++i)
    {
        bool global_time_is_variable_time = false;
        for (size_t j=0; i<par_times.size(); ++j)
        {
            if ( fabs(par_times[j] - global_timeline[j]) < DBL_EPSILON )
            {
                // time is in variable's timeline
                par[i] = old_par[j];
                global_time_is_variable_time = true;
                break;
            }
        }

        // Time is not in variable's own timeline, probability of event here is 0
        if ( !global_time_is_variable_time )
        {
            par[i] = 0.0;
        }
    }

}

/**
 * Takes a par.size() < global_timeline.size() vector and makes it the correct size to work with our global timeline.
 * The parameter has its own reference timeline, which we use to find the rate in the global intervals.
 * This works only for parameters (lambda,mu,phi,r), where the global timeline is simply a finer grid than the variable-specific timelines.
 */
void BirthDeathSamplingTreatmentProcess::expandNonGlobalRateParameterVector(std::vector<double> &par, const std::vector<double> &par_times) const
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



/*
 * Here we calculate all A_i, B_i, C_i, D_i(s_i), and E_i(s_i) for i = 1,...,l
 *
 */
void BirthDeathSamplingTreatmentProcess::prepareProbComputation( void ) const
{
    // TODO: B and C are producing nan values upon initialization, but not when computing tree probabilities, which are working fine

    // re-initialize all the sets
    A_i = std::vector<double>(global_timeline.size(),0.0);
    B_i = std::vector<double>(global_timeline.size(),0.0);
    C_i = std::vector<double>(global_timeline.size(),0.0);

    // E_previous[i] is E_{i-1}(s_i)
    // lnD_previous[i] is log(D_{i-1}(s_i))
    E_previous      = std::vector<double>(global_timeline.size(),0.0);
    lnD_previous    = std::vector<double>(global_timeline.size(),0.0);

    // timeline[0] == 0.0
    double t = global_timeline[0];

    // Compute all starting at 1
    A_i[0] = sqrt( pow(lambda[0] - mu[0] - phi[0],2.0) + 4 * lambda[0] * phi[0]);

    // At the present, only sampling is allowed, no birth/death bursts
    C_i[0] = (1 - phi_event[0]);

    B_i[0] = (1.0 - 2.0 * C_i[0]) * lambda[0] + mu[0] + phi[0];
    B_i[0] /= A_i[0];

    // E_{i-1}(0) = 1, and our E(i,t) function requires i >= 0, so we hard-code this explicitly
    E_previous[0] = 1.0;

    // we always initialize the probability of observing the lineage at the present with the sampling probability
    lnD_previous[0] = 0.0;

    for (size_t i=1; i<global_timeline.size(); ++i)
    {
        t = global_timeline[i];

        // first, we need to compute E and D at the end of the previous interval
        E_previous[i] = E(i-1, t, false);
        lnD_previous[i] = lnD(i-1, t);

        // now we can compute A_i, B_i and C_i at the end of this interval.
        A_i[i] = sqrt( pow(lambda[i] - mu[i] - phi[i],2.0) + 4 * lambda[i] * phi[i]);

        // Only one type of event is allowed
        if ( phi_event[i] >= DBL_EPSILON )
        {
            C_i[i] = (1 - phi_event[i]) * E_previous[i];
        }
        else if ( lambda_event[i] >= DBL_EPSILON )
        {
            C_i[i] = (1 - lambda_event[i]) * E_previous[i] + lambda_event[i] * E_previous[i] * E_previous[i];
        }
        else if ( mu_event[i] >= DBL_EPSILON )
        {
            C_i[i] = (1 - mu_event[i]) * E_previous[i] + mu_event[i];
        }
        else
        {
            C_i[i] = E_previous[i];
        }

        B_i[i] = (1.0 - 2.0 * C_i[i]) * lambda[i] + mu[i] + phi[i];
        B_i[i] /= A_i[i];

    }

    // if we want to condition on survival we need to track versions of A,B,C,E that set phi = 0 and Phi[-0] to 0
    if ( condition == "survival" )
    {
        // timeline[0] == 0.0
        double t = global_timeline[0];

        A_survival_i = std::vector<double>(global_timeline.size(),0.0);
        B_survival_i = std::vector<double>(global_timeline.size(),0.0);
        C_survival_i = std::vector<double>(global_timeline.size(),0.0);
        E_survival_previous = std::vector<double>(global_timeline.size(),0.0);

        // Compute all starting at 1
        A_survival_i[0] = sqrt( pow(lambda[0] - mu[0] - phi[0]*r[0],2.0));

        // At the present, only sampling is allowed, no birth/death bursts
        C_survival_i[0] = (1 - phi_event[0]);

        B_survival_i[0] = (1.0 - 2.0 * C_survival_i[0]) * lambda[0] + mu[0] + phi[0]*r[0];
        B_survival_i[0] /= A_survival_i[0];
        
        E_survival_previous[0] = 1.0;

        for (size_t i=1; i<global_timeline.size(); ++i)
        {
            t = global_timeline[i];

            // first, we need to compute E and D at the end of the previous interval
            E_survival_previous[i] = E(i-1, t, true);

            // now we can compute A_survival_i, B_survival_i and C_survival_i at the end of this interval.
            A_survival_i[i] = lambda[i] - mu[i] - phi[i]*r[i];

            // Only one type of event is allowed
            if ( lambda_event[i] >= DBL_EPSILON )
            {
                C_survival_i[i] =  (1 - lambda_event[i]) * E_survival_previous[i] + lambda_event[i] * E_survival_previous[i] * E_survival_previous[i];
            }
            else if ( mu_event[i] >= DBL_EPSILON )
            {
                C_survival_i[i] = (1 - mu_event[i]) * E_survival_previous[i] + mu_event[i];
            }
            else
            {
                C_survival_i[i] = E_survival_previous[i];
            }

            B_survival_i[i] = (1.0 - 2.0 * C_survival_i[i]) * lambda[i] + mu[i] + phi[i]*r[i];
            B_survival_i[i] /= A_survival_i[i];
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
void BirthDeathSamplingTreatmentProcess::prepareTimeline( void ) const
{
    // clean all the sets
    lambda.clear();
    mu.clear();
    phi.clear();
    r.clear();
    lambda_event.clear();
    mu_event.clear();
    phi_event.clear();
    r_event.clear();

    lambda_times.clear();
    mu_times.clear();
    phi_times.clear();
    r_times.clear();
    lambda_event_times.clear();
    mu_event_times.clear();
    phi_event_times.clear();
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
    if (heterogeneous_phi != NULL)
    {
        phi = heterogeneous_phi->getValue();
    }
    if (heterogeneous_r != NULL)
    {
        r = heterogeneous_r->getValue();
    }
    if (heterogeneous_Lambda != NULL)
    {
        lambda_event = heterogeneous_Lambda->getValue();
    }
    if (heterogeneous_Mu != NULL)
    {
        mu_event = heterogeneous_Mu->getValue();
    }
    if (heterogeneous_Phi != NULL)
    {
        phi_event = heterogeneous_Phi->getValue();
    }
    if (heterogeneous_R != NULL)
    {
        r_event = heterogeneous_R->getValue();
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
        phi_times = interval_times_sampling->getValue();
    }
    if (interval_times_treatment != NULL)
    {
        r_times = interval_times_treatment->getValue();
    }
    if (interval_times_event_speciation != NULL)
    {
        lambda_event_times = interval_times_event_speciation->getValue();
    }
    if (interval_times_event_extinction != NULL)
    {
        mu_event_times = interval_times_event_extinction->getValue();
    }
    if (interval_times_event_sampling != NULL)
    {
        phi_event_times = interval_times_event_sampling->getValue();
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
             interval_times_event_speciation != NULL ||
             interval_times_event_sampling != NULL ||
             interval_times_event_extinction != NULL )
        {
            throw RbException("Both heterogeneous and homogeneous rate change times provided");
        }

        // check that the number of provided parameters matches the global timeline
        // Right now, the global timeline is only the interval times, i.e. breaks between pieces/episodes/windows
        checkVectorSizes(heterogeneous_lambda,interval_times_global,1,spn,true);
        checkVectorSizes(heterogeneous_mu,interval_times_global,1,exn,true);
        checkVectorSizes(heterogeneous_phi,interval_times_global,1,smp,true);
        checkVectorSizes(heterogeneous_r,interval_times_global,1,trt,false);
        checkVectorSizes(heterogeneous_Lambda,interval_times_global,0,spn,false);
        checkVectorSizes(heterogeneous_Mu,interval_times_global,0,exn,false);
        checkVectorSizes(heterogeneous_Phi,interval_times_global,1,smp,false);
        checkVectorSizes(heterogeneous_R,interval_times_global,0,etrt,false);

        sortGlobalTimesAndVectorParameter();

        // we are done with setting up the timeline (i.e., using the provided global timeline) and checking all dimension of parameters
    }
    // We only need to assemble a global timeline if
    else if ( interval_times_speciation != NULL ||
                 interval_times_extinction != NULL ||
                 interval_times_sampling != NULL ||
                 interval_times_treatment != NULL ||
                 interval_times_event_speciation != NULL ||
                 interval_times_event_sampling != NULL ||
                 interval_times_event_extinction != NULL )
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
        if ( heterogeneous_phi == NULL && homogeneous_phi == NULL)
        {
            throw RbException("Sampling rate must be of type RealPos or RealPos[]");
        }
        else if ( heterogeneous_phi != NULL )
        {
            if ( interval_times_sampling == NULL ) throw RbException("No time intervals provided for piecewise constant sampling rates");
            checkVectorSizes(heterogeneous_phi,interval_times_sampling,1,smp,true);
            sortNonGlobalTimesAndVectorParameter(phi_times,phi);
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

        // check if correct number of burst probabilities were provided
        // if provided as a vector, sort to the correct timescale
        if ( heterogeneous_Lambda != NULL )
        {
            if ( interval_times_event_speciation == NULL ) throw RbException("No time intervals provided for speciation bursts");
            checkVectorSizes(heterogeneous_Lambda,interval_times_event_speciation,0,spn,false);
            sortNonGlobalTimesAndVectorParameter(lambda_event_times,lambda_event);
        }

        // check if correct number of mass extinction probabilities were provided
        // if provided as a vector, sort to the correct timescale
        if ( heterogeneous_Mu != NULL )
        {
            if ( interval_times_event_extinction == NULL ) throw RbException("No time intervals provided for piecewise constant extinction rates");
            checkVectorSizes(heterogeneous_Mu,interval_times_event_extinction,0,exn,false);
            sortNonGlobalTimesAndVectorParameter(mu_event_times,mu_event);
        }

        // check if correct number of sampling probabilities were provided
        // if provided as a vector, sort to the correct timescale
        if ( heterogeneous_Phi == NULL && homogeneous_Phi == NULL)
        {
          throw RbException("Event-sampling probabilities must be of type Probability or Probability[]");
        }
        else if ( heterogeneous_Phi != NULL )
        {
            if ( interval_times_event_sampling == NULL ) throw RbException("No time intervals provided for piecewise constant sampling rates");
            checkVectorSizes(heterogeneous_Phi,interval_times_event_sampling,1,smp,false);
            sortNonGlobalTimesAndVectorParameter(phi_event_times,phi_event);
        }

        // check if correct number of treatment probabilities at event-sampling times were provided
        // there must be one fewer value of R than of Phi
        // we default to R = r if both heterogeneous_R and homogeneous_R are NULL
        if ( homogeneous_Phi != NULL )
        {
            // If there's one sampling event, there cannot be any event treatment probabilities
            if ( !(heterogeneous_R == NULL || heterogeneous_R->getValue().size() == 0) )
            {
                throw RbException("Number of event treatment probabilities does not match number of sampling events");
            }
        }
        else if ( heterogeneous_Phi != NULL )
        {
            if ( interval_times_event_sampling == NULL ) throw RbException("No time intervals provided for event sampling probabilities");
            checkVectorSizes(heterogeneous_R,interval_times_event_sampling,0,etrt,false);
            // This should be sorted to match phi_event_times, which is already sorted, so we copy the original value to sort against
            std::vector<double> tmp_phi_event_times = heterogeneous_Phi->getValue();
            if ( heterogeneous_R != NULL )
            {
                sortNonGlobalTimesAndVectorParameter(tmp_phi_event_times,r_event);
            }
        }

        // now we start assembling the global timeline by finding the union of unique intervals for all parameters
        std::set<double> event_times;
        addTimesToGlobalTimeline(event_times,interval_times_speciation);
        addTimesToGlobalTimeline(event_times,interval_times_extinction);
        addTimesToGlobalTimeline(event_times,interval_times_sampling);
        addTimesToGlobalTimeline(event_times,interval_times_treatment);
        addTimesToGlobalTimeline(event_times,interval_times_event_speciation);
        addTimesToGlobalTimeline(event_times,interval_times_event_extinction);
        addTimesToGlobalTimeline(event_times,interval_times_event_sampling);


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
    if ( heterogeneous_phi != NULL )
    {
        if ( phi.size() != global_timeline.size() )
        {
            expandNonGlobalRateParameterVector(phi,phi_times);
        } // else it matches in size and is already sorted and is thus ready to be used
    }
    else
    {
        phi = std::vector<double>(global_timeline.size(),homogeneous_phi->getValue());
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

    // For each parameter vector, we now make sure that its size matches the size of the global vector
    // For Lambda/Mu, there are two cases
    //     1) It is a vector and is is of length global_timeline.size() - 1, in which case we add an event with probability 0.0 at the present, and it is ready to use
    //     2) It is a vector and it DOES NOT match the size of the global timeline, in which case we must expand it to match, which automatically adds an event of P=0.0 at the present

    // Get vector of burst birth probabilities
    if ( heterogeneous_Lambda != NULL )
    {
        // Expand if needed, this will make the first event 0
        if (lambda_event.size() != global_timeline.size() - 1)
        {
            expandNonGlobalProbabilityParameterVector(lambda_event,lambda_event_times);
        }
        else
        {
            // Add first event. lambda_event_0 must be 0 (there can be no burst at the present)
            lambda_event.insert(lambda_event.begin(),0.0);
        }
    }
    else
    {
        // User specified nothing, there are no birth bursts
        lambda_event = std::vector<double>(global_timeline.size(),0.0);
    }

    // Get vector of burst death (mass extinction) probabilities
    if ( heterogeneous_Mu != NULL )
    {
        // Expand if needed, this will make the first event 0
        if (mu_event.size() != global_timeline.size() - 1)
        {
            expandNonGlobalProbabilityParameterVector(mu_event,mu_event_times);
        }
        else
        {
            // mu_event_0 must be 0 (there can be no burst at the present)
            mu_event.insert(mu_event.begin(),0.0);
        }
    }
    else
    {
        // User specified nothing, there are no birth bursts
         mu_event = std::vector<double>(global_timeline.size(),0.0);
    }

    // Get vector of event sampling probabilities
    // For Phi, there are three cases
    //     1) It is a vector and it matches the size of the global timeline, in which case it is already sorted and we can use it
    //     2) It is a vector and it DOES NOT match the size of the global timeline, in which case we must expand it to match
    //     3) It is a scalar, in which case it is Phi[0] and we simply make Phi[>0] all 0.0
    if ( heterogeneous_Phi != NULL )
    {
        // Expand if needed
        if (phi_event.size() != global_timeline.size())
        {
            expandNonGlobalProbabilityParameterVector(phi_event,phi_event_times);
        }
    }
    else
    {
        phi_event = std::vector<double>(global_timeline.size(),0.0);
        if ( homogeneous_Phi != NULL )
        {
            // User specified the sampling fraction at the present
            phi_event[0] = homogeneous_Phi->getValue();
        }
    }

    // Get vector of burst death (mass extinction) probabilities
    // For R, the cases are as follows
    //     1) It is a vector and it is of length phi_event.size() - 1, in which case we simply add R[0] = 0.0 and we can move on
    //     2) It is a vector and it DOES NOT match the size of the global timeline, in which case we must expand it to match
    //     3) It is NULL, in which case we use r in its place
    if ( heterogeneous_R != NULL )
    {
        // Expand if needed, this will make the first event 0
        if (r_event.size() != global_timeline.size() - 1)
        {
            expandNonGlobalProbabilityParameterVector(r_event,phi_event_times);
        }
        else
        {
            // treatment_event_0 must be 0 (there can be no burst at the present)
            r_event.insert(r_event.begin(),0.0);
        }
    }
    else
    {
        // @TODO: @ANDY: Needs revision
        // User specified nothing, there are no birth bursts
        r_event = r;
    }
}

double BirthDeathSamplingTreatmentProcess::pSampling(double start) const
{
  return (1.0 - E(findIndex(start),start,false));
}

double BirthDeathSamplingTreatmentProcess::pSurvival(double start, double end) const
{
    // This computation does not make sense unless there is sampling at the present, and will result in a divide by 0 error
    if ( phi_event[0] < DBL_EPSILON )
    {
        return RbConstants::Double::neginf;
    }
    return (1.0 - E(findIndex(start),start,true))/(1.0 - E(findIndex(end),end,true));
}



void BirthDeathSamplingTreatmentProcess::redrawValue( SimulationCondition condition )
{

    if ( condition == SimulationCondition::MCMC )
    {
        if ( starting_tree == NULL )
        {
            // SH 20221212: The simulateTree functions hangs in certain situations. It's more robust to use the coalescent simulator.
            // simulateTree();
            
            RbVector<Clade> constr;
            // We employ a coalescent simulator to guarantee that the starting tree matches all time constraints
            StartingTreeSimulator simulator;
            RevBayesCore::Tree *my_tree = simulator.simulateTree( taxa, constr );
            
            // When conditioning on the root age, increase the root age of the simulated tree to make it equal to the process age
            if ( use_origin == false )
            {
                double new_root_age = getOriginAge();
                my_tree->getRoot().setAge( new_root_age );
            }
            
            // store the new value
            value = my_tree;
        }
    }
    else if ( condition == SimulationCondition::VALIDATION )
    {
        // update timeline and parameter vectors
        prepareTimeline();
        
        BirthDeathForwardSimulator simulator;
        simulator.setMaxNumLineages(50000);
        
        size_t num_epochs = global_timeline.size();
        std::vector< std::vector<double> > tmp = std::vector< std::vector<double> >( num_epochs, std::vector<double>(1,0) );

        for (size_t i=0; i<num_epochs; ++i) tmp[i][0] = lambda_event[i];
        simulator.setBurstProbability( tmp );
        
        for (size_t i=0; i<num_epochs; ++i) tmp[i][0] = mu[i];
        simulator.setExtinctionRate( tmp );
        
        for (size_t i=0; i<num_epochs; ++i) tmp[i][0] = mu_event[i];
//        for (size_t i=0; i<num_epochs; ++i) tmp[i][0] = 0.0;
        simulator.setMassExtinctionProbability( tmp );
        
        for (size_t i=0; i<num_epochs; ++i) tmp[i][0] = phi_event[i];
        simulator.setSamplingProbability( tmp );
        
        for (size_t i=0; i<num_epochs; ++i) tmp[i][0] = r_event[i];
        simulator.setSamplingExtinctionProbability( tmp );
        
        for (size_t i=0; i<num_epochs; ++i) tmp[i][0] = phi[i];
        simulator.setSamplingRate( tmp );
        
        for (size_t i=0; i<num_epochs; ++i) tmp[i][0] = r[i];
        simulator.setSamplingExtinctionRate( tmp );
        
        for (size_t i=0; i<num_epochs; ++i) tmp[i][0] = lambda[i];
        simulator.setSpeciationRate( tmp );
        
        simulator.setTimeline( global_timeline );
        
        
        simulator.setRootCategoryProbabilities( std::vector<double>(1,1) );
        
        do {
            Tree *my_tree = simulator.simulateTreeConditionTime( getOriginAge(), BirthDeathForwardSimulator::SIM_CONDITION::ROOT);
            
            // store the new value
            delete value;
            value = my_tree;
            
            // to be safe, we copy over the taxa
            taxa = value->getTaxa();

        } while ( RbMath::isFinite( computeLnProbability() ) == false );
    }
    else
    {
        throw RbException("Uknown condition for simulating tree in episodic birth-death-sampling-treatment process.");
    }
}


/**
 * Simulate new speciation times.
 */
double BirthDeathSamplingTreatmentProcess::simulateDivergenceTime(double origin, double present) const
{
    // incorrect placeholder, there is no way to simulate an FBD tree consistent with fossil times, we use a coalescent simulator instead

    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;

    size_t i = findIndex(present);

    // get the parameters
    double age = origin - present;
    double b = lambda[i];
    double d = mu[i];
    double p_e = phi_event[i];
    double x = b - d;

    // make sure age is not negative, otherwise function doesn't work
    assert(age >= 0);

    // get a random draw
    double u = rng->uniform01();

    // compute the time for this draw
    double t = 0.0;
    if ( x > 0 )
    {
        if( p_e > 0.0 )
        {
            t = ( log( ( x / (1 - (u)*(1-(x*exp((-x)*age))/(p_e*b+(b*(1-p_e)-d)*exp((-x)*age) ) ) ) - (b*(1-p_e)-d) ) / (p_e * b) ) )  /  x;
        }
        else
        {
            t = log( (1 - u) * exp(-x * age) + u) / x + age;
        }
    }
    else
    {
        if( p_e > 0.0 )
        {
            t = ( log( ( x / (1 - (u)*(1-x/(p_e*b*exp(x*age)+(b*(1-p_e)-d) ) ) ) - (b*(1-p_e)-d) ) / (p_e * b) ) )  /  x;
        }
        else
        {
            t = log( 1 - u * (1 - exp(age*x))  ) / x;
        }
    }

    // make sure the result is in the right range
    assert(0 <= t and t <= age);

    return present + t;
}


/**
 * Compute the diversity of the tree at time t.
 *
 * \param[in]    t      time at which we want to know the diversity.
 *
 * \return The diversity (number of species in the reconstructed tree).
 */
int BirthDeathSamplingTreatmentProcess::survivors(double t) const
{  
    
    int survivors = 0;

    if ( use_origin )
    {
        if ( t > getOriginAge() )
        {
            return 0;
        }
        survivors = 1;
    } else {
        if ( t > value->getRoot().getAge() )
        {
            return 0;
        }
        survivors = 2;
    }

    for (size_t i=0; i<serial_bifurcation_times.size(); ++i) {
        if (t < serial_bifurcation_times[i])
        {
            survivors++;
        }
    }

    for (size_t i=0; i<serial_tip_ages.size(); ++i) {
        if (t < serial_tip_ages[i])
        {
            survivors--;
        }
    }

    for (size_t i=0; i<global_timeline.size(); ++i)
    {   
        size_t idx = global_timeline.size() - i - 1;
        if ( global_timeline[idx] < t ) {
            break;
        } else if (global_timeline[idx] > t)
        {   
            // by ignoring time = t we implicitly count all tips at a time as survivors
            // This is compatible with the logic in computing event-sampling probabilities but could be changed
            survivors += (int)event_bifurcation_times[idx].size();
            survivors -= (int)event_tip_ages[idx].size();
        }
    }
    return survivors;
}

/**
 * Sorts global times to run from present to past (0->inf) and orders ALL vector parameters to match this.
 * These can only be sorted after the local copies have values in them.
 */
void BirthDeathSamplingTreatmentProcess::sortGlobalTimesAndVectorParameter() const
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
            if (heterogeneous_phi != NULL)
            {
                sort(phi.rbegin(),phi.rend());
            }
            if (heterogeneous_r != NULL)
            {
                sort(r.rbegin(),r.rend());
            }
            if (heterogeneous_Lambda != NULL)
            {
                sort(lambda_event.rbegin(),lambda_event.rend());
            }
            if (heterogeneous_Mu != NULL)
            {
                sort(mu_event.rbegin(),mu_event.rend());
            }
            if (heterogeneous_Phi != NULL)
            {
                sort(phi_event.rbegin(),phi_event.rend());
            }
            if (heterogeneous_R != NULL)
            {
                sort(r_event.rbegin(),r_event.rend());
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
            if (heterogeneous_phi != NULL)
            {
                std::vector<double> old_phi = phi;
                for (size_t i=0; i<global_timeline.size(); ++i)
                {
                    phi[i] = old_phi[ordering[i]];
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
            if (heterogeneous_Lambda != NULL)
            {
                std::vector<double> old_lambda_event = lambda_event;
                for (size_t i=0; i<global_timeline.size(); ++i)
                {
                    lambda_event[i] = old_lambda_event[ordering[i]];
                }
            }
            if (heterogeneous_Mu != NULL)
            {
                std::vector<double> old_mu_event = mu_event;
                for (size_t i=0; i<global_timeline.size(); ++i)
                {
                    mu_event[i] = old_mu_event[ordering[i]];
                }
            }
            if (heterogeneous_Phi != NULL)
            {
                std::vector<double> old_phi_event = phi_event;
                for (size_t i=0; i<global_timeline.size(); ++i)
                {
                    phi_event[i] = old_phi_event[ordering[i]];
                }
            }
            if (heterogeneous_R != NULL)
            {
                std::vector<double> old_r_event = r_event;
                for (size_t i=0; i<global_timeline.size(); ++i)
                {
                    r_event[i] = old_r_event[ordering[i]];
                }
            }

        }
    }

}

/**
 * Sorts times to run from present to past (0->inf) and orders par to match this.
 */
void BirthDeathSamplingTreatmentProcess::sortNonGlobalTimesAndVectorParameter(std::vector<double> &times, std::vector<double> &par) const
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
 * For a time t, determine which if any interval (event) time it corresponds to
 *
 * \param[in]    t      time we want to assess
 *
 * \return -1 if time matches no interval, otherwise the interval
 */
int BirthDeathSamplingTreatmentProcess::whichIntervalTime(double t) const
{
    // Find the i such that s_i <= t < s_{i+1}
    size_t i = findIndex(t);
    
    // Check if s_i == t
    double diff = t - global_timeline[i];
    if ( fabs(diff) < 1E-4 )
    {
        return (int) i;
    }
    else
    {
        return -1;
    }

}

/**
 * Swap the parameters held by this distribution.
 *
 * \param[in]    oldP      Pointer to the old parameter.
 * \param[in]    newP      Pointer to the new parameter.
 */
void BirthDeathSamplingTreatmentProcess::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    // Rate parameters
    if (oldP == heterogeneous_lambda)
    {
        heterogeneous_lambda = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == heterogeneous_mu)
    {
        heterogeneous_mu = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == heterogeneous_phi)
    {
        heterogeneous_phi = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == homogeneous_lambda)
    {
        homogeneous_lambda = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == homogeneous_mu)
    {
        homogeneous_mu = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == homogeneous_phi)
    {
        homogeneous_phi = static_cast<const TypedDagNode<double>* >( newP );
    }
    // Treatment
    else if (oldP == heterogeneous_r)
    {
        heterogeneous_r = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == homogeneous_r)
    {
        homogeneous_r = static_cast<const TypedDagNode<double>* >( newP );
    }
    // Event probability parameters
    if (oldP == heterogeneous_Lambda)
    {
        heterogeneous_Lambda = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == heterogeneous_Mu)
    {
        heterogeneous_Mu = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == heterogeneous_Phi)
    {
        heterogeneous_Phi = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == homogeneous_Phi)
    {
        homogeneous_Phi = static_cast<const TypedDagNode<double>* >( newP );
    }
    else
    {
        // delegate the super-class
        AbstractBirthDeathProcess::swapParameterInternal(oldP, newP);
    }
}

/**
 * Checks if removal probabilities set for this distribution are compatible with sampled ancestors
 * (i.e. removal < 1)
 */
bool BirthDeathSamplingTreatmentProcess::allowsSA() {
    for(auto removal : r) {
        if(removal < 1.0 - DBL_EPSILON) return true;
    }
    for(auto removal : r_event) {
        if(removal < 1.0 - DBL_EPSILON) return true;
    }
    return false;
}

void RevBayesCore::BirthDeathSamplingTreatmentProcess::setValue(Tree * v, bool f)
{
    RevBayesCore::Tree *newv = TreeUtilities::startingTreeInitializer(*v, taxa, age_check_precision);
    AbstractRootedTreeDistribution::setValue(newv, f);
}
