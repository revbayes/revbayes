#include <stddef.h>
#include <algorithm>
#include <cmath>
#include <iterator>
#include <ostream>
#include <set>
#include <string>
#include <vector>

#include "DistributionExponential.h"
#include "PiecewiseConstantOccurrenceBirthDeathProcess.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbMathCombinatorialFunctions.h"
#include "RbMathLogic.h"
#include "AbstractBirthDeathProcess.h"
#include "DagNode.h"
#include "RbException.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "TopologyNode.h"
#include "Tree.h"
#include "TypedDagNode.h"

#include "ComputeLikelihoodsLtMt.h"

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
PiecewiseConstantOccurrenceBirthDeathProcess::PiecewiseConstantOccurrenceBirthDeathProcess(const TypedDagNode<double> *ra,
                                                                                           const TypedDagNode<double> *rho,
                                                                                           const DagNode *inspeciation,
                                                                                           const DagNode *inextinction,
                                                                                           const DagNode *inpsi,
                                                                                           const DagNode *inomega,
                                                                                           const DagNode *inremovalPr,
                                                                                           const TypedDagNode< RbVector<double> > *ht,
                                                                                           const TypedDagNode< RbVector<double> > *lt,
                                                                                           const TypedDagNode< RbVector<double> > *mt,
                                                                                           const TypedDagNode< RbVector<double> > *pt,
                                                                                           const TypedDagNode< RbVector<double> > *ot,
                                                                                           const TypedDagNode< RbVector<double> > *rt,
                                                                                           const TypedDagNode<long> *n,
                                                                                           const std::string &cdt,
                                                                                           const std::vector<Taxon> &tn,
                                                                                           bool uo,
                                                                                           bool Mt,
                                                                                       TypedDagNode<Tree> *t) : AbstractBirthDeathProcess( ra, cdt, tn, uo ),
    start_age(ra),
    homogeneous_timeline(ht),
    homogeneous_rho(rho),
    lambda_timeline(lt),
    mu_timeline(mt),
    psi_timeline(pt),
    omega_timeline(ot),
    removalPr_timeline(rt),
    maxHiddenLin(n),
    useMt(Mt),
    cond(cdt),
    useOrigin (uo),

    ascending(false)
{
    // initialize all the pointers to NULL
    homogeneous_lambda      = NULL;
    homogeneous_mu          = NULL;
    homogeneous_psi         = NULL;
    homogeneous_omega       = NULL;
    homogeneous_removalPr   = NULL;
    heterogeneous_lambda    = NULL;
    heterogeneous_mu        = NULL;
    heterogeneous_psi       = NULL;
    heterogeneous_omega     = NULL;
    heterogeneous_removalPr = NULL;

    if ( homogeneous_timeline != NULL )
    {
        if ( lambda_timeline != NULL || mu_timeline != NULL || psi_timeline != NULL || omega_timeline != NULL || removalPr_timeline != NULL )
        {
            throw(RbException("Both heterogeneous and homogeneous rate change times provided"));
        }

        std::vector<double> times = homogeneous_timeline->getValue();
        std::vector<double> times_sorted_ascending = times;
        std::vector<double> times_sorted_descending = times;

        sort(times_sorted_ascending.begin(), times_sorted_ascending.end() );
        sort(times_sorted_descending.rbegin(), times_sorted_descending.rend() );

        if( times == times_sorted_ascending )
        {
            ascending = true;
        }
        else if ( times != times_sorted_ascending )
        {
            throw(RbException("Rate change times must be provided in sorted order"));
        }

        addParameter( homogeneous_timeline );
    }
////////
    heterogeneous_lambda = dynamic_cast<const TypedDagNode<RbVector<double> >*>(inspeciation);
    homogeneous_lambda = dynamic_cast<const TypedDagNode<double >*>(inspeciation);

    addParameter( homogeneous_lambda );
    addParameter( heterogeneous_lambda );

    if ( heterogeneous_lambda == NULL && homogeneous_lambda == NULL)
    {
        throw(RbException("Speciation rate must be of type RealPos or RealPos[]"));
    }
    else if( heterogeneous_lambda != NULL )
    {
        if( homogeneous_timeline == NULL && lambda_timeline == NULL ) throw( RbException("No time intervals provided for piecewise constant speciation rates") );

        if ( lambda_timeline != NULL && heterogeneous_lambda->getValue().size() != lambda_timeline->getValue().size() + 1 )
        {
            std::stringstream ss;
            ss << "Number of speciation rates (" << heterogeneous_lambda->getValue().size() << ") does not match number of time intervals (" << lambda_timeline->getValue().size() + 1 << ")";
            throw(RbException(ss.str()));
        }
        else if (homogeneous_timeline != NULL && heterogeneous_lambda->getValue().size() != homogeneous_timeline->getValue().size() + 1)
        {
            std::stringstream ss;
            ss << "Number of speciation rates (" << heterogeneous_lambda->getValue().size() << ") does not match number of time intervals (" << homogeneous_timeline->getValue().size() + 1 << ")";
            throw(RbException(ss.str()));
        }

        if (lambda_timeline != NULL)
        {
            std::vector<double> times = lambda_timeline->getValue();
            std::vector<double> times_sorted_ascending = times;
            std::vector<double> times_sorted_descending = times;

            sort(times_sorted_ascending.begin(), times_sorted_ascending.end() );
            sort(times_sorted_descending.rbegin(), times_sorted_descending.rend() );

            if( times == times_sorted_ascending )
            {
                ascending = true;
            }
            else if ( times != times_sorted_ascending )
            {
                throw(RbException("Speciation rate change times must be provided in order"));
            }
        }
    }
    /////////////

    heterogeneous_mu = dynamic_cast<const TypedDagNode<RbVector<double> >*>(inextinction);
    homogeneous_mu = dynamic_cast<const TypedDagNode<double >*>(inextinction);

    addParameter( homogeneous_mu );
    addParameter( heterogeneous_mu );

    if ( heterogeneous_mu == NULL && homogeneous_mu == NULL)
    {
        throw(RbException("Extinction rate must be of type RealPos or RealPos[]"));
    }
    else if( heterogeneous_mu != NULL )
    {
        if( homogeneous_timeline == NULL && mu_timeline == NULL ) throw( RbException("No time intervals provided for piecewise constant extinction rates") );

        if ( mu_timeline != NULL && heterogeneous_mu->getValue().size() != mu_timeline->getValue().size() + 1 )
        {
            std::stringstream ss;
            ss << "Number of extinction rates (" << heterogeneous_mu->getValue().size() << ") does not match number of time intervals (" << mu_timeline->getValue().size() + 1 << ")";
            throw(RbException(ss.str()));
        }
        else if (homogeneous_timeline != NULL && heterogeneous_mu->getValue().size() != homogeneous_timeline->getValue().size() + 1)
        {
            std::stringstream ss;
            ss << "Number of extinction rates (" << heterogeneous_mu->getValue().size() << ") does not match number of time intervals (" << homogeneous_timeline->getValue().size() + 1 << ")";
            throw(RbException(ss.str()));
        }

        if (mu_timeline != NULL)
        {
            std::vector<double> times = mu_timeline->getValue();
            std::vector<double> times_sorted = times;

            if( ascending )
            {
                sort(times_sorted.begin(), times_sorted.end() );
            }
            else
            {
                sort(times_sorted.rbegin(), times_sorted.rend() );
            }

            if ( times != times_sorted )
            {
                throw(RbException("Extinction rate change times must be provided in order"));
            }
        }
    }
   /////////

    heterogeneous_psi = dynamic_cast<const TypedDagNode<RbVector<double> >*>(inpsi);
    homogeneous_psi = dynamic_cast<const TypedDagNode<double >*>(inpsi);

    addParameter( homogeneous_psi );
    addParameter( heterogeneous_psi );

    if ( heterogeneous_psi == NULL && homogeneous_psi == NULL)
    {
        throw(RbException("Fossil sampling rate must be of type RealPos or RealPos[]"));
    }
    else if( heterogeneous_psi != NULL )
    {
        if( homogeneous_timeline == NULL && psi_timeline == NULL ) throw( RbException("No time intervals provided for piecewise constant serial sampling rates") );

        if ( psi_timeline != NULL && heterogeneous_psi->getValue().size() != psi_timeline->getValue().size() + 1 )
        {
            std::stringstream ss;
            ss << "Number of fossil sampling rates (" << heterogeneous_psi->getValue().size() << ") does not match number of time intervals (" << psi_timeline->getValue().size() + 1 << ")";
            throw(RbException(ss.str()));
        }
        else if (homogeneous_timeline != NULL && heterogeneous_psi->getValue().size() != homogeneous_timeline->getValue().size() + 1)
        {
            std::stringstream ss;
            ss << "Number of fossil sampling rates (" << heterogeneous_psi->getValue().size() << ") does not match number of time intervals (" << homogeneous_timeline->getValue().size() + 1 << ")";
            throw(RbException(ss.str()));
        }

        if (psi_timeline != NULL)
        {
            std::vector<double> times = psi_timeline->getValue();
            std::vector<double> times_sorted = times;

            if( ascending )
            {
                sort(times_sorted.begin(), times_sorted.end() );
            }
            else
            {
                sort(times_sorted.rbegin(), times_sorted.rend() );
            }

            if ( times != times_sorted )
            {
                throw(RbException("Fossil sampling rate change times must be provided in order"));
            }
        }
    }
 ///////

    heterogeneous_omega = dynamic_cast<const TypedDagNode<RbVector<double> >*>(inomega);
    homogeneous_omega = dynamic_cast<const TypedDagNode<double >*>(inomega);

    addParameter( homogeneous_omega );
    addParameter( heterogeneous_omega );

    if ( heterogeneous_omega == NULL && homogeneous_omega == NULL)
    {
        throw(RbException("Occurrence sampling rate must be of type RealPos or RealPos[]"));
    }
    else if( heterogeneous_omega != NULL )
    {
        if( homogeneous_timeline == NULL && omega_timeline == NULL ) throw( RbException("No time intervals provided for piecewise constant occurrence sampling") );

        if ( omega_timeline != NULL && heterogeneous_omega->getValue().size() != omega_timeline->getValue().size() + 1 )
        {
            std::stringstream ss;
            ss << "Number of occurrence sampling rates (" << heterogeneous_omega->getValue().size() << ") does not match number of time intervals (" << omega_timeline->getValue().size() + 1 << ")";
            throw(RbException(ss.str()));
        }
        else if (homogeneous_timeline != NULL && heterogeneous_omega->getValue().size() != homogeneous_timeline->getValue().size() + 1)
        {
            std::stringstream ss;
            ss << "Number of occurrence sampling rates (" << heterogeneous_omega->getValue().size() << ") does not match number of time intervals (" << homogeneous_timeline->getValue().size() + 1 << ")";
            throw(RbException(ss.str()));
        }

        if (omega_timeline != NULL)
        {
            std::vector<double> times = omega_timeline->getValue();
            std::vector<double> times_sorted = times;

            if( ascending )
            {
                sort(times_sorted.begin(), times_sorted.end() );
            }
            else
            {
                sort(times_sorted.rbegin(), times_sorted.rend() );
            }

            if ( times != times_sorted )
            {
                throw(RbException("Occurrence sampling rates change times must be provided in order"));
            }
        }
    }

    ////////
        heterogeneous_removalPr = dynamic_cast<const TypedDagNode<RbVector<double> >*>(inspeciation);
        homogeneous_removalPr = dynamic_cast<const TypedDagNode<double >*>(inspeciation);

        addParameter( homogeneous_removalPr );
        addParameter( heterogeneous_removalPr );

        if ( heterogeneous_removalPr == NULL && homogeneous_removalPr == NULL)
        {
            throw(RbException("Removal probability must be of type RealPos or RealPos[]"));
        }
        else if( heterogeneous_removalPr != NULL )
        {
            if( homogeneous_timeline == NULL && removalPr_timeline == NULL ) throw( RbException("No time intervals provided for piecewise constant removal probabilities") );

            if ( removalPr_timeline != NULL && heterogeneous_removalPr->getValue().size() != removalPr_timeline->getValue().size() + 1 )
            {
                std::stringstream ss;
                ss << "Number of removal probabilities (" << heterogeneous_removalPr->getValue().size() << ") does not match number of time intervals (" << removalPr_timeline->getValue().size() + 1 << ")";
                throw(RbException(ss.str()));
            }
            else if (homogeneous_timeline != NULL && heterogeneous_removalPr->getValue().size() != homogeneous_timeline->getValue().size() + 1)
            {
                std::stringstream ss;
                ss << "Number of removal probabilities (" << heterogeneous_removalPr->getValue().size() << ") does not match number of time intervals (" << homogeneous_timeline->getValue().size() + 1 << ")";
                throw(RbException(ss.str()));
            }

            if (removalPr_timeline != NULL)
            {
                std::vector<double> times = removalPr_timeline->getValue();
                std::vector<double> times_sorted_ascending = times;
                std::vector<double> times_sorted_descending = times;

                sort(times_sorted_ascending.begin(), times_sorted_ascending.end() );
                sort(times_sorted_descending.rbegin(), times_sorted_descending.rend() );

                if( times == times_sorted_ascending )
                {
                    ascending = true;
                }
                else if ( times != times_sorted_ascending )
                {
                    throw(RbException("Speciation rate change times must be provided in order"));
                }
            }
        }
    /////////////
    //prepareProbComputation();

    if (t != NULL)
    {
      delete value;
      value = &(t->getValue());
    }
    else
    {
      std::cout << "Simulate tree" <<std::endl;
      simulateTree();
      std::cout << "Simulate tree" <<std::endl;

    }
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
PiecewiseConstantOccurrenceBirthDeathProcess* PiecewiseConstantOccurrenceBirthDeathProcess::clone( void ) const
{
    return new PiecewiseConstantOccurrenceBirthDeathProcess( *this );
}


/**
 * Compute the log-transformed probability of the current value under the current parameter values.
 *
 */
double PiecewiseConstantOccurrenceBirthDeathProcess::computeLnProbabilityDivergenceTimes( void ) const
{
    // prepare the probability computation
    std::cout << "PreparProbComputation" << std::endl;
    prepareProbComputation();

    // variable declarations and initialization
    const RevBayesCore::Tree tree(*value);
    std::vector<double> occAges = std::vector<double>();
    const std::vector<double> time_points_Mt( 1, 0.0 );
    std::cout << "Compute likelihood Mt piecewise" << std::endl;
    MatrixReal B_Mt = RevBayesCore::ComputeLikelihoodsForwardsMtPiecewise(start_age, timeline, lambda, mu, psi, omega, homogeneous_rho, removalPr, maxHiddenLin, cond, time_points_Mt, useOrigin, occAges, tree);

    const long N = maxHiddenLin->getValue();
    const size_t k = value->getNumberOfExtantTips();

    double likelihood = B_Mt[0][0];
    // std::cout << "B_Mt[0][0] : " << B_Mt[0][0] << std::endl;
    for(int i = 1; i < N+1; i++){
        // std::cout << "B_Mt[0][" << i << "] : " << B_Mt[0][i] << std::endl;
        // std::cout << "B_Mt[0][" << i << "] * pow(rh,k) * pow(1.0 - rh,i) : " << B_Mt[0][i] * pow(rh,k) * pow(1.0 - rh,i) << std::endl;
        likelihood += B_Mt[0][i] * pow(rho,k) * pow(1.0 - rho,i);
    }
    double lnProbTimes = likelihood;
    return lnProbTimes;
}

/**
 * Compute the log probability of the current value under the current parameter values.
 * Tip-dating (Theorem 1, Stadler et al. 2013 PNAS)
 */
double PiecewiseConstantOccurrenceBirthDeathProcess::computeLnProbabilityTimes( void ) const
{
    // variable declarations and initialization
    double lnProbTimes = 0.0;
    double process_time = getOriginAge();
    size_t num_initial_lineages = 0;
    TopologyNode* root = &value->getRoot();

    if ( use_origin == true )
    {
        // If we are conditioning on survival from the origin,
        // then we must divide by 2 the log survival probability computed by AbstractBirthDeathProcess
        // TODO: Generalize AbstractBirthDeathProcess to allow conditioning on the origin
        num_initial_lineages = 1;
    }
    // if conditioning on root, root node must be a "true" bifurcation event
    else
    {
        if (root->getChild(0).isSampledAncestor() || root->getChild(1).isSampledAncestor())
        {
            return RbConstants::Double::neginf;
        }

        num_initial_lineages = 2;
    }

    // get node/time variables
    size_t num_nodes = value->getNumberOfNodes();

    // classify nodes
    int num_extant_taxa = 0;

    // retrieved the speciation times
    std::vector<double> serial_tip_ages = std::vector<double>();
    std::vector<double> sampled_ancestor_ages = std::vector<double>();
    std::vector<double> internal_node_ages = std::vector<double>();

    for (size_t i = 0; i < num_nodes; i++)
    {
        const TopologyNode& n = value->getNode( i );

        if ( n.isTip() && n.isFossil() && n.isSampledAncestor() )
        {
            // node is sampled ancestor
            sampled_ancestor_ages.push_back(n.getAge());
        }
        else if ( n.isTip() && n.isFossil() && !n.isSampledAncestor() )
        {
            // node is serial leaf
            serial_tip_ages.push_back( n.getAge() );
        }
        else if ( n.isTip() && !n.isFossil() )
        {
            // node is extant leaf
            num_extant_taxa++;
        }
        else if ( n.isInternal() && !n.getChild(0).isSampledAncestor() && !n.getChild(1).isSampledAncestor() )
        {
            if ( n.isRoot() == false || use_origin == true )
            {
                // node is bifurcation event (a "true" node)
                internal_node_ages.push_back( n.getAge() );
            }
        }
    }

    // add the serial tip age terms
    for (size_t i = 0; i < serial_tip_ages.size(); ++i)
    {
        if ( RbMath::isFinite(lnProbTimes) == false )
        {
            return RbConstants::Double::nan;
        }

        double t = serial_tip_ages[i];
        size_t index = l(t);

        // add the log probability for the serial sampling events
        if ( psi[index - 1] == 0.0 )
        {
            return RbConstants::Double::neginf;
            //std::stringstream ss;
            //ss << "The serial sampling rate in interval " << index << " is zero, but the tree has serial sampled tips in this interval.";
            //throw RbException(ss.str());
        }
        else
        {
            lnProbTimes += log( psi[index - 1] ) + log( p(index, t) / q(index, t) );
        }
    }

        // add the extant tip age term
    if (num_extant_taxa > 0)
    {
        lnProbTimes += num_extant_taxa * log( homogeneous_rho->getValue() );
    }

    // add the sampled ancestor age terms
    for (size_t i = 0; i < sampled_ancestor_ages.size(); ++i)
    {
        if ( RbMath::isFinite(lnProbTimes) == false )
        {
            return RbConstants::Double::nan;
        }

        double t = sampled_ancestor_ages[i];
        size_t index = l(t);
        lnProbTimes += log( psi[index - 1] );
    }

    // add the single lineage propagation terms
    for (size_t i = 1; i <= timeline.size(); ++i)
    {
        if ( RbMath::isFinite(lnProbTimes) == false )
        {
            return RbConstants::Double::nan;
        }

        double t = timeline[i-1];
        int div = survivors(t);
        lnProbTimes += div * log( q(i+1, t) );
    }

    // add the bifurcation age terms
    for (size_t i = 0; i < internal_node_ages.size(); ++i)
    {
        if ( RbMath::isFinite(lnProbTimes) == false )
        {
            return RbConstants::Double::nan;
        }

        double t = internal_node_ages[i];
        size_t index = l(t);
        lnProbTimes += log( q(index, t) * lambda[index - 1] );
    }

    // add the initial age term
    size_t index = l(process_time);
    lnProbTimes += num_initial_lineages * log( q(index, process_time ) );

    // condition on survival
    if ( condition == "survival" )
    {
        lnProbTimes -= num_initial_lineages * log( pSurvival(process_time,0) );
    }
    // condition on nTaxa
    else if ( condition == "nTaxa" )
    {
        lnProbTimes -= lnProbNumTaxa( value->getNumberOfTips(), 0, process_time, true );
    }

    if ( RbMath::isFinite(lnProbTimes) == false )
    {
        return RbConstants::Double::nan;
    }


    return lnProbTimes;
}


/**
 * return the index i so that t_{i-1} > t >= t_i
 * where t_i is the instantaneous sampling time (i = 0,...,l)
 * t_0 is origin
 * t_l = 0.0
 */
size_t PiecewiseConstantOccurrenceBirthDeathProcess::l(double t) const
{
    return timeline.rend() - std::upper_bound( timeline.rbegin(), timeline.rend(), t) + 1;
}


double PiecewiseConstantOccurrenceBirthDeathProcess::lnProbTreeShape(void) const
{
    // the birth death divergence times density is derived for a (ranked) unlabeled oriented tree
    // so we convert to a (ranked) labeled non-oriented tree probability by multiplying by 2^{n+m-1} / n!
    // where n is the number of extant tips, m is the number of extinct tips

    int num_taxa = (int)value->getNumberOfTips();
    int num_extinct = (int)value->getNumberOfExtinctTips();
    int num_sa = (int)value->getNumberOfSampledAncestors();

    return (num_taxa - num_sa - 1) * RbConstants::LN2 - RbMath::lnFactorial(num_taxa - num_extinct);
}


/**
 * Compute the probability of survival if the process starts with one species at time start and ends at time end.
 *
 * \param[in]    start      Start time of the process.
 * \param[in]    end        End/stopping time of the process.
 *
 * \return Probability of survival.
 */
double PiecewiseConstantOccurrenceBirthDeathProcess::pSurvival(double start, double end) const
{
    double t = start;

    std::vector<double> serial_bak = psi;

    std::fill(psi.begin(), psi.end(), 0.0);

    double p0 = p(l(t), t);

    psi = serial_bak;

    return 1.0 - p0;
}


/**
 * p_i(t)
 */
double PiecewiseConstantOccurrenceBirthDeathProcess::p( size_t i, double t ) const
{
    if ( t == 0)
        return 1.0;

    // get the parameters
    double b = lambda[i-1];
    double d = mu[i-1];
    double f = psi[i-1];
    double ti = timeline[i-1];
    double r = rho ;

    double diff = b - d - f;
    double bp   = b*f;
    double dt   = t - ti;

    double A = sqrt( diff*diff + 4.0*bp);
    double B = ( (1.0 - 2.0*(1.0-r)*p(i+1,ti) )*b + d + f ) / A;

    double e = exp(-A*dt);
    double tmp = b + d + f - A * ((1.0+B)-e*(1.0-B))/((1.0+B)+e*(1.0-B));

    return tmp / (2.0*b);
}


/**
 *
 *
 */
void PiecewiseConstantOccurrenceBirthDeathProcess::prepareProbComputation( void ) const
{
    // clean all the sets
    lambda.clear();
    mu.clear();
    psi.clear();
    omega.clear();
    removalPr.clear();


    if (homogeneous_timeline != NULL)
    {
        timeline = homogeneous_timeline->getValue();

        //timeline.push_back(0.0);

        for (size_t i = 0; i < timeline.size() + 1; i++)
        {
            lambda.push_back( getSpeciationRate(i) );

            mu.push_back( getExtinctionRate(i) );

            psi.push_back( getFossilSamplingRate(i) );

            omega.push_back( getOccurrenceSamplingRate(i) );

            removalPr.push_back( getRemovalProbability(i) );
        }
    }
    else
    {
        timeline.clear();

        std::set<double> eventTimes;

        if ( lambda_timeline != NULL )
        {
            const std::vector<double>& birthTimes = lambda_timeline->getValue();
            for (std::vector<double>::const_iterator it = birthTimes.begin(); it != birthTimes.end(); ++it)
            {
                eventTimes.insert( *it );
            }
        }

        if ( mu_timeline != NULL )
        {
            const std::vector<double>& deathTimes = mu_timeline->getValue();
            for (std::vector<double>::const_iterator it = deathTimes.begin(); it != deathTimes.end(); ++it)
            {
                eventTimes.insert( *it );
            }
        }

        if ( psi_timeline != NULL )
        {
            const std::vector<double>& fossilTimes = psi_timeline->getValue();
            for (std::vector<double>::const_iterator it = fossilTimes.begin(); it != fossilTimes.end(); ++it)
            {
                eventTimes.insert( *it );
            }
        }

        if ( omega_timeline != NULL )
        {
            const std::vector<double>& occurrenceTimes = omega_timeline->getValue();
            for (std::vector<double>::const_iterator it = occurrenceTimes.begin(); it != occurrenceTimes.end(); ++it)
            {
                eventTimes.insert( *it );
            }
        }

        if ( removalPr_timeline != NULL )
        {
            const std::vector<double>& removalPrTimes = removalPr_timeline->getValue();
            for (std::vector<double>::const_iterator it = removalPrTimes.begin(); it != removalPrTimes.end(); ++it)
            {
                eventTimes.insert( *it );
            }
        }



        for (std::set<double>::reverse_iterator it = eventTimes.rbegin(); it != eventTimes.rend(); it++)
        {
            timeline.push_back(*it);
        }
        //timeline.push_back(0.0);


        size_t birth_index = 0;
        size_t death_index = 0;
        size_t fossil_index = 0;
        size_t occurrence_index = 0;
        size_t removalPr_index = 0;

        for (size_t i = 0; i < timeline.size()+1; i++)
        {
            if ( lambda_timeline != NULL )
            {
                size_t index = ascending ? heterogeneous_lambda->getValue().size() - 1 - birth_index : birth_index;

                lambda.push_back( heterogeneous_lambda->getValue()[index] );

                if ( birth_index < lambda_timeline->getValue().size() && lambda_timeline->getValue()[index] == timeline[i] )
                {
                    birth_index++;
                }
            }
            else
            {
                lambda.push_back( homogeneous_lambda->getValue() );
            }

            if ( mu_timeline != NULL )
            {
                size_t index = ascending ? heterogeneous_mu->getValue().size() - 1 - death_index : death_index;

                mu.push_back( heterogeneous_mu->getValue()[index] );

                if ( death_index < mu_timeline->getValue().size() && mu_timeline->getValue()[index] == timeline[i] )
                {
                    death_index++;
                }
            }
            else
            {
                mu.push_back( homogeneous_mu->getValue() );
            }

            if ( psi_timeline != NULL )
            {
                size_t index = ascending ? heterogeneous_psi->getValue().size() - 1 - fossil_index : fossil_index;

                psi.push_back( heterogeneous_psi->getValue()[index] );

                if ( fossil_index < psi_timeline->getValue().size() && psi_timeline->getValue()[index] == timeline[i] )
                {
                    fossil_index++;
                }
            }
            else
            {
                psi.push_back( homogeneous_psi->getValue() );
            }

            if ( omega_timeline != NULL )
            {
                size_t index = ascending ? heterogeneous_omega->getValue().size() - 1 - occurrence_index : occurrence_index;

                omega.push_back( heterogeneous_omega->getValue()[index] );

                if ( occurrence_index < omega_timeline->getValue().size() && omega_timeline->getValue()[index] == timeline[i] )
                {
                    occurrence_index++;
                }
            }
            else
            {
                omega.push_back( homogeneous_omega->getValue() );
            }

            if ( removalPr_timeline != NULL )
            {
                size_t index = ascending ? heterogeneous_removalPr->getValue().size() - 1 - removalPr_index : removalPr_index;

                removalPr.push_back( heterogeneous_removalPr->getValue()[index] );

                if ( removalPr_index < removalPr_timeline->getValue().size() && removalPr_timeline->getValue()[index] == timeline[i] )
                {
                    removalPr_index++;
                }
            }
            else
            {
                removalPr.push_back( homogeneous_removalPr->getValue() );
            }
        }
    }
}


/**
 * q_i(t)
 */
double PiecewiseConstantOccurrenceBirthDeathProcess::q( size_t i, double t ) const
{

    if ( t == 0.0 )
    {
        return 1.0;
    }

    // get the parameters
    double b = lambda[i-1];
    double d = mu[i-1];
    double f = psi[i-1];
    double r = rho;
    double ti = timeline[i-1];

    double diff = b - d - f;
    double bp   = b*f;
    double dt   = t - ti;

    double A = sqrt( diff*diff + 4.0*bp);
    double B = ( (1.0 - 2.0*(1.0-r)*p(i+1,ti) )*b + d + f ) / A;

    double e = exp(-A*dt);
    double tmp = (1.0+B) + e*(1.0-B);

    return 4.0*e / (tmp*tmp);
}






double PiecewiseConstantOccurrenceBirthDeathProcess::getSpeciationRate( size_t index ) const
{

    // remove the old parameter first
    if ( homogeneous_lambda != NULL )
    {
        return homogeneous_lambda->getValue();
    }
    else
    {
        size_t num = heterogeneous_lambda->getValue().size();

        if (index >= num)
        {
            throw(RbException("Speciation rate index out of bounds"));
        }
        return ascending ? heterogeneous_lambda->getValue()[num - 1 - index] : heterogeneous_lambda->getValue()[index];
    }
}

double PiecewiseConstantOccurrenceBirthDeathProcess::getExtinctionRate( size_t index ) const
{

    // remove the old parameter first
    if ( homogeneous_mu != NULL )
    {
        return homogeneous_mu->getValue();
    }
    else
    {
        size_t num = heterogeneous_mu->getValue().size();

        if (index >= num)
        {
            throw(RbException("Extinction rate index out of bounds"));
        }
        return ascending ? heterogeneous_mu->getValue()[num - 1 - index] : heterogeneous_mu->getValue()[index];
    }
}

double PiecewiseConstantOccurrenceBirthDeathProcess::getFossilSamplingRate( size_t index ) const
{

    // remove the old parameter first
    if ( homogeneous_psi != NULL )
    {
        return homogeneous_psi->getValue();
    }
    else
    {
        size_t num = heterogeneous_psi->getValue().size();

        if (index >= num)
        {
            throw(RbException("Extinction rate index out of bounds"));
        }
        return ascending ? heterogeneous_psi->getValue()[num - 1 - index] : heterogeneous_psi->getValue()[index];
    }
}

double PiecewiseConstantOccurrenceBirthDeathProcess::getOccurrenceSamplingRate( size_t index ) const
{

    // remove the old parameter first
    if ( homogeneous_omega != NULL )
    {
        return homogeneous_omega->getValue();
    }
    else
    {
        size_t num = heterogeneous_omega->getValue().size();

        if (index >= num)
        {
            throw(RbException("Extinction rate index out of bounds"));
        }
        return ascending ? heterogeneous_omega->getValue()[num - 1 - index] : heterogeneous_omega->getValue()[index];
    }
}

double PiecewiseConstantOccurrenceBirthDeathProcess::getRemovalProbability( size_t index ) const
{

    // remove the old parameter first
    if ( homogeneous_removalPr != NULL )
    {
        return homogeneous_removalPr->getValue();
    }
    else
    {
        size_t num = heterogeneous_removalPr->getValue().size();

        if (index >= num)
        {
            throw(RbException("Extinction rate index out of bounds"));
        }
        return ascending ? heterogeneous_removalPr->getValue()[num - 1 - index] : heterogeneous_removalPr->getValue()[index];
    }
}





/**
 * Simulate new speciation times.
 */
double PiecewiseConstantOccurrenceBirthDeathProcess::simulateDivergenceTime(double origin, double present) const
{
    // incorrect placeholder for constant SSBDP


    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;

    size_t i = l(present);

    // get the parameters
    double age = origin - present;
    double b = lambda[i];
    double d = mu[i];
    double r = rho;


    // get a random draw
    double u = rng->uniform01();

    // compute the time for this draw
    // see Hartmann et al. 2010 and Stadler 2011
    double t = 0.0;
    if ( b > d )
    {
        if( r > 0.0 )
        {
            t = ( log( ( (b-d) / (1 - (u)*(1-((b-d)*exp((d-b)*age))/(r*b+(b*(1-r)-d)*exp((d-b)*age) ) ) ) - (b*(1-r)-d) ) / (r * b) ) )  /  (b-d);
        }
        else
        {
            t = log( 1 - u * (exp(age*(d-b)) - 1) / exp(age*(d-b)) ) / (b-d);
        }
    }
    else
    {
        if( r > 0.0 )
        {
            t = ( log( ( (b-d) / (1 - (u)*(1-(b-d)/(r*b*exp((b-d)*age)+(b*(1-r)-d) ) ) ) - (b*(1-r)-d) ) / (r * b) ) )  /  (b-d);
        }
        else
        {
            t = log( 1 - u * (1 - exp(age*(b-d)))  ) / (b-d);
        }
    }

    return present + t;
}


/**
 * Compute the diversity of the tree at time t.
 *
 * \param[in]    t      time at which we want to know the diversity.
 *
 * \return The diversity (number of species in the reconstructed tree).
 */
int PiecewiseConstantOccurrenceBirthDeathProcess::survivors(double t) const
{

    const std::vector<TopologyNode*>& nodes = value->getNodes();

    int survivors = 0;
    for (std::vector<TopologyNode*>::const_iterator it = nodes.begin(); it != nodes.end(); ++it)
    {
        TopologyNode* n = *it;
        if ( n->getAge() < t )
        {
            if ( n->isRoot() || n->getParent().getAge() > t )
            {
                survivors++;
            }
        }
    }

    return survivors;
}

/**
 * Swap the parameters held by this distribution.
 *
 * \param[in]    oldP      Pointer to the old parameter.
 * \param[in]    newP      Pointer to the new parameter.
 */
void PiecewiseConstantOccurrenceBirthDeathProcess::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == heterogeneous_lambda)
    {
        heterogeneous_lambda = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == heterogeneous_mu)
    {
        heterogeneous_mu = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == heterogeneous_psi)
    {
        heterogeneous_psi = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == heterogeneous_omega)
    {
        heterogeneous_omega = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == heterogeneous_removalPr)
    {
        heterogeneous_removalPr = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == homogeneous_lambda)
    {
        homogeneous_lambda = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == homogeneous_mu)
    {
        homogeneous_mu = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == homogeneous_psi)
    {
        homogeneous_psi = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == homogeneous_rho)
    {
        homogeneous_rho = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == homogeneous_omega)
    {
        homogeneous_omega = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == homogeneous_removalPr)
    {
        homogeneous_removalPr = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == homogeneous_timeline)
    {
        homogeneous_timeline = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == lambda_timeline)
    {
        lambda_timeline = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == mu_timeline)
    {
        mu_timeline = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == psi_timeline)
    {
        psi_timeline = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == omega_timeline)
    {
        omega_timeline = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == removalPr_timeline)
    {
        removalPr_timeline = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else
    {
        // delegate the super-class
        AbstractBirthDeathProcess::swapParameterInternal(oldP, newP);
    }
}
