#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <ostream>
#include <string>
#include <type_traits>
#include <vector>

#include "AbstractBirthDeathProcess.h"
#include "AbstractFossilizedBirthDeathRangeProcess.h"
#include "DistributionExponential.h"
#include "FossilizedBirthDeathSpeciationProcess.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbMathCombinatorialFunctions.h"
#include "RbException.h"
#include "RbMathFunctions.h"
#include "StringUtilities.h"
#include "StochasticNode.h"
#include "Taxon.h"
#include "TimeInterval.h"
#include "TopologyNode.h"
#include "Tree.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { template <class valueType> class RbVector; }

using namespace RevBayesCore;

/**
 * Constructor. 
 * We delegate most parameters to the base class and initialize the members.
 *
 * \param[in]    s              Speciation rates.
 * \param[in]    e              Extinction rates.
 * \param[in]    p              Fossil sampling rates.
 * \param[in]    c              Fossil observation counts.
 * \param[in]    r              Instantaneous sampling probabilities.
 * \param[in]    a              Anagenetic speciation rates.
 * \param[in]    b              Symmetric speciation probability.
 * \param[in]    t              Rate change times.
 * \param[in]    cdt            Condition of the process (time/sampling/survival).
 * \param[in]    tn             Taxa.
 * \param[in]    c              Complete sampling?
 * \param[in]    re             Augmented age resampling weight.
 */
FossilizedBirthDeathSpeciationProcess::FossilizedBirthDeathSpeciationProcess(const TypedDagNode<double> *ra,
                                                           const DagNode *inspeciation,
                                                           const DagNode *inextinction,
                                                           const DagNode *inpsi,
                                                           const TypedDagNode<double> *inrho,
                                                           const DagNode *inlambda_a,
                                                           const DagNode *inbeta,
                                                           const TypedDagNode< RbVector<double> > *intimes,
                                                           const std::string &incondition,
                                                           const std::vector<Taxon> &intaxa,
                                                           bool c,
                                                           bool re) :
    AbstractBirthDeathProcess(ra, incondition, intaxa, true, NULL),
    AbstractFossilizedBirthDeathRangeProcess(inspeciation, inextinction, inpsi, inrho, intimes, incondition, intaxa, c, re)
{
    for(std::vector<const DagNode*>::iterator it = range_parameters.begin(); it != range_parameters.end(); it++)
    {
        addParameter(*it);
    }

    homogeneous_lambda_a             = NULL;
    homogeneous_beta                 = NULL;
    heterogeneous_lambda_a           = NULL;
    heterogeneous_beta               = NULL;

    heterogeneous_lambda_a = dynamic_cast<const TypedDagNode<RbVector<double> >*>(inlambda_a);
    homogeneous_lambda_a   = dynamic_cast<const TypedDagNode<double >*>(inlambda_a);
    heterogeneous_beta     = dynamic_cast<const TypedDagNode<RbVector<double> >*>(inbeta);
    homogeneous_beta       = dynamic_cast<const TypedDagNode<double >*>(inbeta);

    addParameter( homogeneous_lambda_a );
    addParameter( heterogeneous_lambda_a );
    addParameter( homogeneous_beta );
    addParameter( heterogeneous_beta );

    RbException no_timeline_err = RbException("No time intervals provided for heterogeneous fossilized birth death process");

    RbException inconsistent_rates = RbException("Inconsistent number of rates in fossilized birth death process.");

    if ( heterogeneous_lambda_a != NULL )
    {
        if ( timeline == NULL ) throw no_timeline_err;

        if ( heterogeneous_lambda_a->getValue().size() != num_intervals ) throw inconsistent_rates;
    }
    if ( heterogeneous_beta != NULL )
    {
        if ( timeline == NULL ) throw no_timeline_err;

        if ( heterogeneous_beta->getValue().size() != num_intervals ) throw inconsistent_rates;
    }

    I             = std::vector<bool>(taxa.size(), false);

    anagenetic    = std::vector<double>(num_intervals, 0.0);
    symmetric     = std::vector<double>(num_intervals, 0.0);

    
    redrawValue();
    updateStartEndTimes(this->getValue().getRoot());
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself 
 */
FossilizedBirthDeathSpeciationProcess* FossilizedBirthDeathSpeciationProcess::clone( void ) const
{
    return new FossilizedBirthDeathSpeciationProcess( *this );
}


/**
 * Compute the log-transformed probability of the current value under the current parameter values.
 *
 */
double FossilizedBirthDeathSpeciationProcess::computeLnProbabilityDivergenceTimes( void )
{
    double lnProb = computeLnProbabilityRanges();

    lnProb += computeLnProbabilityTimes();

    return lnProb;
}


/**
 * Compute the log-transformed probability of the current value under the current parameter values.
 *
 */
double FossilizedBirthDeathSpeciationProcess::computeLnProbabilityTimes( void ) const
{
    double lnProb = 0.0;

    for (size_t i = 0; i < taxa.size(); i++)
    {
        // include the anagenetic speciation density for descendants of sampled ancestors
        if ( I[i] == true )
        {
            double y_a   = b_i[i];
            double d     = d_i[i];

            size_t y_ai = findIndex(y_a);
            size_t di   = findIndex(d);

            // offset speciation density
            lnProb -= log( birth[y_ai] );
            // offset the extinction density for the ancestor
            lnProb -= log( death[y_ai] );
            // include anagenetic speciation density
            lnProb += log( anagenetic[y_ai] );
        }
    }

    return lnProb;
}


double FossilizedBirthDeathSpeciationProcess::getMaxTaxonAge( const TopologyNode& node ) const
{
    if( node.isTip() )
    {
        return age[node.getIndex()];
    }
    else
    {
        double max = 0;

        for( size_t i = 0; i < node.getNumberOfChildren(); i++)
        {
            max = std::max( getMaxTaxonAge( node.getChild(i) ), max );
        }

        return max;
    }
}


double FossilizedBirthDeathSpeciationProcess::lnProbTreeShape(void) const
{
    // the fossilized birth death range divergence times density is derived for an unlabeled oriented tree
    // so we convert to a labeled oriented tree probability by multiplying by 1 / n!
    // where n is the number of extant tips

    return - RbMath::lnFactorial( value->getNumberOfExtantTips() );
}


/**
 * Compute the probability of survival if the process starts with one species at time start and ends at time end.
 *
 * \param[in]    start      Start time of the process.
 * \param[in]    end        End/stopping time of the process.
 *
 * \return Probability of survival.
 */
double FossilizedBirthDeathSpeciationProcess::pSurvival(double start, double end) const
{
    return AbstractFossilizedBirthDeathRangeProcess::p(findIndex(start), start, true);
}


/**
 * q_i(t)
 */
double FossilizedBirthDeathSpeciationProcess::q( size_t i, double t, bool tilde ) const
{

    if ( t == times[i] ) return 0.0;

    // get the parameters
    double b = birth[i];
    double d = death[i];
    double f = fossil[i];
    double r = (i == 0 ? homogeneous_rho->getValue() : 0.0);
    double ti = times[i];

    double diff = b - d - f;
    double dt   = t - ti;

    double A = sqrt( diff*diff + 4.0*b*f);
    double B = ( (1.0 - 2.0*(1.0-r)*p_i[i] )*b + d + f ) / A;

    double ln_e = -A*dt;

    double tmp = (1.0 + B) + exp(ln_e)*(1.0 - B);

    double q = log(4.0) + ln_e - 2.0*log(tmp);

    if (tilde)
    {
        q = 0.5 * (q - (b+d+f)*dt);

        double a = anagenetic[i];
        double s = symmetric[i];

        q = - a - s * (b + d + f) * dt + (1.0 - s) * q;
    }

    return q;
}


/**
 *
 */
void FossilizedBirthDeathSpeciationProcess::redrawValue(void)
{
    AbstractBirthDeathProcess::redrawValue();

    const std::vector<TopologyNode*> nodes = this->getValue().getNodes();

    for( size_t i = 0; i < this->getValue().getNumberOfTips(); i++)
    {
        size_t j = find(taxa.begin(), taxa.end(), nodes[i]->getTaxon()) - taxa.begin();

        nodes[i]->setIndex(j);
    }

    this->getValue().orderNodesByIndex();
}


/**
 *
 */
void FossilizedBirthDeathSpeciationProcess::simulateClade(std::vector<TopologyNode *> &n, double age, double present, bool alwaysReturn)
{

    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;

    // get the minimum birth age
    std::vector<double> first_occurrences;

    double current_age = RbConstants::Double::inf;
    double minimum_age = 0.0;
    double max_age = getOriginAge();

    for (size_t i = 0; i < n.size(); ++i)
    {
        // make sure the tip age is equal to the last occurrence
        if( n[i]->isTip() )
        {
            bool extinct = n[i]->getTaxon().isExtinct();

            double present = times.front();

            n[i]->setAge( extinct * rng->uniform01() * (n[i]->getTaxon().getMinAge() - present) + present );

            size_t j = find(taxa.begin(), taxa.end(), n[i]->getTaxon()) - taxa.begin();

            double minmax = std::max(o_i[i], n[i]->getAge());
            this->age[j] = GLOBAL_RNG->uniform01()*(std::min(max_age, taxa[i].getMaxAge()) - minmax) + minmax;
        }

        double first_occurrence = getMaxTaxonAge( *n[i] );

        if( first_occurrence > minimum_age )
        {
            minimum_age = first_occurrence;
        }

        first_occurrences.push_back( first_occurrence );

        if ( current_age > n[i]->getAge() )
        {
            current_age = n[i]->getAge();
        }
    }

    if( minimum_age > max_age )
    {
        throw RbException() << "Tree age is " << max_age << " but minimum fossil origin is " << minimum_age;
    }


    if ( age <= minimum_age )
    {
        age = rng->uniform01() * ( max_age - minimum_age ) + minimum_age;
    }


    std::vector<double> ages;
    while ( n.size() > 2 && current_age < age )
    {

        // get all the nodes with first occurrences younger than the current age
        std::vector<TopologyNode*> active_nodes;
        std::vector<TopologyNode*> active_right_nodes;
        for (size_t i = 0; i < n.size(); ++i)
        {
            if ( current_age >= n[i]->getAge() )
            {
                active_nodes.push_back( n[i] );

                if( current_age >= first_occurrences[i] )
                {
                    active_right_nodes.push_back( n[i] );
                }
            }

        }

        // we need to get the next node age older than the current age
        double next_node_age = age;
        for (size_t i = 0; i < n.size(); ++i)
        {
            if ( current_age < n[i]->getAge() && n[i]->getAge() < next_node_age )
            {
                next_node_age = n[i]->getAge();
            }
            if ( current_age < first_occurrences[i] && first_occurrences[i] < next_node_age )
            {
                next_node_age = first_occurrences[i];
            }

        }

        // only simulate if there are at least two valid/active nodes and one active right node
        if ( active_nodes.size() <= 2 || active_right_nodes.empty() )
        {
            current_age = next_node_age;
        }
        else
        {
            // now we simulate new ages
            double next_sim_age = simulateNextAge(active_nodes.size()-2, age, present, current_age, alwaysReturn);

            if ( next_sim_age < next_node_age )
            {
                // randomly pick two nodes
                size_t index_left = static_cast<size_t>( floor(rng->uniform01()*active_nodes.size()) );
                TopologyNode* left_child = active_nodes[index_left];

                size_t index_right = static_cast<size_t>( floor(rng->uniform01()*active_right_nodes.size()) );
                TopologyNode* right_child = active_right_nodes[index_right];

                while( left_child == right_child )
                {
                    index_left = static_cast<size_t>( floor(rng->uniform01()*active_nodes.size()) );
                    left_child = active_nodes[index_left];

                    index_right = static_cast<size_t>( floor(rng->uniform01()*active_right_nodes.size()) );
                    right_child = active_right_nodes[index_right];
                }

                // erase the nodes also from the origin nodes vector
                std::vector<TopologyNode *>::iterator child_it_left = std::find( n.begin(), n.end(), left_child );
                std::vector<double>::iterator fa_it_left = first_occurrences.begin() + std::distance( n.begin(), child_it_left );
                double fa_left = *fa_it_left;

                first_occurrences.erase(fa_it_left);
                n.erase(child_it_left);

                std::vector<TopologyNode *>::iterator child_it_right = std::find( n.begin(), n.end(), right_child );
                std::vector<double>::iterator fa_it_right = first_occurrences.begin() + std::distance( n.begin(), child_it_right );
                double fa_right = *fa_it_right;

                first_occurrences.erase(fa_it_right);
                n.erase(child_it_right);


                // create a parent for the two
                TopologyNode *parent = new TopologyNode();
                parent->addChild( left_child );
                parent->addChild( right_child );
                left_child->setParent( parent );
                right_child->setParent( parent );
                parent->setAge( next_sim_age );

                // insert the parent to our list
                n.push_back( parent );
                first_occurrences.push_back( std::max(fa_left, fa_right) );

                current_age = next_sim_age;
                ages.push_back( next_sim_age );
            }
            else
            {
                current_age = next_node_age;
            }

        }

        if ( n.size() > 2 && current_age >= age  )
            throw RbException() << "Unexpected number of taxa (remaining #taxa was " << n.size() << " and age was " << current_age << " with maximum age of " << age << ") in tree simulation";

    }


    if ( n.size() == 2 )
    {

        // pick two nodes
        TopologyNode* left_child = n[0];
        TopologyNode* right_child = n[1];

        // make sure the speciation event is older than the new species first occurrence
        if( first_occurrences[1] > age )
        {
            if( first_occurrences[0] > age )
            {
                throw RbException() << "Cannot simulate clade of age " << age << ", minimum age is " << minimum_age;
            }
            else
            {
                std::swap( left_child, right_child );
            }
        }
        else if( age > first_occurrences[0] )
        {
            if( rng->uniform01() < 0.5 )
            {
                std::swap( left_child, right_child );
            }
        }

        // erase the nodes also from the origin nodes vector
        n.clear();

        // create a parent for the two
        TopologyNode *parent = new TopologyNode();
        parent->addChild( left_child );
        parent->addChild( right_child );
        left_child->setParent( parent );
        right_child->setParent( parent );
        parent->setAge( age );

        // insert the parent to our list
        n.push_back( parent );
    }
    else
    {
        throw RbException() << "Unexpected number of taxa (" << n.size() << ") in tree simulation";
    }


}

/**
 * @param alwaysReturn whether the simulation can return times which are not valid draws from the distribution (for initial values)
*/
std::vector<double> FossilizedBirthDeathSpeciationProcess::simulateDivergenceTimes(size_t n, double origin, double present, double min, bool alwaysReturn) const
{
    if(!alwaysReturn) throw RbException("Impossible to simulate under a true FossilizedBirthDeathSpeciation process");

    std::vector<double> t(n, 0.0);

    for (size_t i = 0; i < n; ++i)
    {
        t[i] = simulateDivergenceTime(origin, min);
    }

    // finally sort the times
    std::sort(t.begin(), t.end());

    return t;
}

/**
 * Simulate new speciation times.
 */
double FossilizedBirthDeathSpeciationProcess::simulateDivergenceTime(double origin, double present) const
{
    // incorrect placeholder for constant SSBDP
    // direct forward simulation under FBDRP is not feasible

    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;

    size_t i = findIndex(present);

    // get the parameters
    double age = origin - present;
    double b = birth[i];
    double d = death[i];
    double p_e = i == 0 ? homogeneous_rho->getValue() : 0.0;


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


int FossilizedBirthDeathSpeciationProcess::updateStartEndTimes( const TopologyNode& node )
{
    if( node.isTip() )
    {
        return node.getIndex();
    }

    int species = -1;

    std::vector<TopologyNode* > children = node.getChildren();

    bool sa = node.isSampledAncestorTipOrParent();

    for(int c = 0; c < children.size(); c++)
    {
        const TopologyNode& child = *children[c];

        int i = updateStartEndTimes(child);

        // if child is a tip, set the species/end time
        if( child.isTip() )
        {
            double age = child.getAge();

            if ( age != d_i[i] )
            {
                d_i[i] = age;
                dirty_psi[i] = true;
                dirty_taxa[i] = true;
                // resample augmented age
                if ( resampling == true && resampled == false )
                {
                    resampleAge(i);
                }
            }
        }

        // is child a new species?
        // set start time at this node
        if( ( sa == false && c > 0 ) || ( sa && !child.isSampledAncestorTip() ) )
        {
            double age = node.getAge(); // y_{a(i)}

            if ( age != b_i[i] )
            {
                b_i[i] = age;
                dirty_psi[i] = true;
                dirty_taxa[i] = true;
                // resample augmented age
                if ( resampling == true && resampled == false )
                {
                    resampleAge(i);
                }
            }

            I[i] = sa;
        }
        // child is the ancestral species
        else
        {
            // propagate species index
            species = i;

            // if this is the root
            // set the start time to the origin
            if( node.isRoot() )
            {
                double age = getOriginAge();

                if ( age != b_i[i] )
                {
                    b_i[i] = age;
                    origin = age;
                    dirty_psi[i] = true;
                    dirty_taxa[i] = true;
                    // resample augmented age
                    if ( resampling == true && resampled == false )
                    {
                        resampleAge(i);
                    }
                }
            }
        }
    }

    return species;
}

/**
 *
 *
 */
void FossilizedBirthDeathSpeciationProcess::prepareProbComputation( void ) const 
{
    AbstractFossilizedBirthDeathRangeProcess::prepareProbComputation();

    if ( homogeneous_lambda_a != NULL )
    {
        anagenetic = std::vector<double>(num_intervals, homogeneous_lambda_a->getValue() );
    }
    else
    {
        anagenetic = heterogeneous_lambda_a->getValue();
    }
    if ( homogeneous_beta != NULL )
    {
        symmetric = std::vector<double>(num_intervals, homogeneous_beta->getValue() );
    }
    else
    {
        symmetric = heterogeneous_beta->getValue();
    }

    for (size_t i = 0; i < num_intervals; i++)
    {
        if ( i < num_intervals-1 )
        {
            double dt = times[i+1] - times[i];

            q_tilde_i[i] = - anagenetic[i] - symmetric[i] * (birth[i] + death[i] + fossil[i]) * dt + (1.0 - symmetric[i]) * q_tilde_i[i];
        }
    }
}


/**
 * Compute the log-transformed probability of the current value under the current parameter values.
 *
 */
void FossilizedBirthDeathSpeciationProcess::updateStartEndTimes( void )
{
    updateStartEndTimes(getValue().getRoot());
}


void FossilizedBirthDeathSpeciationProcess::keepSpecialization(const DagNode *toucher)
{
    AbstractFossilizedBirthDeathRangeProcess::keepSpecialization(toucher);
}


void FossilizedBirthDeathSpeciationProcess::restoreSpecialization(const DagNode *toucher)
{
    AbstractFossilizedBirthDeathRangeProcess::restoreSpecialization(toucher);
}


void FossilizedBirthDeathSpeciationProcess::touchSpecialization(const DagNode *toucher, bool touchAll)
{
    if ( toucher == dag_node )
    {
        if ( touched == false )
        {
            stored_likelihood = partial_likelihood;
            stored_Psi = Psi;

            std::set<size_t> touched_indices = dag_node->getTouchedElementIndices();

            for ( std::set<size_t>::iterator it = touched_indices.begin(); it != touched_indices.end(); it++)
            {
                size_t i = (*it) / taxa.size();

                dirty_psi[i]  = true;
                dirty_taxa[i] = true;
            }
        }

        touched = true;
    }
    else
    {
        AbstractBirthDeathProcess::touchSpecialization(toucher, touchAll);
        AbstractFossilizedBirthDeathRangeProcess::touchSpecialization(toucher, touchAll);
    }
}


/**
 * Swap the parameters held by this distribution.
 *
 * \param[in]    oldP      Pointer to the old parameter.
 * \param[in]    newP      Pointer to the new parameter.
 */
void FossilizedBirthDeathSpeciationProcess::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == heterogeneous_lambda_a)
    {
        heterogeneous_lambda_a = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == heterogeneous_beta)
    {
        heterogeneous_beta = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == homogeneous_lambda_a)
    {
        homogeneous_lambda_a = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == homogeneous_beta)
    {
        homogeneous_beta = static_cast<const TypedDagNode<double>* >( newP );
    }
    else
    {
        AbstractBirthDeathProcess::swapParameterInternal(oldP, newP);
        AbstractFossilizedBirthDeathRangeProcess::swapParameterInternal(oldP, newP);
    }
}
