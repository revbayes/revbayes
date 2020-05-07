#include <stddef.h>
#include <cmath>
#include <iomanip>
#include <ostream>
#include <vector>

#include "DistributionBeta.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RateAgeBetaShift.h"
#include "TopologyNode.h"
#include "AbstractMove.h"
#include "DagNode.h"
#include "RbException.h"
#include "RbOrderedSet.h"
#include "StochasticNode.h"
#include "Tree.h"

using namespace RevBayesCore;

RateAgeBetaShift::RateAgeBetaShift(StochasticNode<Tree> *tr, std::vector<StochasticNode<double> *> v, StochasticNode<RbVector<double> > *sv, double d, bool t, double w) : AbstractMove( w, t),
    tree( tr ),
    rates_vec( v ),
    rates( sv ),
    delta( d ),
    stored_node( NULL ),
    stored_age( 0.0 ),
    stored_rates( (rates == NULL ? rates_vec.size() : rates->getValue().size()), 0.0 ),
    num_accepted_current_period( 0 ),
    num_accepted_total( 0 )
{
    
    addNode( tree );
    addNode( rates );
    for (std::vector<StochasticNode<double>* >::iterator it = rates_vec.begin(); it != rates_vec.end(); ++it)
    {
        // get the pointer to the current node
        DagNode* the_node = *it;
        
        addNode( the_node );
    }

    
}


/**
 * Basic destructor doing nothing.
 */
RateAgeBetaShift::~RateAgeBetaShift( void )
{
    // nothing special to do
    // everything should be taken care of in the base class
    
}




/* Clone object */
RateAgeBetaShift* RateAgeBetaShift::clone( void ) const
{
    
    return new RateAgeBetaShift( *this );
}



const std::string& RateAgeBetaShift::getMoveName( void ) const
{
    
    static std::string name = "RateAgeBetaShift";
    
    return name;
}


double RateAgeBetaShift::getMoveTuningParameter( void ) const
{
    return delta;
}


size_t RateAgeBetaShift::getNumberAcceptedCurrentPeriod( void ) const
{
    return num_accepted_current_period;
}


size_t RateAgeBetaShift::getNumberAcceptedTotal( void ) const
{
    return num_accepted_total;
}


/** Perform the move */
void RateAgeBetaShift::performMcmcMove( double prHeat, double lHeat, double pHeat )
{
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    Tree& tau = tree->getValue();
    RbOrderedSet<DagNode*> affected;
    tree->initiateGetAffectedNodes( affected );
    
    double oldLnLike = 0.0;
    bool check_likelihood_shortcuts = rng->uniform01() < 0.001;
    if ( check_likelihood_shortcuts == true )
    {
        for (RbOrderedSet<DagNode*>::iterator it = affected.begin(); it != affected.end(); ++it)
        {
            (*it)->touch();
            oldLnLike += (*it)->getLnProbability();
        }
    }
    
    // pick a random node which is not the root and neithor the direct descendant of the root
    TopologyNode* node;
    size_t nide_idx = 0;
    do {
        double u = rng->uniform01();
        nide_idx = size_t( std::floor(tau.getNumberOfNodes() * u) );
        node = &tau.getNode(nide_idx);
    } while ( node->isRoot() || node->isTip() ); 
    
    TopologyNode& parent = node->getParent();
    
    // we need to work with the times
    double parent_age  = parent.getAge();
    double my_age      = node->getAge();
    double child_Age   = node->getChild( 0 ).getAge();
    if ( child_Age < node->getChild( 1 ).getAge())
    {
        child_Age = node->getChild( 1 ).getAge();
    }
    
    // now we store all necessary values
    stored_node = node;
    stored_age = my_age;
    
    
    stored_rates[nide_idx] = (rates == NULL ? rates_vec[nide_idx]->getValue() : rates->getValue()[nide_idx]);
    for (size_t i = 0; i < node->getNumberOfChildren(); i++)
    {
        size_t child_idx = node->getChild(i).getIndex();
        stored_rates[child_idx] = (rates == NULL ? rates_vec[child_idx]->getValue() : rates->getValue()[child_idx]);
    }
    
    
    // draw new ages and compute the hastings ratio at the same time
    double m = (my_age-child_Age) / (parent_age-child_Age);
    double a = delta * m + 1.0;
    double b = delta * (1.0-m) + 1.0;
    double new_m = RbStatistics::Beta::rv(a, b, *rng);
    double my_new_age = (parent_age-child_Age) * new_m + child_Age;
    
    // compute the Hastings ratio
    double forward = RbStatistics::Beta::lnPdf(a, b, new_m);
    double new_a = delta * new_m + 1.0;
    double new_b = delta * (1.0-new_m) + 1.0;
    double backward = RbStatistics::Beta::lnPdf(new_a, new_b, m);
    
    // set the age
    tau.getNode(nide_idx).setAge( my_new_age );
    
    // touch the tree so that the likelihoods are getting stored
    tree->touch();
    
    // get the probability ratio of the tree
    double tree_prob_ratio = tree->getLnProbabilityRatio();
    
    
    // set the rates
    double pa = node->getParent().getAge();
    double my_new_rate = (pa - my_age) * stored_rates[nide_idx] / (pa - my_new_age);
    
    // now we set the new value
    // this will automatically call a touch
    if ( rates == NULL )
    {
        rates_vec[nide_idx]->setValue( new double( my_new_rate ) );
    }
    else
    {
        rates->getValue()[nide_idx] = my_new_rate;
        rates->touch();
    }
    // get the probability ratio of the new rate
    double rates_prob_ratio = ( rates == NULL ? rates_vec[nide_idx]->getLnProbabilityRatio() : 0.0 );
    double jacobian = log((pa - my_age) / (pa - my_new_age));
    
    for (size_t i = 0; i < node->getNumberOfChildren(); i++)
    {
        size_t child_idx = node->getChild(i).getIndex();
        double a = node->getChild(i).getAge();
        double child_new_rate = (my_age - a) * stored_rates[child_idx] / (my_new_age - a);
        
        // now we set the new value
        // this will automatically call a touch
        if ( rates == NULL )
        {
            rates_vec[child_idx]->setValue( new double( child_new_rate ) );
        }
        else
        {
            rates->getValue()[child_idx] = child_new_rate;
        }

        // get the probability ratio of the new rate
        if ( rates == NULL )
        {
            rates_prob_ratio += rates_vec[child_idx]->getLnProbabilityRatio();
        }
        jacobian += log((my_age - a) / (my_new_age - a));
        
    }
    if ( rates != NULL )
    {
        rates_prob_ratio = rates->getLnProbabilityRatio();
    }
    
    if ( check_likelihood_shortcuts == true )
    {
        tree->touch();
        if ( rates == NULL )
        {
            rates_vec[nide_idx]->touch();
        }
        else
        {
            rates->touch();
        }
        for (size_t i = 0; i < node->getNumberOfChildren(); i++)
        {
            size_t childIdx = node->getChild(i).getIndex();
            if ( rates == NULL )
            {
                rates_vec[childIdx]->touch();
            }
        }
        
        double ln_prob_ratio = 0;
        double new_ln_like = 0;
        for (RbOrderedSet<DagNode*>::iterator it = affected.begin(); it != affected.end(); ++it)
        {

            double tmp = (*it)->getLnProbabilityRatio();
            ln_prob_ratio += tmp;
            new_ln_like += (*it)->getLnProbability();
        }
    
        if ( RbMath::isFinite(ln_prob_ratio) && fabs(ln_prob_ratio) > 1E-2 )
        {
            double ln_prob_ratio_2 = 0;
            double new_ln_like_2 = 0;
            for (RbOrderedSet<DagNode*>::iterator it = affected.begin(); it != affected.end(); ++it)
            {
                
                double tmp2 = (*it)->getLnProbabilityRatio();
                ln_prob_ratio_2 += tmp2;
                new_ln_like_2 += (*it)->getLnProbability();
                std::cerr << "\tP(" << (*it)->getName() << ") = " << tmp2 << std::endl;
            }
            
            throw RbException("Likelihood shortcut computation failed in rate-age-proposal.");
        }
        
    }
    
    double hastings_ratio = backward - forward + jacobian;
    double ln_acceptance_ratio = pHeat * prHeat * (tree_prob_ratio + rates_prob_ratio) + hastings_ratio;
    
    if (ln_acceptance_ratio >= 0.0)
    {
        num_accepted_total++;
        num_accepted_current_period++;
        
        tree->touch();
        tree->keep();
        if ( rates == NULL )
        {
            rates_vec[nide_idx]->touch();
            rates_vec[nide_idx]->keep();
        }
        else
        {
            rates->touch();
            rates->keep();
        }
        for (size_t i = 0; i < node->getNumberOfChildren(); i++)
        {
            size_t childIdx = node->getChild(i).getIndex();
            if ( rates == NULL )
            {
                rates_vec[childIdx]->touch();
                rates_vec[childIdx]->keep();
            }
        }
    }
    else if (ln_acceptance_ratio < -300.0)
    {
        reject();
        tree->restore();
        if ( rates == NULL )
        {
            rates_vec[nide_idx]->restore();
        }
        else
        {
            rates->restore();
        }
        for (size_t i = 0; i < node->getNumberOfChildren(); i++)
        {
            size_t childIdx = node->getChild(i).getIndex();
            if ( rates == NULL )
            {
                rates_vec[childIdx]->restore();
            }
        }
    }
    else
    {
        double r = exp(ln_acceptance_ratio);
        // Accept or reject the move
        double u = GLOBAL_RNG->uniform01();
        if (u < r)
        {
            num_accepted_total++;
            num_accepted_current_period++;
            
            //keep
            tree->touch();
            tree->keep();
            if ( rates == NULL )
            {
                rates_vec[nide_idx]->touch();
                rates_vec[nide_idx]->keep();
            }
            else
            {
                rates->touch();
                rates->keep();
            }
            for (size_t i = 0; i < node->getNumberOfChildren(); i++)
            {
                size_t childIdx = node->getChild(i).getIndex();
                if ( rates == NULL )
                {
                    rates_vec[childIdx]->touch();
                    rates_vec[childIdx]->keep();
                }
            }
        }
        else
        {
            reject();
            tree->restore();
            if ( rates == NULL )
            {
                rates_vec[nide_idx]->restore();
            }
            else
            {
                rates->restore();
            }
            for (size_t i = 0; i < node->getNumberOfChildren(); i++)
            {
                size_t childIdx = node->getChild(i).getIndex();
                if ( rates == NULL )
                {
                    rates_vec[childIdx]->restore();
                }
            }
        }
    }

}


void RateAgeBetaShift::printSummary(std::ostream &o, bool current_period) const
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
    const std::string &dn_name = (*nodes.begin())->getName();
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
    size_t num_accepted = num_accepted_total;
    if (current_period == true)
    {
        num_tried = num_tried_current_period;
        num_accepted = num_accepted_current_period;
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
    if (num_accepted > 0) a_length -= (int)log10(num_accepted);
    
    for (int i = 0; i < a_length; ++i)
    {
        o << " ";
    }
    o << num_accepted;
    o << " ";
    
    // print the acceptance ratio
    double ratio = num_accepted / (double)num_tried;
    if (num_tried == 0) ratio = 0;
    int r_length = 5;
    
    for (int i = 0; i < r_length; ++i)
    {
        o << " ";
    }
    o << ratio;
    o << " ";
    
//    proposal->printParameterSummary( o );
    o << "delta = " << delta;
    
    o << std::endl;
    
    o.setf(previousFlags);
    o.precision(previousPrecision);
    
}


void RateAgeBetaShift::reject( void )
{
    
    // undo the proposal
    tree->getValue().getNode( stored_node->getIndex() ).setAge( stored_age );
    
    // undo the rates
    size_t node_idx = stored_node->getIndex();
    if ( rates == NULL )
    {
        rates_vec[node_idx]->setValue(new double(stored_rates[node_idx]));
    }
    else
    {
        rates->getValue()[ node_idx ] = stored_rates[node_idx];
        rates->touch();
    }
    for (size_t i = 0; i < stored_node->getNumberOfChildren(); i++)
    {
        size_t childIdx = stored_node->getChild(i).getIndex();
        if ( rates == NULL )
        {
            rates_vec[childIdx]->setValue(new double(stored_rates[childIdx]));
        }
        else
        {
            rates->getValue()[ childIdx ] = stored_rates[childIdx];
            rates->touch();
        }
    }

    
#ifdef ASSERTIONS_TREE
    if ( fabs(storedAge - storedNode->getAge()) > 1E-8 )
    {
        throw RbException("Error while rejecting RateAgeBetaShift proposal: Node ages were not correctly restored!");
    }
#endif
    
}

/**
 * Reset the move counters. Here we only reset the counter for the number of accepted moves.
 *
 */
void RateAgeBetaShift::resetMoveCounters( void )
{
    num_accepted_current_period = 0;
}


void RateAgeBetaShift::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    if (oldN == tree)
    {
        tree = static_cast<StochasticNode<Tree>* >(newN) ;
    }
    else
    {
        if ( rates == NULL )
        {
            for (size_t i = 0; i < rates_vec.size(); i++)
            {
                if (oldN == rates_vec[i])
                {
                    rates_vec[i] = static_cast<StochasticNode<double>* >(newN);
                }
            }
        }
        else
        {
            if (oldN == rates)
            {
                rates = static_cast<StochasticNode< RbVector<double> >* >(newN);
            }
        }
    }
    
}


void RateAgeBetaShift::setMoveTuningParameter(double tp)
{
    delta = tp;
}


void RateAgeBetaShift::setNumberAcceptedCurrentPeriod( size_t na )
{
    num_accepted_current_period = na;
}


void RateAgeBetaShift::setNumberAcceptedTotal( size_t na )
{
    num_accepted_total = na;
}


void RateAgeBetaShift::tune( void )
{
    
    if ( num_tried_current_period > 2 )
    {
        double rate = num_accepted_current_period / double(num_tried_current_period);
        
        if ( rate > 0.44 )
        {
            delta /= (1.0 + ((rate-0.44)/0.56) );
        }
        else
        {
            delta *= (2.0 - rate/0.44 );
        }
    }

}


