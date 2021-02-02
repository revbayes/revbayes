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

RateAgeBetaShift::RateAgeBetaShift(StochasticNode<Tree> *tr, std::vector<StochasticNode<double> *> v, double d, bool t, double w) : AbstractMove( w, t),
    tree( tr ),
    rates( v ),
    delta( d ),
    storedNode( NULL ),
    storedAge( 0.0 ),
    storedRates( rates.size(), 0.0 ),
    num_accepted_current_period( 0 ),
    num_accepted_total( 0 )
{
    
    addNode( tree );
    for (std::vector<StochasticNode<double>* >::iterator it = rates.begin(); it != rates.end(); ++it)
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
    bool checkLikelihoodShortcuts = rng->uniform01() < 0.001;
    if ( checkLikelihoodShortcuts == true )
    {
        for (RbOrderedSet<DagNode*>::iterator it = affected.begin(); it != affected.end(); ++it)
        {
            (*it)->touch();
            oldLnLike += (*it)->getLnProbability();
        }
    }
    
    // pick a random node which is not the root and neithor the direct descendant of the root
    TopologyNode* node;
    size_t nodeIdx = 0;
    do {
        double u = rng->uniform01();
        nodeIdx = size_t( std::floor(tau.getNumberOfNodes() * u) );
        node = &tau.getNode(nodeIdx);
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
    storedNode = node;
    storedAge = my_age;
    
    
    storedRates[nodeIdx] = rates[nodeIdx]->getValue();
    for (size_t i = 0; i < node->getNumberOfChildren(); i++)
    {
        size_t childIdx = node->getChild(i).getIndex();
        storedRates[childIdx] = rates[childIdx]->getValue();
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
    tau.getNode(nodeIdx).setAge( my_new_age );
    
    // touch the tree so that the likelihoods are getting stored
    tree->touch();
    
    // get the probability ratio of the tree
    double treeProbRatio = tree->getLnProbabilityRatio();
    
    
    // set the rates
    double pa = node->getParent().getAge();
    double my_new_rate = (pa - my_age) * storedRates[nodeIdx] / (pa - my_new_age);
    
    // now we set the new value
    // this will automatically call a touch
    rates[nodeIdx]->setValue( new double( my_new_rate ) );
    
    // get the probability ratio of the new rate
    double ratesProbRatio = rates[nodeIdx]->getLnProbabilityRatio();
    double jacobian = log((pa - my_age) / (pa - my_new_age));
    
    for (size_t i = 0; i < node->getNumberOfChildren(); i++)
    {
        size_t childIdx = node->getChild(i).getIndex();
        double a = node->getChild(i).getAge();
        double child_new_rate = (my_age - a) * storedRates[childIdx] / (my_new_age - a);
        
        // now we set the new value
        // this will automatically call a touch
        rates[childIdx]->setValue( new double( child_new_rate ) );

        // get the probability ratio of the new rate
        ratesProbRatio += rates[childIdx]->getLnProbabilityRatio();
        jacobian += log((my_age - a) / (my_new_age - a));
        
    }
    
    if ( checkLikelihoodShortcuts == true )
    {
        double lnProbRatio = 0;
        double newLnLike = 0;
        for (RbOrderedSet<DagNode*>::iterator it = affected.begin(); it != affected.end(); ++it)
        {

            double tmp = (*it)->getLnProbabilityRatio();
            lnProbRatio += tmp;
            newLnLike += (*it)->getLnProbability();
        }
    
        if ( fabs(lnProbRatio) > 1E-8 )
        {
            double lnProbRatio2 = 0;
            double newLnLike2 = 0;
            for (RbOrderedSet<DagNode*>::iterator it = affected.begin(); it != affected.end(); ++it)
            {
                
                double tmp2 = (*it)->getLnProbabilityRatio();
                lnProbRatio2 += tmp2;
                newLnLike2 += (*it)->getLnProbability();
            }
            
            throw RbException("Likelihood shortcut computation failed in rate-age-proposal.");
        }
        
    }
    
    double hastingsRatio = backward - forward + jacobian;
    double ln_acceptance_ratio = pHeat * prHeat * (treeProbRatio + ratesProbRatio) + hastingsRatio;
    
    if (ln_acceptance_ratio >= 0.0)
    {
        num_accepted_total++;
        num_accepted_current_period++;
        
        tree->keep();
        rates[nodeIdx]->keep();
        for (size_t i = 0; i < node->getNumberOfChildren(); i++)
        {
            size_t childIdx = node->getChild(i).getIndex();
            rates[childIdx]->keep();
        }
    }
    else if (ln_acceptance_ratio < -300.0)
    {
        reject();
        tree->restore();
        rates[nodeIdx]->restore();
        for (size_t i = 0; i < node->getNumberOfChildren(); i++)
        {
            size_t childIdx = node->getChild(i).getIndex();
            rates[childIdx]->restore();
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
            tree->keep();
            rates[nodeIdx]->keep();
            for (size_t i = 0; i < node->getNumberOfChildren(); i++)
            {
                size_t childIdx = node->getChild(i).getIndex();
                rates[childIdx]->keep();
            }
        }
        else
        {
            reject();
            tree->restore();
            rates[nodeIdx]->restore();
            for (size_t i = 0; i < node->getNumberOfChildren(); i++)
            {
                size_t childIdx = node->getChild(i).getIndex();
                rates[childIdx]->restore();
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
    tree->getValue().getNode( storedNode->getIndex() ).setAge( storedAge );
    
    // undo the rates
    size_t nodeIdx = storedNode->getIndex();
    rates[nodeIdx]->setValue(new double(storedRates[nodeIdx]));
    for (size_t i = 0; i < storedNode->getNumberOfChildren(); i++)
    {
        size_t childIdx = storedNode->getChild(i).getIndex();
        rates[childIdx]->setValue(new double(storedRates[childIdx]));
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
        for (size_t i = 0; i < rates.size(); i++)
        {
            if (oldN == rates[i])
            {
                rates[i] = static_cast<StochasticNode<double>* >(newN);
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


