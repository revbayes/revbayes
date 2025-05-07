#include "BirthDeathForwardSimulator.h"

#include <math.h>
#include <iterator>
#include <set>
#include <string>

#include "DistributionExponential.h"
#include "RandomNumberGenerator.h"
#include "RandomNumberFactory.h"
#include "TopologyNode.h"
#include "RbException.h"
#include "StringUtilities.h"
#include "Tree.h"

using namespace RevBayesCore;


BirthDeathForwardSimulator::BirthDeathForwardSimulator( void ) :
    MAX_NUM_LINEAGES(50000),
    complete_tree(false)
{

}


bool BirthDeathForwardSimulator::checkParameters( void ) const
{

    size_t NUM_TIME_INTERVALS = timeline.size();

    if ( lambda.size() != NUM_TIME_INTERVALS )
    {
        throw RbException("You need to provide the same number of time intervals as speciation rates.");
    }

    if ( !(mu.size() == NUM_TIME_INTERVALS || mu.size() == 1) )
    {
        throw RbException("You need to provide the same number of time intervals as extinction rates.");
    }

    if ( !(phi.size() == NUM_TIME_INTERVALS || phi.size() == 1) )
    {
        throw RbException("You need to provide the same number of time intervals as fossilization rates.");
    }

    if ( !(r.size() == NUM_TIME_INTERVALS || r.size() == 1) )
    {
        throw RbException("You need to provide the same number of time intervals as treatment probabilities.");
    }

    if ( !(Lambda.size() == NUM_TIME_INTERVALS || Lambda.size() == 1) )
    {
        throw RbException("You need to provide the same number of time intervals as speciation probabilities.");
    }

    if ( !(Mu.size() == NUM_TIME_INTERVALS || Mu.size() == 1) )
    {
        throw RbException("You need to provide the same number of time intervals as extinction probabilities.");
    }

    if ( !(Phi.size() == NUM_TIME_INTERVALS || Phi.size() == 1) )
    {
        throw RbException("You need to provide the same number of time intervals as fossilization probabilities.");
    }

    if ( !(R.size() == NUM_TIME_INTERVALS || R.size() == 1) )
    {
        throw RbException("You need to provide the same number of time intervals as event-treatment probabilities.");
    }

    return true;
}

size_t BirthDeathForwardSimulator::getNumberOfCategories( void ) const
{

    if ( lambda.size() > 0 )
    {
        return lambda[0].size();
    }
    else
    {
        return 0;
    }

}


double BirthDeathForwardSimulator::getLambdaProbability( size_t index, size_t n ) const
{
    // There is no burst speciation at the present day
    if ( index == 0 )
    {
        return 0.0;
    }

    if ( Lambda.size() > index )
    {
        if ( Lambda[index].size() > n )
        {
            return Lambda[index][n];
        }
        else if ( Lambda[index].size() == 1 )
        {
            return Lambda[index][0];
        }
        else
        {
            throw RbException("Problem in the birth-death simulation: We could not get the speciation event probabilities.");
        }
    }
    else if ( Lambda.size() == 1 )
    {
        if ( Lambda[0].size() > n )
        {
            return Lambda[0][n];
        }
        else if ( Lambda[0].size() == 1 )
        {
            return Lambda[0][0];
        }
        else
        {
            throw RbException("Problem in the birth-death simulation: We could not get the speciation event probabilities.");
        }
    }
    else
    {
        throw RbException("Problem in the birth-death simulation: We could not get the speciation event probabilities.");
    }
}


std::vector<double> BirthDeathForwardSimulator::getLambdaRate( size_t index, size_t n ) const
{
    if ( lambda.size() > index )
    {
        if ( lambda[index].size() > 1 || lambda[index].size() == n )
        {
            return lambda[index];
        }
        else if ( lambda[index].size() == 1 )
        {
            return std::vector<double>(n, lambda[index][0] );
        }
        else
        {
            throw RbException("Problem in the birth-death simulation: We could not get the speciation rates.");
        }
    }
    else if ( lambda.size() == 1 )
    {
        if ( lambda[0].size() > 1 || lambda[index].size() == n )
        {
            return lambda[0];
        }
        else if ( lambda[0].size() == 1 )
        {
            return std::vector<double>(n, lambda[0][0] );
        }
        else
        {
            throw RbException("Problem in the birth-death simulation: We could not get the speciation rates.");
        }
    }
    else
    {
        throw RbException("Problem in the birth-death simulation: We could not get the speciation rates.");
    }

}


double BirthDeathForwardSimulator::getMuProbability( size_t index, size_t n ) const
{
    
    // std::cout << "getting mu probability for index " << index << " with n = " << n << std::endl;

    // There is no mass extinction at the present day
    if ( index == 0 )
    {
        return 0.0;
    }

    if ( Mu.size() > index )
    {
        if ( Mu[index].size() > n )
        {
            return Mu[index][n];
        }
        else if ( Mu[index].size() == 1 )
        {
            return Mu[index][0];
        }
        else
        {
            throw RbException("Problem in the birth-death simulation: We could not get the extinction event probabilities.");
        }
    }
    else if ( Mu.size() == 1 )
    {
        if ( Mu[0].size() > n )
        {
            return Mu[0][n];
        }
        else if ( Mu[0].size() == 1 )
        {
            return Mu[0][0];
        }
        else
        {
            throw RbException("Problem in the birth-death simulation: We could not get the extinction event probabilities.");
        }
    }
    else
    {
        throw RbException("Problem in the birth-death simulation: We could not get the extinction event probabilities.");
    }

}


std::vector<double> BirthDeathForwardSimulator::getMuRate( size_t index, size_t n ) const
{

    if ( mu.size() > index )
    {
        if ( mu[index].size() > 1 || mu[index].size() == n )
        {
            return mu[index];
        }
        else if ( mu[index].size() == 1 )
        {
            return std::vector<double>(n, mu[index][0] );
        }
        else
        {
            throw RbException("Problem in the birth-death simulation: We could not get the extinction rates.");
        }
    }
    else if ( mu.size() == 1 )
    {
        if ( mu[0].size() > 1 || mu[index].size() == n )
        {
            return mu[0];
        }
        else if ( mu[0].size() == 1 )
        {
            return std::vector<double>(n, mu[0][0] );
        }
        else
        {
            throw RbException("Problem in the birth-death simulation: We could not get the extinction rates.");
        }
    }
    else
    {
        throw RbException("Problem in the birth-death simulation: We could not get the extinction rates.");
    }

}


double BirthDeathForwardSimulator::getPhiProbability( size_t time_index, size_t cat_index ) const
{
    if ( Phi.size() > time_index )
    {
        if ( Phi[time_index].size() > cat_index )
        {
            return Phi[time_index][cat_index];
        }
        else if ( Phi[time_index].size() == 1 )
        {
            return Phi[time_index][0];
        }
        else
        {
            throw RbException("Problem in the birth-death simulation: We could not get the sampling event probabilities.");
        }
    }
    else if ( time_index > 0 )
    {
        return 0.0;     // we return by default 0.0
    }
    else if ( Phi.size() == 1 )
    {
        if ( Phi[0].size() > cat_index )
        {
            return Phi[0][cat_index];
        }
        else if ( Phi[0].size() == 1 )
        {
            return Phi[0][0];
        }
        else
        {
            throw RbException("Problem in the birth-death simulation: We could not get the sampling event probabilities.");
        }
    }
    else
    {
        throw RbException("Problem in the birth-death simulation: We could not get the sampling event probabilities.");
    }

}


std::vector<double> BirthDeathForwardSimulator::getPhiRate( size_t index, size_t n ) const
{

    if ( phi.size() > index )
    {
        if ( phi[index].size() > 1 || phi[index].size() == n )
        {
            return phi[index];
        }
        else if ( phi[index].size() == 1 )
        {
            return std::vector<double>(n, phi[index][0] );
        }
        else
        {
            throw RbException("Problem in the birth-death simulation: We could not get the sampling rate.");
        }
    }
    else if ( phi.size() == 1 )
    {
        if ( phi[0].size() > 1 || phi[index].size() == n )
        {
            return phi[0];
        }
        else if ( phi[0].size() == 1 )
        {
            return std::vector<double>(n, phi[0][0] );
        }
        else
        {
            throw RbException("Problem in the birth-death simulation: We could not get the sampling rate.");
        }
    }
    else
    {
        throw RbException("Problem in the birth-death simulation: We could not get the sampling rate.");
    }

}


double BirthDeathForwardSimulator::getRProbability( size_t index, size_t n ) const
{

    // There is no death-after-sampling at the present day
    if ( index == 0 )
    {
        return 0.0;
    }

    if ( R.size() > index )
    {
        if ( R[index].size() > n )
        {
            return R[index][n];
        }
        else if ( R[index].size() == 1 )
        {
            return R[index][0];
        }
        else
        {
            throw RbException("Problem in the birth-death simulation: We could not get the treatment event probabilities.");
        }
    }
    else if ( R.size() == 1 )
    {
        if ( R[0].size() > n )
        {
            return R[0][n];
        }
        else if ( R[0].size() == 1 )
        {
            return R[0][0];
        }
        else
        {
            throw RbException("Problem in the birth-death simulation: We could not get the treatment event probabilities.");
        }
    }
    else
    {
        throw RbException("Problem in the birth-death simulation: We could not get the treament event probabilities.");
    }

}


std::vector<double> BirthDeathForwardSimulator::getRRate( size_t index, size_t n ) const
{

    if ( r.size() > index )
    {
        if ( r[index].size() > 1 || r[index].size() == n )
        {
            return r[index];
        }
        else if ( r[index].size() == 1 )
        {
            return std::vector<double>(n, r[index][0] );
        }
        else
        {
            throw RbException("Problem in the birth-death simulation: We could not get the treatment rate.");
        }
    }
    else if ( r.size() == 1 )
    {
        if ( r[0].size() > 1 || r[index].size() == n )
        {
            return r[0];
        }
        else if ( r[0].size() == 1 )
        {
            return std::vector<double>(n, r[0][0] );
        }
        else
        {
            throw RbException("Problem in the birth-death simulation: We could not get the treatment rate.");
        }
    }
    else
    {
        throw RbException("Problem in the birth-death simulation: We could not get the treatment rate.");
    }

}


std::vector<double> BirthDeathForwardSimulator::getRootCategoryProbabilities( size_t n ) const
{
    if ( root_cat_probability.size() == n )
    {
        return root_cat_probability;
    }
    else if ( root_cat_probability.size() == 1 )
    {
        return std::vector<double>(n, root_cat_probability[0] );
    }
    else
    {
        throw RbException("Problem in the birth-death simulation: We could not get the root category probabilities.");
    }

}


bool BirthDeathForwardSimulator::hasExtantSurvivor(const TopologyNode &n) const
{
    if ( n.isTip() == true )
    {
        return n.getAge() == 0.0;
    }
    else
    {
        bool survivor = false;
        for (size_t i=0; i<n.getNumberOfChildren(); ++i)
        {
            survivor |= hasExtantSurvivor( n.getChild(i) );
        }
        return survivor;
    }
}


Tree* BirthDeathForwardSimulator::simulateTreeConditionTime(double start_age, SIM_CONDITION condition) const
{

    if ( checkParameters() == false )
    {

    }
    
    size_t NUM_TIME_INTERVALS = timeline.size();
    size_t NUM_CATEGORIES = getNumberOfCategories();

    RandomNumberGenerator *rng = GLOBAL_RNG;

    TopologyNode *root = NULL;
    do {

        // we define here a set of all sampled nodes
        // this set will be necessary for sampled fossils
        std::set<TopologyNode*> sampled_nodes;
        std::set<TopologyNode*> extinct_not_sampled_nodes;

        std::vector< std::set<TopologyNode*> > active_nodes_in_actegories = std::vector< std::set<TopologyNode*> >( NUM_CATEGORIES, std::set<TopologyNode*>() );

        int current_time_index = int(NUM_TIME_INTERVALS) - 1;
        while ( current_time_index > 0 && timeline[current_time_index] > start_age )
        {
            --current_time_index;
        }

        // draw the index of the root category
        size_t root_cat_index = 0;
        std::vector<double> root_cat_probs = getRootCategoryProbabilities( NUM_CATEGORIES );
        double u = rng->uniform01();
        while ( root_cat_index < NUM_CATEGORIES && u > root_cat_probs[root_cat_index] )
        {
            u -= root_cat_probs[root_cat_index];
            ++root_cat_index;
        }

        // for now we assume that we always start with the root
        // we should make the an argument
        size_t current_num_active_nodes = 2;
        root = new TopologyNode();
        root->setAge( start_age );
        // left child
        TopologyNode *left = new TopologyNode();
        root->addChild( left );
        left->setParent( root );
        // right child
        TopologyNode *right = new TopologyNode();
        root->addChild( right );
        right->setParent( root );
        
        root->setSerialSpeciation( true );

        // now add the two children to the active nodes
        active_nodes_in_actegories[root_cat_index].insert( left  );
        active_nodes_in_actegories[root_cat_index].insert( right );

        std::vector<double> current_lambda   = getLambdaRate(current_time_index, NUM_CATEGORIES);
        std::vector<double> current_mu       = getMuRate(current_time_index, NUM_CATEGORIES);
        std::vector<double> current_phi      = getPhiRate(current_time_index, NUM_CATEGORIES);
        std::vector<double> current_r        = getRRate(current_time_index, NUM_CATEGORIES);

        double current_lambda_total     = 0;
        double current_mu_total         = 0;
        double current_phi_total        = 0;
        for ( size_t i=0; i<NUM_CATEGORIES; ++i)
        {
            double num_lineages_in_category = double( active_nodes_in_actegories[i].size() );
            current_lambda_total    += num_lineages_in_category * current_lambda[i];
            current_mu_total        += num_lineages_in_category * current_mu[i];
            current_phi_total       += num_lineages_in_category * current_phi[i];
        }
        double current_rate_total = current_lambda_total + current_mu_total + current_phi_total;
        

        double current_age = start_age;
        while ( current_age > 0.0 && current_num_active_nodes > 0 )
        {
            // draw the next event
            double next_age = current_age - RbStatistics::Exponential::rv(current_rate_total, *rng);

            // check if the next event is passed the current epoch
            if ( next_age < timeline[current_time_index] )
            {
                double event_age = timeline[current_time_index];
                // now we perform the different events

                for ( size_t i=0; i<NUM_CATEGORIES; ++i )
                {
                    // (1) a burst event
                    double this_burst_prob = getLambdaProbability(current_time_index, i);
                    if ( this_burst_prob > 0.0 )
                    {
                        std::set<TopologyNode*> old_active_nodes_in_category = active_nodes_in_actegories[i];
                        for ( std::set<TopologyNode*>::const_iterator it=old_active_nodes_in_category.begin(); it!=old_active_nodes_in_category.end(); ++it)
                        {
                            u = rng->uniform01();
                            if ( u < this_burst_prob )
                            {
                                TopologyNode *this_node = *it;
                                this_node->setAge( event_age );
                                // left child
                                TopologyNode *left = new TopologyNode();
                                this_node->addChild( left );
                                left->setParent( this_node );
                                // right child
                                TopologyNode *right = new TopologyNode();
                                this_node->addChild( right );
                                right->setParent( this_node );

                                // update the active node set
                                active_nodes_in_actegories[i].erase( this_node );
                                active_nodes_in_actegories[i].insert( left );
                                active_nodes_in_actegories[i].insert( right );
                                
                                this_node->setBurstSpeciation( true );
                                

                            } // end-if there was a burst speciation event for this node

                        } // end-for all nodes in this rate category

                    } // end-if this category has a non-zero burst probability


                    // (2) a mass-extinction event
                    double this_mass_extinction_prob = getMuProbability(current_time_index, i);
                    if ( this_mass_extinction_prob > 0.0 )
                    {
                        std::set<TopologyNode*> old_active_nodes_in_category = active_nodes_in_actegories[i];
                        for ( std::set<TopologyNode*>::const_iterator it=old_active_nodes_in_category.begin(); it!=old_active_nodes_in_category.end(); ++it)
                        {
                            u = rng->uniform01();
                            if ( u < this_mass_extinction_prob )
                            {
                                TopologyNode *this_node = *it;
                                this_node->setAge( event_age );

                                // update the active node set
                                active_nodes_in_actegories[i].erase( this_node );

                                // remember this node to prune later
                                extinct_not_sampled_nodes.insert( this_node );

                            } // end-if there was a mass-extinction event for this node

                        } // end-for all nodes in this rate category

                    } // end-if this category has a non-zero mass-extinction probability


                    // (3) a sampling event
                    double this_sampling_prob = getPhiProbability(current_time_index, i);
                    if ( this_sampling_prob > 0.0 )
                    {
                        double this_sampling_extinction_prob = getRProbability(current_time_index, i);
                        std::set<TopologyNode*> old_active_nodes_in_category = active_nodes_in_actegories[i];
                        for ( std::set<TopologyNode*>::const_iterator it=old_active_nodes_in_category.begin(); it!=old_active_nodes_in_category.end(); ++it)
                        {
                            u = rng->uniform01();
                            if ( u < this_sampling_prob )
                            {
                                TopologyNode *this_node = *it;
                                this_node->setAge( event_age );
                                
//                                this_node->setSamplingEvent( true );
                                
                                // store this node
                                sampled_nodes.insert( this_node );
                                
                                if ( event_age > 0.0 )
                                {

                                    // simulate also if this sampling lead to a death event (i.e., treatment)
                                    u = rng->uniform01();
                                    if ( u < this_sampling_extinction_prob )
                                    {
                                        // update the active node set
                                        active_nodes_in_actegories[i].erase( this_node );
                                    }
                                    else
                                    {
                                        // left child (sampled ancestor)
                                        TopologyNode *left = new TopologyNode();
                                        this_node->addChild( left );
                                        left->setParent( this_node );
                                    left->setAge( event_age );
                                    left->setSampledAncestor( true );
                                    left->setSamplingEvent( true );
                                    sampled_nodes.insert( left );

                                    // right child
                                    TopologyNode *right = new TopologyNode();
                                    this_node->addChild( right );
                                    right->setParent( this_node );

                                    // update the active node set
                                    active_nodes_in_actegories[i].erase( this_node );
                                    active_nodes_in_actegories[i].insert( right );
                                    
                                }
                                }

                            } // end-if there was a sampling event for this node

                        } // end-for all nodes in this rate category

                    } // end-if this category has a non-zero sampling probability

                } // end-for over all categories

                // update the age
                current_age = event_age;

                // update the interval index
                --current_time_index;

//                if ( current_time_index < NUM_TIME_INTERVALS )
                if ( current_time_index >= 0 )
                {
                    // update all the rate summaries
                    current_lambda   = getLambdaRate(current_time_index, NUM_CATEGORIES);
                    current_mu       = getMuRate(current_time_index, NUM_CATEGORIES);
                    current_phi      = getPhiRate(current_time_index, NUM_CATEGORIES);
                    current_r        = getRRate(current_time_index, NUM_CATEGORIES);

                    current_num_active_nodes = 0;

                    current_lambda_total     = 0;
                    current_mu_total         = 0;
                    current_phi_total        = 0;
                    for ( size_t i=0; i<NUM_CATEGORIES; ++i)
                    {
                        current_num_active_nodes += active_nodes_in_actegories[i].size();

                        double num_lineages_in_category = double( active_nodes_in_actegories[i].size() );
                        current_lambda_total    += num_lineages_in_category * current_lambda[i];
                        current_mu_total        += num_lineages_in_category * current_mu[i];
                        current_phi_total       += num_lineages_in_category * current_phi[i];
                    }
                    current_rate_total = current_lambda_total + current_mu_total + current_phi_total;
                }

            } // end-if the time passed the epoch boundary
            else
            {

                // check which type of event we had
                u = rng->uniform01();
                if ( u < (current_lambda_total)/current_rate_total )
                {
                    // it was a speciation event

                    u = rng->uniform01();
                    for ( size_t i=0; i<NUM_CATEGORIES; ++i)
                    {
                        double num_lineages_in_category = double( active_nodes_in_actegories[i].size() );
                        double this_cat_lambda = current_lambda[i];
                        if ( u < (num_lineages_in_category*this_cat_lambda)/current_lambda_total )
                        {
                            u = rng->uniform01();
                            size_t this_node_index = size_t( floor(u*num_lineages_in_category) );
                            std::set<TopologyNode*>::const_iterator it = active_nodes_in_actegories[i].begin();
                            std::advance(it, this_node_index);

                            TopologyNode *this_node = *it;
                            this_node->setAge( next_age );
                            // left child
                            TopologyNode *left = new TopologyNode();
                            this_node->addChild( left );
                            left->setParent( this_node );
                            // right child
                            TopologyNode *right = new TopologyNode();
                            this_node->addChild( right );
                            right->setParent( this_node );

                            // update the active node set
                            active_nodes_in_actegories[i].erase( this_node );
                            active_nodes_in_actegories[i].insert( left );
                            active_nodes_in_actegories[i].insert( right );

                            // update counters
                            ++current_num_active_nodes;
                            current_lambda_total    += this_cat_lambda;
                            current_mu_total        += current_mu[i];
                            current_phi_total       += current_phi[i];
                            current_rate_total      += this_cat_lambda + current_mu[i] + current_phi[i];

                            this_node->setSerialSpeciation( true );
                            

                            // now stop the loop
                            break;
                        }
                        else
                        {
                            u -= (num_lineages_in_category*this_cat_lambda)/current_lambda_total;
                        } // end-if-else the speciation event was in this category

                    } // end-for loop over all categories

                } // end-if it was a speciation event
                else if ( u < (current_lambda_total+current_mu_total)/current_rate_total )
                {
                    // it was an extinction event

                    u = rng->uniform01();
                    for ( size_t i=0; i<NUM_CATEGORIES; ++i)
                    {
                        double num_lineages_in_category = double( active_nodes_in_actegories[i].size() );
                        double this_cat_mu = current_mu[i];
                        if ( u < (num_lineages_in_category*this_cat_mu)/current_mu_total )
                        {
                            size_t this_node_index = size_t( floor(u*num_lineages_in_category) );
                            std::set<TopologyNode*>::const_iterator it = active_nodes_in_actegories[i].begin();
                            std::advance(it, this_node_index);

                            TopologyNode *this_node = *it;
                            this_node->setAge( next_age );

                            // update the active node set
                            active_nodes_in_actegories[i].erase( this_node );

                            // remember this node to prune later
                            extinct_not_sampled_nodes.insert( this_node );

                            // update counters
                            --current_num_active_nodes;
                            current_lambda_total    -= current_lambda[i];
                            current_mu_total        -= this_cat_mu;
                            current_phi_total       -= current_phi[i];
                            current_rate_total      -= (current_lambda[i] + this_cat_mu + current_phi[i]);

                            // now stop the loop
                            break;
                        }
                        else
                        {
                            u -= (num_lineages_in_category*this_cat_mu)/current_mu_total;
                        } // end-if-else the speciation event was in this category

                    } // end-for loop over all categories

                } // end-if it was an extinction event
                else if ( u < (current_lambda_total+current_mu_total+current_phi_total)/current_rate_total )
                {
                    // it was a sampling event
                    u = rng->uniform01();
                    for ( size_t i=0; i<NUM_CATEGORIES; ++i)
                    {
                        double num_lineages_in_category = double( active_nodes_in_actegories[i].size() );
                        double this_cat_phi = current_phi[i];
                        double this_cat_r = current_r[i];
                        if ( u < (num_lineages_in_category*this_cat_phi)/current_phi_total )
                        {
                            size_t this_node_index = size_t( floor(u*num_lineages_in_category) );
                            std::set<TopologyNode*>::const_iterator it = active_nodes_in_actegories[i].begin();
                            std::advance(it, this_node_index);

                            TopologyNode *this_node = *it;
                            this_node->setAge( next_age );
                            
                            // // the current node to our set of sampled nodes
                            // sampled_nodes.insert( this_node );

                            // now draw if this node died during sampling (i.e., treatment)
                            u = rng->uniform01();
                            if ( u < this_cat_r )
                            {
                                // update the active node set
                                active_nodes_in_actegories[i].erase( this_node );

                                // the current node to our set of sampled nodes
                                sampled_nodes.insert( this_node );

                                // update counters
                                --current_num_active_nodes;
                                current_lambda_total    -= current_lambda[i];
                                current_mu_total        -= current_mu[i];
                                current_phi_total       -= this_cat_phi;
                                current_rate_total      -= (current_lambda[i] + current_mu[i] + this_cat_phi);
                            } // end-of this node died during sampling
                            else
                            {
                                // left child (sampled ancestor)
                                TopologyNode *left = new TopologyNode();
                                this_node->addChild( left );
                                left->setParent( this_node );
                                left->setAge( next_age );
                                left->setSampledAncestor( true );
                                left->setSerialSampling( true );
                                sampled_nodes.insert( left );

                                // right child
                                TopologyNode *right = new TopologyNode();
                                this_node->addChild( right );
                                right->setParent( this_node );

                                // update the active node set
                                active_nodes_in_actegories[i].erase( this_node );
                                active_nodes_in_actegories[i].insert( right );
                                
                            }
                            // now stop the loop
                            break;
                        }
                        else
                        {
                            u -= (num_lineages_in_category*this_cat_phi)/current_phi_total;
                        } // end-if-else the speciation event was in this category

                    } // end-for loop over all categories

                } // end-if it was a sampling event
                else
                {
                    // this should not happen
                    throw RbException("Problem in birth-death forward simulator: could not determine which type of event happened");
                }

                current_age = next_age;
            }

            // Ensure tree doesn't get too big to manage
            if ( current_num_active_nodes > MAX_NUM_LINEAGES )
            {
                delete root;
                root = NULL;
                break;
            }

        } // end-while we still have time and nodes to simulate forward


        // If we have not already given up because the tree got too big, we can now clean up the tree and check if it meets conditioning standards
        if ( root != NULL ) 
        {
            // set the ages of all remaining active nodes to 0.0
            // add nodes that were not sampled at the present to the extinct pile
            for ( size_t i=0; i<NUM_CATEGORIES; ++i)
            {
                const std::set<TopologyNode*> &active_nodes_in_category = active_nodes_in_actegories[i];
                for ( std::set<TopologyNode*>::const_iterator it=active_nodes_in_category.begin(); it!=active_nodes_in_category.end(); ++it)
                {
                    TopologyNode *this_node = *it;
                    this_node->setAge( 0.0 );

                    bool found = sampled_nodes.find( this_node ) != sampled_nodes.end();
                    if ( found == false )
                    {
                        extinct_not_sampled_nodes.insert(this_node);
                    }
                }
            }
            
            // Complete trees should still meet the user's condition choice,
            // but checking requires destructive pruning
            TopologyNode *complete_root = NULL;
            if ( complete_tree == true )
            {
                complete_root = root->clone();
            }

            bool prune = true;
            // now prune away all the extinct nodes
            if ( prune == true )
            {
                for ( std::set<TopologyNode*>::const_iterator it=extinct_not_sampled_nodes.begin(); it!=extinct_not_sampled_nodes.end(); ++it)
                {
                    TopologyNode *this_node = *it;

                    if ( this_node->isRoot() == false ) // this should truly never happen, only if there was no survivor
                    {
                        TopologyNode &parent = this_node->getParent();
                        
                        if ( parent.getNumberOfChildren() == 1 )
                        {

                            parent.removeChild( this_node );
                            delete this_node;
                            
                        }
                        else
                        {
                            TopologyNode *sibling = &parent.getChild( 0 );
                            if ( sibling == this_node ) sibling = &parent.getChild( 1 );

                            if ( parent.isRoot() == false )
                            {
                                TopologyNode &grandparent = parent.getParent();
                                TopologyNode *uncle = &grandparent.getChild( 0 );
                                if ( uncle == &parent ) uncle = &grandparent.getChild( 1 );

                                parent.removeChild( sibling );
                                grandparent.removeChild( &parent );
                                parent.setParent( NULL );

                                sibling->setParent( &grandparent );
                                grandparent.addChild( sibling );

                                sibling->setSampledAncestor(false);
                                
                                delete &parent; // this will also delete this child node

                            }
                            else
                            {
                                parent.removeChild( sibling );
                                sibling->setParent( NULL );
                                root = sibling;
                                
                                delete &parent;

                            }

                        }

                    }

                }

            }

            if ( condition == SURVIVAL && current_num_active_nodes <= 1 )
            {

                delete root;
                delete complete_root;
                root = NULL;
                
            }
            // next we check that both sides of the root node were sampled ()
            else if ( condition == ROOT && (root->getAge() < start_age ||
                                            current_num_active_nodes <= 1 ||
                                            hasExtantSurvivor(root->getChild(0)) == false ||
                                            hasExtantSurvivor(root->getChild(1)) == false ) )
            {

                delete root;
                delete complete_root;
                root = NULL;
                
            }

            // The sampled portion of the complete tree will meet the user's condition choice
            if ( root != NULL && complete_tree == true )
            {

                // we need to free the old memory
                delete root;
                root = complete_root;

            }

        }

    } while ( root == NULL );

    // build and return the tree
    Tree *my_tree = new Tree();
    my_tree->setRoot(root, true);
    my_tree->setRooted( true );
    
    size_t num_nodes = my_tree->getNumberOfNodes();

    for (size_t i=0; i<num_nodes; ++i)
    {
        TopologyNode &n = my_tree->getNode(i);
        if ( n.isTip() == true )
        {
            n.setName("Tip_"+ StringUtilities::toString(i+1) );
        }
        else if ( n.isSampledAncestorTip() == true )
        {
            n.setName("Tip_"+ StringUtilities::toString(i+1) );
        }
        
    }

    return my_tree;
}


void BirthDeathForwardSimulator::setBurstProbability( const std::vector<std::vector< double > > &l )
{
    Lambda = l;
}

void BirthDeathForwardSimulator::setCompleteTree( bool c )
{
    complete_tree = c;
}

void BirthDeathForwardSimulator::setExtinctionRate( const std::vector<std::vector< double > > &m )
{
    mu = m;
}

void BirthDeathForwardSimulator::setMaxNumLineages( size_t m )
{
    MAX_NUM_LINEAGES = m;
}


void BirthDeathForwardSimulator::setMassExtinctionProbability( const std::vector<std::vector< double > > &m )
{
    Mu = m;
}


void BirthDeathForwardSimulator::setRootCategoryProbabilities( const std::vector<double> &p )
{
    root_cat_probability = p;
}



void BirthDeathForwardSimulator::setSamplingProbability( const std::vector<std::vector< double > > &p )
{
    Phi = p;
}


void BirthDeathForwardSimulator::setSamplingExtinctionProbability( const std::vector<std::vector< double > > &r )
{
    R = r;
}


void BirthDeathForwardSimulator::setSamplingRate( const std::vector<std::vector< double > > &p )
{
    phi = p;
}


void BirthDeathForwardSimulator::setSamplingExtinctionRate( const std::vector<std::vector< double > > &r )
{
    this->r = r;
}


void BirthDeathForwardSimulator::setSpeciationRate( const std::vector<std::vector< double > > &l )
{
    lambda = l;
}


void BirthDeathForwardSimulator::setTimeline( const std::vector< double > &t )
{
    timeline = t;
}
