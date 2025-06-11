#include <cstdlib>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "MultispeciesCoalescentMigration.h"
#include "ConstantNode.h"
#include "DistributionExponential.h"
#include "MultispeciesCoalescentMigrationODE.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbMathCombinatorialFunctions.h"
#include "RbMathLogic.h"
#include "TopologyNode.h"
#include "RateGenerator.h"
#include "RbException.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "Taxon.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

#include "boost/numeric/odeint.hpp"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

MultispeciesCoalescentMigration::MultispeciesCoalescentMigration(const TypedDagNode<Tree> *sp, const std::vector<Taxon> &t, const TypedDagNode<RateGenerator>* q, const TypedDagNode<double >* d) : TypedDistribution<Tree>( NULL ),
    taxa(t),
    species_tree( sp ),
    Ne( new ConstantNode<double>("", new double(1.0) ) ),
    Nes( NULL ),
    Q(q),
    delta( d ),
    num_taxa( taxa.size() ),
    log_tree_topology_prob (0.0)
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( species_tree );
    addParameter( Q );
    addParameter( delta );

    std::set<std::string> species_names;
    for (std::vector<Taxon>::const_iterator it=taxa.begin(); it!=taxa.end(); ++it)
    {
        species_names.insert( it->getSpeciesName() );
    }
    
    double ln_fact = RbMath::lnFactorial((int)(num_taxa));
    
    log_tree_topology_prob = (num_taxa - 1) * RbConstants::LN2 - 2.0 * ln_fact - std::log( num_taxa ) ;
    
    redrawValue();
    
}


MultispeciesCoalescentMigration::~MultispeciesCoalescentMigration()
{
    
}



void MultispeciesCoalescentMigration::attachTimes(Tree *psi, std::vector<TopologyNode *> &tips, size_t index, const std::vector<double> &times) {
    
    if (index < num_taxa-1)
    {
        // Get the rng
        RandomNumberGenerator* rng = GLOBAL_RNG;
        
        // randomly draw one node from the list of tips
        size_t tip_index = static_cast<size_t>( floor(rng->uniform01()*tips.size()) );
        
        // get the node from the list
        TopologyNode* parent = tips.at(tip_index);
        psi->getNode( parent->getIndex() ).setAge( times[num_taxa - index - 2] );
        
        // remove the randomly drawn node from the list
        tips.erase(tips.begin()+tip_index);
        
        // add a left child
        TopologyNode* leftChild = &parent->getChild(0);
        if ( !leftChild->isTip() )
        {
            tips.push_back(leftChild);
        }
        
        // add a right child
        TopologyNode* rightChild = &parent->getChild(1);
        if ( !rightChild->isTip() )
        {
            tips.push_back(rightChild);
        }
        
        // recursive call to this function
        attachTimes(psi, tips, index+1, times);
    }
}


void MultispeciesCoalescentMigration::buildRandomBinaryTree(std::vector<TopologyNode*> &tips)
{
    
    if (tips.size() < num_taxa)
    {
        // Get the rng
        RandomNumberGenerator* rng = GLOBAL_RNG;
        
        // randomly draw one node from the list of tips
        size_t index = static_cast<size_t>( floor(rng->uniform01()*tips.size()) );
        
        // get the node from the list
        TopologyNode* parent = tips.at(index);
        
        // remove the randomly drawn node from the list
        tips.erase(tips.begin()+index);
        
        // add a left child
        TopologyNode* leftChild = new TopologyNode();
        parent->addChild(leftChild);
        leftChild->setParent(parent);
        tips.push_back(leftChild);
        
        // add a right child
        TopologyNode* rightChild = new TopologyNode();
        parent->addChild(rightChild);
        rightChild->setParent(parent);
        tips.push_back(rightChild);
        
        // recursive call to this function
        buildRandomBinaryTree(tips);
    }
}


MultispeciesCoalescentMigration* MultispeciesCoalescentMigration::clone( void ) const
{
    
    return new MultispeciesCoalescentMigration( *this );
}


double MultispeciesCoalescentMigration::computeLnProbability( void )
{
    
    const Tree &sp = species_tree->getValue();
    const std::vector< TopologyNode* >& species_tree_nodes = sp.getNodes();
    const std::vector< TopologyNode* >& individual_tree_nodes = this->getValue().getNodes();
    
    size_t num_populations = sp.getNumberOfTips()*2-1;
    size_t num_individuals = num_taxa*2-1;
    
    // first let's create a map from species names to the nodes of the species tree
    std::map<std::string, TopologyNode * > species_names_2_species_nodes;
    for (std::vector< TopologyNode *>::const_iterator it = species_tree_nodes.begin(); it != species_tree_nodes.end(); ++it)
    {
        if ( (*it)->isTip() == true )
        {
            const std::string &name = (*it)->getName();
            species_names_2_species_nodes[name] = *it;
        }
    }
    
    // we need to get all the species ages
    std::multimap<double, TopologyNode * > species_ages_2_nodes;
    for (std::vector< TopologyNode *>::const_iterator it = species_tree_nodes.begin(); it != species_tree_nodes.end(); ++it)
    {
        TopologyNode *the_node = *it;
        if ( the_node->isTip() == false )
        {
            double age = the_node->getAge();
            species_ages_2_nodes.insert( std::pair<double, TopologyNode *>(age,the_node) );
            
        }
    }
    
    // create a vector of probabilities
    // the first #pop * #ind are the probs that an individual i is in population j, which are all initialized to 0.0
    std::vector<double> probabilities = std::vector<double>(num_populations*num_individuals + 1, 0.0);
    // the final probability is the probability of no coalescent event, which is at present time P(no coal|t=0) = 1.0
    probabilities[num_populations*num_individuals] = 1.0;
    
    // we need to get all the species ages
    std::multimap<double, TopologyNode * > individual_ages_2_nodes;
//    std::set<TopologyNode * > current_individuals;
    for (std::vector< TopologyNode *>::const_iterator it = individual_tree_nodes.begin(); it != individual_tree_nodes.end(); ++it)
    {
        TopologyNode *the_individual_node = *it;
        if ( the_individual_node->isTip() == true )
        {
            // this node was a tip node
            
            // store the node into our set of currently active lineages
            // const std::string &name = the_node->getName();
//            current_individuals.insert( the_individual_node );
            
            // we need to set the initial probability of the individual being in population j
            size_t individual_node_index = the_individual_node->getIndex();
            const std::string &species_name = the_individual_node->getSpeciesName();
            TopologyNode *species_node = species_names_2_species_nodes.at(species_name);
            size_t species_node_index = species_node->getIndex();
            probabilities[num_populations*individual_node_index+species_node_index] = 1.0;

        }
        else
        {
            // the node represents a coalescent event
            // therefore, we store the node in our map which sorts the coalescent events
            double age = the_individual_node->getAge();
            individual_ages_2_nodes.insert( std::pair<double, TopologyNode *>(age,the_individual_node) );
        }
    }
    
    // get the parameters of the model (Q, theta, delta)
    std::vector<double> thetas = std::vector<double>(num_populations, 0.0);
    for (size_t i=0; i<num_populations; ++i)
    {
        thetas[i] = 1.0 / getNe( i );
    }
    const RateGenerator* migration_rate_matrix = &Q->getValue();
    double overall_migration_rate = delta->getValue();

    double ln_probability = 0.0;
    std::multimap<double, TopologyNode * >::const_iterator it_species_ages      = species_ages_2_nodes.begin();
    std::multimap<double, TopologyNode * >::const_iterator it_individual_ages   = individual_ages_2_nodes.begin();
    double next_speciation_age = RbConstants::Double::inf;
    double next_coalescent_age = RbConstants::Double::inf;
    double previous_event_age = 0.0;
    while ( it_individual_ages != individual_ages_2_nodes.end() ) // && it_species_ages != species_ages_2_nodes.end()
    {
        
        // get the next speciation age
        if ( it_species_ages != species_ages_2_nodes.end() )
        {
            // there is at least one more speciation event coming (back in the past), so we get it
            next_speciation_age = it_species_ages->first;
        }
        else
        {
            // there are no more speciation events, which means we are in the root branch
            next_speciation_age = RbConstants::Double::inf;
        }
        
        // get the next speciation age
        if ( it_individual_ages != individual_ages_2_nodes.end() )
        {
            // there is at least one more coalescent event, so we take it
            next_coalescent_age = it_individual_ages->first;
        }
        else
        {
            // there are no more coalescent events, which means that all lineages have already coalesced
            // TODO: we could probably stop the computation here, because the only thing a single lineage can do is to migrate.
            next_coalescent_age = RbConstants::Double::inf;
        }
        
        // TODO: this might now be necessary but I'm doing this to be sure
        // we need to reset the probability of no coalescent event
        probabilities[num_populations*num_individuals] = 1.0;
        
        double next_event_age = RbMath::min(next_coalescent_age, next_speciation_age);
        double dt = 1E-8;
        
        MultispeciesCoalescentMigrationODE ode = MultispeciesCoalescentMigrationODE(thetas, migration_rate_matrix, overall_migration_rate, num_individuals, num_populations);

        typedef boost::numeric::odeint::runge_kutta_dopri5< std::vector< double > > stepper_type;
        boost::numeric::odeint::integrate_adaptive( make_controlled( 1E-7, 1E-7, stepper_type() ), ode, probabilities, previous_event_age, next_event_age, dt );
        
        // add the probability that there was no coalescent event until then
        ln_probability += log( probabilities[num_populations*num_individuals] );
        
        // we need to reset the probability of no coalescent event
        probabilities[num_populations*num_individuals] = 1.0;
        
        // check if this was a coalescent event
        if ( next_coalescent_age < next_speciation_age )
        {
            // get the individual node that represent the coalescent event
            TopologyNode* this_individual = it_individual_ages->second;
            
            // get the index of the left and right coalescing individuals
            size_t left_coalescing_individual  = this_individual->getChild(0).getIndex();
            size_t right_coalescing_individual = this_individual->getChild(1).getIndex();
                       
            // now compute the probability of a coalescent event
            // this could have happened in any population, so we need to integrate over all possible populations
            // and take the product of the probabilities that both indvididuals were present in that population
            double prob_coalescent = 0.0;
            for ( size_t pop_index=0; pop_index<num_populations; ++pop_index )
            {

                prob_coalescent += ( thetas[pop_index] *
                                     probabilities[num_populations*left_coalescing_individual+pop_index] *
                                     probabilities[num_populations*right_coalescing_individual+pop_index]);
                
            }
            ln_probability += log( prob_coalescent );
            
            for ( size_t pop_index=0; pop_index<num_populations; ++pop_index )
            {
                
                double prob_coal_this_pop = ( thetas[pop_index] *
                                              probabilities[num_populations*left_coalescing_individual+pop_index] *
                                              probabilities[num_populations*right_coalescing_individual+pop_index]);
                probabilities[num_populations*this_individual->getIndex()+pop_index] = prob_coal_this_pop / prob_coalescent;
                probabilities[num_populations*left_coalescing_individual+pop_index] = 0.0;
                probabilities[num_populations*right_coalescing_individual+pop_index] = 0.0;
                
            }
            
            // move the iterater of the individual ages
            ++it_individual_ages;
        }
        else if ( RbMath::isFinite(next_speciation_age) )
        {
            // this was a speciation event
            // we need to reorder the probabilities and lineages
            
            // get the individual node that represent the coalescent event
            TopologyNode* this_species = it_species_ages->second;
            
            
            // get the index of the left and right coalescing individuals
            size_t left_coalescing_pop  = this_species->getChild(0).getIndex();
            size_t right_coalescing_pop = this_species->getChild(1).getIndex();
            size_t new_pop_index        = this_species->getIndex();
            
            // "move" all the individuals into the new population
            // that means, we need to add the probabilities of the two descendant populations
            for ( size_t ind_index=0; ind_index<num_individuals; ++ind_index )
            {
                probabilities[num_populations*ind_index+new_pop_index] = probabilities[num_populations*ind_index+left_coalescing_pop] + probabilities[num_populations*ind_index+right_coalescing_pop];
                // to be safe, we set all the old probabilities to 0
                probabilities[num_populations*ind_index+left_coalescing_pop]  = 0.0;
                probabilities[num_populations*ind_index+right_coalescing_pop] = 0.0;
            }
            
            // move the iterater of the speciation ages
            ++it_species_ages;
        }
        
        // update the event age
        previous_event_age = next_event_age;
        
    }
    
    return ln_probability; // + logTreeTopologyProb;
    
}


double MultispeciesCoalescentMigration::drawNe( size_t index )
{
    
    return getNe( index );
}



double MultispeciesCoalescentMigration::getNe(size_t index) const
{
    
    if ( Ne != NULL )
    {
        return Ne->getValue();
    }
    else if ( Nes != NULL )
    {
        if ( index >= Nes->getValue().size() )
        {
            exit(-1);
        }
        return Nes->getValue()[index];
    }
    else
    {
        exit(-1);
    }
}


double MultispeciesCoalescentMigration::getMigrationRate(size_t from, size_t to) const
{
    return Q->getValue().getRate(from, to, 0, delta->getValue());
}


void MultispeciesCoalescentMigration::redrawValue( void )
{
    
    simulateTree();
    
}


void MultispeciesCoalescentMigration::resetTipAllocations( void )
{
    
    const Tree &sp = species_tree->getValue();
    const std::vector< TopologyNode* > &species_tree_nodes = sp.getNodes();
    
    // first let's create a map from species names to the nodes of the species tree
    std::map<std::string, TopologyNode * > species_names_2_species_nodes;
    for (std::vector< TopologyNode *>::const_iterator it = species_tree_nodes.begin(); it != species_tree_nodes.end(); ++it)
    {
        if ( (*it)->isTip() == true )
        {
            const std::string &name = (*it)->getName();
            species_names_2_species_nodes[name] = *it;
        }
    }
    
    // create a map for the individuals to branches
    individuals_per_branch = std::vector< std::set< const TopologyNode* > >(sp.getNumberOfNodes(), std::set< const TopologyNode* >() );
    for (size_t i=0; i<num_taxa; ++i)
    {
        //        const std::string &tip_name = it->getName();
        const TopologyNode &n = value->getNode( i );
        const std::string &species_name = n.getSpeciesName();
        TopologyNode *species_node = species_names_2_species_nodes.at(species_name);
        individuals_per_branch[ species_node->getIndex() ].insert( &n );
    }
    
    
}


void MultispeciesCoalescentMigration::setNes(TypedDagNode< RbVector<double> >* input_nes)
{
    
    removeParameter( Nes );
    removeParameter( Ne );
    
    Nes = input_nes;
    Ne  = NULL;
    
    addParameter( Nes );
    
}


void MultispeciesCoalescentMigration::setNe(TypedDagNode<double>* input_ne)
{
    
    removeParameter( Ne );
    removeParameter( Nes );
    
    Ne  = input_ne;
    Nes = NULL;
    
    addParameter( Ne );
}


/**
 * Set the current value.
 */
void MultispeciesCoalescentMigration::setValue(Tree *v, bool f )
{
    
    // delegate to super class
    TypedDistribution<Tree>::setValue(v, f);
    
    resetTipAllocations();
    
}




void MultispeciesCoalescentMigration::simulateTree( void )
{
    
    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    const Tree &sp = species_tree->getValue();
    const std::vector< TopologyNode* > &species_tree_nodes = sp.getNodes();
    // first let's create a map from species names to the nodes of the species tree
    std::map<std::string, TopologyNode * > species_names_2_nodes;
    std::multimap<double, TopologyNode * > species_ages_2_nodes;

    for (std::vector< TopologyNode *>::const_iterator it = species_tree_nodes.begin(); it != species_tree_nodes.end(); ++it)
    {
        TopologyNode *the_species_node = *it;
        if ( the_species_node->isTip() == true )
        {
            const std::string &name = the_species_node->getName();
            species_names_2_nodes[name] = the_species_node;
        }
        else
        {
            double age = the_species_node->getAge();
            species_ages_2_nodes.insert( std::pair<double, TopologyNode *>(age,the_species_node) );
            
        }
    }
    
    
    std::map< size_t, std::vector< TopologyNode* > > individuals_per_branch;
    std::map< TopologyNode*, size_t> individuals_2_population;
    std::vector<size_t> num_individuals_in_population = std::vector<size_t>(species_tree_nodes.size(), 0);

    for (std::vector< Taxon >::iterator it = taxa.begin(); it != taxa.end(); ++it)
    {
        TopologyNode *n = new TopologyNode( *it );
        const std::string &species_name = n->getSpeciesName();
        
        if ( species_name == "" )
        {
            throw RbException("Cannot match a taxon without species to a tip in the species tree. The taxon map is probably wrong.");
        }
        
        TopologyNode *species_node = species_names_2_nodes.at(species_name);
        
        if ( species_node == NULL )
        {
            throw RbException() << "Could not match a taxon with name" << species_name << " to any of the tips in the species tree.";
        }
        
        // initialize the age to be the present time
        n->setAge( 0.0 );
        
        // add the current individual to our map with a index to the population
        individuals_2_population.insert( std::pair<TopologyNode *, size_t>(n,species_node->getIndex()) );
        num_individuals_in_population[species_node->getIndex()]++;
        
//        individuals_per_branch;
    }
    
    TopologyNode *root = NULL;
    double current_time = 0.0;
    std::multimap<double, TopologyNode *>::const_iterator it = species_ages_2_nodes.begin();
    while ( individuals_2_population.size() > 1 )
    {
        
        TopologyNode *sp_node = NULL;
        double species_age = RbConstants::Double::inf;
        if ( it != species_ages_2_nodes.end() )
        {
            sp_node = it->second;
            species_age = it->first;
        }
        
        enum EVENT_TYPE { COALESCENT, MIGRATION, NO_EVENT };
        
        EVENT_TYPE next_event_type = NO_EVENT;
        size_t next_event_node = 0;
        double next_event_time = RbConstants::Double::inf;
        for (size_t this_population=0; this_population<species_tree_nodes.size(); ++this_population)
        {
            
            // first we simulate coalescent events
            size_t num_ind_this_pop = num_individuals_in_population[ this_population ];
            
            // we can only coalesce if there are more than 1 individuals
            if ( num_ind_this_pop > 1 )
            {
                double branch_ne = drawNe( this_population );
                double theta = 1.0 / branch_ne;

                double n_pairs = num_ind_this_pop * (num_ind_this_pop-1) / 2.0;
                double lambda = n_pairs * theta;
                double next_coalescent_time = current_time + RbStatistics::Exponential::rv( lambda, *rng);
            
                // check if the next coalescent time is smaller than the next event time
                // in that case, the next event could be a coalescent event
                if ( next_event_time > next_coalescent_time )
                {
                    next_event_time = next_coalescent_time;
                    next_event_type = COALESCENT;
                    next_event_node = this_population;
                }
            
            } // end-if there are >1 individuals in this population
            
            
            // compute the overall migration rate
            double migration_rate = 0.0;
            for (size_t alternative_population=0; alternative_population<species_tree_nodes.size()-1; ++alternative_population)
            {
                
                if ( this_population != alternative_population && this_population != species_tree_nodes.size()-1 )
                {
                    migration_rate += getMigrationRate(this_population, alternative_population);
                }
                
            }
            // we need to multiply the rate with the number of individuals
            migration_rate *= num_individuals_in_population[ this_population ];
            
            // draw the next migration event time
            double next_migration_time = current_time + RbStatistics::Exponential::rv( migration_rate, *rng);
            
            // if the next migration event time is small than the so far next event time
            // then this will be a migration event
            if ( next_event_time > next_migration_time )
            {
                next_event_time = next_migration_time;
                next_event_type = MIGRATION;
                next_event_node = this_population;
            }
            
        } // end-for over all populations
        
        
//        const TopologyNode *sp_parent_node = NULL;
////        double branch_length = RbConstants::Double::inf;
//        if ( sp_node->isRoot() == false )
//        {
//            sp_parent_node = &sp_node->getParent();
//            branch_length = sp_parent_node->getAge() - sp_node->getAge();
//        }
        
        if ( next_event_time < species_age  )
        {
            
            if ( next_event_type == COALESCENT )
            {
                // first we simulate coalescent events
                size_t num_ind_this_pop = num_individuals_in_population[ next_event_node ];

                std::vector<TopologyNode*> individuals_at_current_branch;
                for ( std::multimap<TopologyNode *, size_t>::const_iterator it = individuals_2_population.begin(); it != individuals_2_population.end(); ++it)
                {
                    if ( it->second == next_event_node )
                    {
                        individuals_at_current_branch.push_back( it->first );
                    }
                }
                // randomly coalesce two lineages
                size_t index = static_cast<size_t>( floor(rng->uniform01()*individuals_at_current_branch.size()) );
                TopologyNode *left = individuals_at_current_branch[index];
                individuals_at_current_branch.erase( individuals_at_current_branch.begin() + index);
            
                index = static_cast<size_t>( floor(rng->uniform01()*individuals_at_current_branch.size()) );
                TopologyNode *right = individuals_at_current_branch[index];
                individuals_at_current_branch.erase( individuals_at_current_branch.begin() + index);
            
                TopologyNode *new_parent = new TopologyNode();
                new_parent->addChild(left);
                left->setParent(new_parent);
                new_parent->addChild(right);
                right->setParent(new_parent);
            
                new_parent->setAge( next_event_time );
                root = new_parent;
            
//                individuals_at_current_branch.push_back( new_parent );
                individuals_2_population.insert( std::pair<TopologyNode *, size_t>(new_parent,next_event_node) );
                
                individuals_2_population.erase( left );
                individuals_2_population.erase( right );
                --num_individuals_in_population[ next_event_node ];
                
            }
            else if ( next_event_type == MIGRATION )
            {
                // second, we simulate migration events
                
                std::vector<TopologyNode*> individuals_at_current_branch;
                for ( std::multimap<TopologyNode *, size_t>::const_iterator it = individuals_2_population.begin(); it != individuals_2_population.end(); ++it)
                {
                    if ( it->second == next_event_node )
                    {
                        individuals_at_current_branch.push_back( it->first );
                    }
                }
                
                size_t this_population = next_event_node;
                
                double total_migration_rate = 0.0;
                for (size_t alternative_population=0; alternative_population<individuals_per_branch.size(); ++alternative_population)
                {
                    
                    if ( this_population != alternative_population )
                    {
                        total_migration_rate += getMigrationRate(this_population, alternative_population);
                    }
                    
                }
                
                size_t target_population = 0;
                double remaining_migration_rate = total_migration_rate * rng->uniform01();
                for (size_t alternative_population=0; alternative_population<individuals_per_branch.size(); ++alternative_population)
                {
                    
                    if ( this_population != alternative_population )
                    {
                        remaining_migration_rate -= getMigrationRate(this_population, alternative_population);
                    }
                    
                    if ( remaining_migration_rate < 0 )
                    {
                        target_population = alternative_population;
                        break;
                    }
                    
                }
                
                // randomly choose a lineage for migration
                size_t index = static_cast<size_t>( floor(rng->uniform01()*individuals_at_current_branch.size()) );
                TopologyNode *individual = individuals_at_current_branch[index];
                
                --num_individuals_in_population[this_population];
                ++num_individuals_in_population[target_population];
                individuals_2_population[ individual ] = target_population;

            }
            else
            {
                throw RbException("No such event.");
            }
            
            
        }
        else if ( sp_node != NULL )
        {
            size_t new_pop      = sp_node->getIndex();
            size_t left_pop     = sp_node->getChild(0).getIndex();
            size_t right_pop    = sp_node->getChild(1).getIndex();
            
            for ( std::multimap<TopologyNode *, size_t>::const_iterator it = individuals_2_population.begin(); it != individuals_2_population.end(); ++it)
            {
                if ( it->second == left_pop || it->second == right_pop )
                {
                    individuals_2_population[it->first] = new_pop;
                }
            }
            
            num_individuals_in_population[ new_pop ]    = num_individuals_in_population[ left_pop ] + num_individuals_in_population[ right_pop ];
            num_individuals_in_population[ left_pop ]   = 0;
            num_individuals_in_population[ right_pop ]  = 0;
            
            ++it;
            next_event_time = species_age;
        }
        
        current_time = next_event_time;
    }
    
    // the time tree object (topology + times)
    Tree *psi = new Tree();
    
    // internally we treat unrooted topologies the same as rooted
    psi->setRooted( true );
    
    // initialize the topology by setting the root
    psi->setRoot(root, true);
    
//    for ( std::map<TopologyNode*, double>::iterator it = nodes_2_ages.begin(); it != nodes_2_ages.end(); ++it)
//    {
//        TopologyNode *node = it->first;
//        node->setAge( it->second );
//    }
        
    // finally store the new value
    value = psi;
    
    resetTipAllocations();
}


/** Swap a parameter of the distribution */
void MultispeciesCoalescentMigration::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if ( oldP == species_tree )
    {
        species_tree = static_cast<const TypedDagNode< Tree >* >( newP );
    }
    
    if ( oldP == Nes )
    {
        Nes = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    
    if ( oldP == Ne )
    {
        Ne = static_cast<const TypedDagNode< double >* >( newP );
    }
    
    if ( oldP == Q )
    {
        Q = static_cast<const TypedDagNode< RateGenerator >* >( newP );
    }
    
    if ( oldP == delta )
    {
        delta = static_cast<const TypedDagNode< double >* >( newP );
    }
    
}

