#include "Clade.h"
#include "ConstantNode.h"
#include "MultispeciesCoalescentMigration.h"
#include "MultispeciesCoalescentMigrationODE.h"
#include "DistributionExponential.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbMathCombinatorialFunctions.h"
#include "TopologyNode.h"

#include <algorithm>
#include <cmath>

using namespace RevBayesCore;

MultispeciesCoalescentMigration::MultispeciesCoalescentMigration(const TypedDagNode<Tree> *sp, const std::vector<Taxon> &t) : TypedDistribution<Tree>( NULL ),
    taxa(t),
    species_tree( sp ),
    num_taxa( taxa.size() ),
    log_tree_topology_prob (0.0)
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( species_tree );
    
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
        TopologyNode* leftChild = new TopologyNode(0);
        parent->addChild(leftChild);
        leftChild->setParent(parent);
        tips.push_back(leftChild);
        
        // add a right child
        TopologyNode* rightChild = new TopologyNode(0);
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
    const std::vector< TopologyNode* > &species_tree_nodes = sp.getNodes();
    const std::vector< TopologyNode* > &individual_tree_nodes = this->getValue().getNodes();
    
    size_t num_populations = sp.getNumberOfTips();
    size_t num_individuals = num_taxa;
    
    // first let's create a map from species names to the nodes of the species tree
    std::map<std::string, TopologyNode * > species_names_2_species_nodes;
    for (std::vector< TopologyNode *>::const_iterator it = species_tree_nodes.begin(); it != species_tree_nodes.end(); ++it)
    {
        if ( (*it)->isTip() )
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
        if ( the_node->isTip() == true )
        {
//            const std::string &name = the_node->getName();
//            species_names_2_nodes[name] = the_node;
        }
        else
        {
            double age = the_node->getAge();
            species_ages_2_nodes.insert( std::pair<double, TopologyNode *>(age,the_node) );
            
        }
    }
    
    // create a vector of probabilities
    std::vector<double> probabilities = std::vector<double>(num_populations*num_individuals + 1, 0.0);
    probabilities[num_populations*num_individuals] = 1.0;
    
    // we need to get all the species ages
    std::multimap<double, TopologyNode * > individual_ages_2_nodes;
    std::set<TopologyNode * > current_individuals;
    for (std::vector< TopologyNode *>::const_iterator it = individual_tree_nodes.begin(); it != individual_tree_nodes.end(); ++it)
    {
        TopologyNode *the_node = *it;
        if ( the_node->isTip() == true )
        {
            // const std::string &name = the_node->getName();
            current_individuals.insert( the_node );
            
            
            // we also need to set the initial probability of the individual being in population j
            size_t node_index = the_node->getIndex();
            const std::string &species_name = the_node->getSpeciesName();
            TopologyNode *species_node = species_names_2_species_nodes[species_name];
            size_t species_node_index = species_node->getIndex();
            probabilities[num_populations*node_index+species_node_index] = 1.0;

        }
        else
        {
            double age = the_node->getAge();
            individual_ages_2_nodes.insert( std::pair<double, TopologyNode *>(age,the_node) );
            
        }
    }

    double ln_probability = 0.0;
    // double current_age = 0.0;
    std::multimap<double, TopologyNode * >::const_iterator it_species_ages      = species_ages_2_nodes.begin();
    std::multimap<double, TopologyNode * >::const_iterator it_individual_ages   = individual_ages_2_nodes.begin();
    double next_speciation_age = RbConstants::Double::inf;
    double next_coalescent_age = RbConstants::Double::inf;
    while ( current_individuals.empty() == false )
    {
        
        // get the next speciation age
        if ( it_species_ages != species_ages_2_nodes.end() )
        {
            next_speciation_age = it_species_ages->first;
        }
        else
        {
            next_speciation_age = RbConstants::Double::inf;
        }
        
        // get the next speciation age
        if ( it_individual_ages != individual_ages_2_nodes.end() )
        {
            next_coalescent_age = it_individual_ages->first;
        }
        else
        {
            next_coalescent_age = RbConstants::Double::inf;
        }
        
//        MultispeciesCoalescentMigrationODE ode = MultispeciesCoalescentMigrationODE(<#const std::vector<double> &t#>, <#const RevBayesCore::RateGenerator *q#>, <#double r#>, <#size_t n_ind#>, <#size_t n_pop#>)
    }
    
    // create a map for the individuals to branches
    individuals_per_branch = std::vector< std::set< const TopologyNode* > >(sp.getNumberOfNodes(), std::set< const TopologyNode* >() );
    for (size_t i=0; i<num_taxa; ++i)
    {
        //        const std::string &tip_name = it->getName();
        const TopologyNode &n = value->getNode( i );
        const std::string &species_name = n.getSpeciesName();
        TopologyNode *species_node = species_names_2_species_nodes[species_name];
        individuals_per_branch[ species_node->getIndex() ].insert( &n );
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
        return Nes->getValue()[index];
    }
    else
    {
        std::cerr << "Error: Null Pointers for Ne and Nes." << std::endl;
        exit(-1);
    }
}


double MultispeciesCoalescentMigration::getMigrationRate(size_t from, size_t to) const
{
    return Q->getValue().getRate(from, to, 0, delta->getValue());
}


//double MultispeciesCoalescentMigration::recursivelyComputeLnProbability( const RevBayesCore::TopologyNode &species_node )
//{
//    
//    double ln_prob_coal = 0;
//    
//    if ( species_node.isTip() == false )
//    {
//        individuals_per_branch[ species_node.getIndex() ].clear();
//        
//        for (size_t i=0; i<species_node.getNumberOfChildren(); ++i)
//        {
//            ln_prob_coal += recursivelyComputeLnProbability( species_node.getChild(i) );
//        }
//    }
//    
//    //    const TopologyNode *species_parent_node = NULL;
//    double species_age = species_node.getAge();
//    double parent_species_age = RbConstants::Double::inf;
//    
//    //        double branchLength = RbConstants::Double::inf;
//    if ( species_node.isRoot() == false )
//    {
//        const TopologyNode &species_parent_node = species_node.getParent();
//        parent_species_age = species_parent_node.getAge();
//    }
//    
//    // create a local copy of the individuals per branch
//    const std::set<const TopologyNode*> &initial_individuals = individuals_per_branch[species_node.getIndex()];
//    std::set<const TopologyNode*> remaining_individuals = initial_individuals;
//    
//    // get all coalescent events among the individuals
//    std::vector<double> coal_times;
//    
//    
//    std::map<double, const TopologyNode *> coal_times_2_nodes;
//    for ( std::set<const TopologyNode*>::iterator it = remaining_individuals.begin(); it != remaining_individuals.end(); ++it)
//    {
//        const TopologyNode *ind = (*it);
//        if ( ind->isRoot() == false )
//        {
//            const TopologyNode &parent = ind->getParent();
//            double parent_age = parent.getAge();
//            coal_times_2_nodes[ parent_age ] = &parent;
//        }
//    }
//    
//    double current_time = species_age;
//    while ( current_time < parent_species_age && coal_times_2_nodes.size() > 0 )
//    {
//        
//        const TopologyNode *parent = coal_times_2_nodes.begin()->second;
//        double parent_age = parent->getAge();
//        current_time = parent_age;
//        
//        if ( parent_age < parent_species_age )
//        { //Coalescence in the species tree branch
//            
//            // get the left and right child of the parent
//            const TopologyNode *left = &parent->getChild( 0 );
//            const TopologyNode *right = &parent->getChild( 1 );
//            if ( remaining_individuals.find( left ) == remaining_individuals.end() || remaining_individuals.find( right ) == remaining_individuals.end() )
//            {
//                // one of the children does not belong to this species tree branch
//                return RbConstants::Double::neginf;
//            }
//            
//            //We remove the coalescent event and the coalesced lineages
//            coal_times_2_nodes.erase( coal_times_2_nodes.begin() );
//            remaining_individuals.erase( remaining_individuals.find( left ) );
//            remaining_individuals.erase( remaining_individuals.find( right ) );
//            
//            //We insert the parent in the vector of lineages in this branch
//            remaining_individuals.insert( parent );
//            if ( parent->isRoot() == false )
//            {
//                const TopologyNode *grand_parent = &parent->getParent();
//                coal_times_2_nodes[ grand_parent->getAge() ] = grand_parent;
//            }
//            
//            coal_times.push_back( parent_age );
//            
//            
//        } //End if coalescence in the species tree branch
//        else
//        { //No more coalescences in this species tree branch
//            
//            // jump out of the while loop
//            //                currentTime = speciesAge;
//            break;
//        }
//        
//        
//    } // end of while loop
//    
//    if ( initial_individuals.size() > 1 )
//    {
////        ln_prob_coal += computeLnCoalescentProbability(initial_individuals.size(), coal_times, species_age, parent_species_age, species_node.getIndex(), species_node.isRoot() == false);
//    }
//    
//    
//    // merge the two sets of individuals that go into the next species
//    if ( species_node.isRoot() == false )
//    {
//        std::set<const TopologyNode *> &incoming_lineages = individuals_per_branch[ species_node.getParent().getIndex() ];
//        incoming_lineages.insert( remaining_individuals.begin(), remaining_individuals.end());
//    }
//    
//    return ln_prob_coal;
//}


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
        if ( (*it)->isTip() )
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
        TopologyNode *species_node = species_names_2_species_nodes[species_name];
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
        TopologyNode *the_node = *it;
        if ( the_node->isTip() == true )
        {
            const std::string &name = the_node->getName();
            species_names_2_nodes[name] = the_node;
        }
        else
        {
            double age = the_node->getAge();
            species_ages_2_nodes.insert( std::pair<double, TopologyNode *>(age,the_node) );
            
        }
    }
    
    
    std::map< size_t, std::vector< TopologyNode* > > individuals_per_branch;
    
    for (std::vector< Taxon >::iterator it = taxa.begin(); it != taxa.end(); ++it)
    {
        TopologyNode *n = new TopologyNode( *it );
        const std::string &species_name = n->getSpeciesName();
        
        if ( species_name == "" )
        {
            throw RbException("Cannot match a taxon without species to a tip in the species tree. The taxon map is probably wrong.");
        }
        
        TopologyNode *species_node = species_names_2_nodes[species_name];
        
        if ( species_node == NULL )
        {
            throw RbException("Could not match a taxon with name" + species_name + " to any of the tips in the species tree.");
        }
        
        n->setAge( 0.0 );
        std::vector< TopologyNode * > &nodes_at_node = individuals_per_branch[ species_node->getIndex() ];
        nodes_at_node.push_back( n );
    }
    
//    std::map<TopologyNode *, double> nodes_2_ages;
    TopologyNode *root = NULL;
    double current_time = 0.0;
    // we loop over the nodes of the species tree in phylogenetic traversal
    for (std::multimap<double, TopologyNode *>::const_iterator it = species_ages_2_nodes.begin(); it != species_ages_2_nodes.end(); ++it)
    {
        TopologyNode *sp_node = it->second;
        double species_age = it->first;
        
        enum EVENT_TYPE { COALESCENT, MIGRATION, NO_EVENT };
        
        EVENT_TYPE next_event_type = NO_EVENT;
        size_t next_event_node;
        double next_event_time = RbConstants::Double::inf;
        for (size_t this_population=0; this_population<individuals_per_branch.size(); ++this_population)
        {
            
            // first we simulate coalescent events
            const std::vector<TopologyNode*> &initial_individuals_at_branch = individuals_per_branch[ this_population ];
            double branch_ne = drawNe( this_population );
            
            double theta = 1.0 / branch_ne;
            size_t j = initial_individuals_at_branch.size();
            double n_pairs = j * (j-1) / 2.0;
            double lambda = n_pairs * theta;
            double next_coalescent_time = current_time + RbStatistics::Exponential::rv( lambda, *rng);
            
            if ( next_event_time < next_coalescent_time )
            {
                next_event_time = next_coalescent_time;
                next_event_type = COALESCENT;
                next_event_node = this_population;
            }

            double migration_rate = 0.0;
            for (size_t alternative_population=0; alternative_population<individuals_per_branch.size(); ++alternative_population)
            {
                
                if ( this_population != alternative_population )
                {
                    migration_rate += getMigrationRate(this_population, alternative_population);
                }
                
            }
            migration_rate *= initial_individuals_at_branch.size();
            double next_migration_time = current_time + RbStatistics::Exponential::rv( migration_rate, *rng);
            
            if ( next_event_time < next_migration_time )
            {
                next_event_time = next_migration_time;
                next_event_type = MIGRATION;
                next_event_node = this_population;
            }
            
        }
        
        
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
                std::vector<TopologyNode*> &individuals_at_current_branch = individuals_per_branch[ next_event_node ];

                
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
            
                root = new_parent;
            
                individuals_at_current_branch.push_back( new_parent );
                
                new_parent->setAge( next_event_time );
            
//                nodes_2_ages[new_parent] = next_event_time + sp_node->getAge();
                
            }
            else if ( next_event_type == MIGRATION )
            {
                // first we simulate coalescent events
                std::vector<TopologyNode*> &individuals_at_current_branch = individuals_per_branch[ next_event_node ];
                
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
                individuals_at_current_branch.erase( individuals_at_current_branch.begin() + index);
                
                individuals_per_branch[ target_population ].push_back( individual );

            }
            else
            {
                throw RbException("No such event.");
            }
            
            
        }
        else if ( sp_node != NULL )
        {
            std::vector<TopologyNode *> &incoming_lineages = individuals_per_branch[sp_node->getIndex()];
            const std::vector<TopologyNode*> &initial_individuals_at_left_branch  = individuals_per_branch[ sp_node->getChild(0).getIndex() ];
            const std::vector<TopologyNode*> &initial_individuals_at_right_branch = individuals_per_branch[ sp_node->getChild(1).getIndex() ];

            incoming_lineages.insert(incoming_lineages.end(), initial_individuals_at_left_branch.begin(), initial_individuals_at_left_branch.end());
            incoming_lineages.insert(incoming_lineages.end(), initial_individuals_at_right_branch.begin(), initial_individuals_at_right_branch.end());
            
            individuals_per_branch[ sp_node->getChild(0).getIndex() ] = std::vector<TopologyNode*>();
            individuals_per_branch[ sp_node->getChild(1).getIndex() ] = std::vector<TopologyNode*>();
        }
        
        
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

