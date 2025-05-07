#include <cstddef>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "DuplicationLossProcess.h"
#include "DistributionExponential.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbMathCombinatorialFunctions.h"
#include "TopologyNode.h"
#include "RbException.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "Taxon.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;




DuplicationLossProcess::DuplicationLossProcess(const TypedDagNode<Tree> *it, const std::vector<Taxon> &t) : TypedDistribution<Tree>( NULL ),
    taxa(t),
    individual_tree( it ),
    num_taxa( taxa.size() ),
    log_tree_topology_prob (0.0)
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( individual_tree );
    
    double ln_fact = RbMath::lnFactorial((int)(num_taxa));
    
    log_tree_topology_prob = (num_taxa - 1) * RbConstants::LN2 - 2.0 * ln_fact - std::log( num_taxa ) ;
    
    redrawValue();
    
}


DuplicationLossProcess::~DuplicationLossProcess()
{
    
}



void DuplicationLossProcess::attachTimes(Tree *psi, std::vector<TopologyNode *> &tips, size_t index, const std::vector<double> &times)
{
    
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


void DuplicationLossProcess::buildRandomBinaryTree(std::vector<TopologyNode*> &tips)
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


DuplicationLossProcess* DuplicationLossProcess::clone(void) const
{
    
    return new DuplicationLossProcess(*this);
}

double DuplicationLossProcess::computeLnProbability( void )
{
    resetTipAllocations();
    
    // variable declarations and initialization
    double ln_prob_coal = 0;
    
    const Tree &it = individual_tree->getValue();
    
    ln_prob_coal = recursivelyComputeLnProbability( it.getRoot() );
    
    
    return ln_prob_coal; // + logTreeTopologyProb;
    
}


double DuplicationLossProcess::computeLnDuplicationLossProbability(size_t num_genes_recent, const std::vector<double> &dupl_ages, double age_recent, double age_ancient, const TopologyNode &node_individual, bool f)
{
    
    double ln_prob = 0.0;
    
    size_t index_individual = node_individual.getIndex();
    double dupl_rate = duplication_rate->getValue();
    
    double current_age = age_recent;
    double current_ext_prob = 0.0;
    
    // if this is a tip
    if ( node_individual.isTip() )
    {
        current_ext_prob = 1.0 - gene_sampling_probability->getValue()[index_individual];
    }
    else
    {
        size_t index_left = node_individual.getChild(0).getIndex();
        size_t index_right = node_individual.getChild(1).getIndex();
        
        double ext_prob_left  = extinction_probs[index_left];
        double ext_prob_right = extinction_probs[index_right];
        
        current_ext_prob = ext_prob_left * ext_prob_right;
    }

    for ( size_t i=0; i<dupl_ages.size(); ++i )
    {
        
        double this_dupl_age = dupl_ages[i];
        double dt = this_dupl_age - current_age;
        
        ln_prob += (num_genes_recent-i) * log( computeD(dt, current_ext_prob) );
        ln_prob += log( dupl_rate );
        ln_prob += log( num_genes_recent - i - 1 );
        
        current_age = this_dupl_age;
        current_ext_prob = computeE( dt, current_ext_prob );
        
    }
    
    // final branch segment before coalescent of individual tree
    // i.e., there are no more duplications
    double dt = age_ancient - current_age;
    ln_prob += (num_genes_recent-dupl_ages.size()) * log( computeD(dt, current_ext_prob) );

    extinction_probs[index_individual] = computeE( dt, current_ext_prob );
    
    return ln_prob;
}

double DuplicationLossProcess::computeD(double dt, double e)
{
    // TODO: Some more testing and cleaning up @Dominik
    double lambda = duplication_rate->getValue();
    double mu = loss_rate->getValue();

    double A = lambda - mu;
    double B = ( (1.0 - 2*e)*lambda+mu ) / A;
    
    double et = exp(-A * dt);
    
    double D = 4.0 * et;
    D /= (1.0+B+et*(1-B));
    D /= (1.0+B+et*(1-B));

    return D;
}

double DuplicationLossProcess::computeE(double dt, double e)
{
    // TODO: Some more testing and cleaning up @Dominik
    double lambda = duplication_rate->getValue();
    double mu = loss_rate->getValue();
    
    double A = lambda - mu;
    double B = ( (1.0 - 2*e)*lambda+mu ) / A;
    
    double et = exp(A * dt);
    
    double E = lambda + mu - A *(1.0+B-et*(1.0-B))/(1.0+B+et*(1.0-B));
    E /= (2*lambda);
    
    return E;
}



double DuplicationLossProcess::recursivelyComputeLnProbability( const RevBayesCore::TopologyNode &individual_node )
{
    
    double ln_prob_dupl_loss = 0;
    
    double individual_age = individual_node.getAge();
    double parent_individual_age = RbConstants::Double::inf;
    
    if ( individual_node.isRoot() == false )
    {
        const TopologyNode &individual_parent_node = individual_node.getParent();
        parent_individual_age = individual_parent_node.getAge();
    }
    
    if ( individual_node.isTip() == false )
    {
        genes_per_branch_recent[ individual_node.getIndex() ].clear();
        
        for (size_t i=0; i<individual_node.getNumberOfChildren(); ++i)
        {
            ln_prob_dupl_loss += recursivelyComputeLnProbability( individual_node.getChild(i) );
        }
    }
    else
    {
        size_t left_index  = individual_node.getChild(0).getIndex();
        size_t right_index = individual_node.getChild(1).getIndex();
        std::set< const TopologyNode* > &genes_for_this_individual = genes_per_branch_recent[ individual_node.getIndex() ];
        const std::set< const TopologyNode* > &genes_for_left_descendant  = genes_per_branch_ancient[ left_index ];
//        const std::set< const TopologyNode* > &genes_for_right_descendant = genes_per_branch_ancient[ right_index ];
        for (std::set<const TopologyNode*>::iterator it=genes_for_this_individual.begin(); it!=genes_for_this_individual.end(); ++it)
        {
            const TopologyNode *this_gene_node = *it;
            double this_age = this_gene_node->getParent().getAge();
            if ( fabs( this_age - individual_age) < 1E-10 )
            {
                // coalescent event
            }
            else
            {
                // check if this gene comes from the left or right descendant
                if ( genes_for_left_descendant.find(this_gene_node) != genes_for_left_descendant.end() )
                {
                    // found in the left descendant
                    // multiply with probability of loss in right descendant
                    ln_prob_dupl_loss += log(extinction_probs[right_index]);
                }
                else
                {
                    // not found in the left descendant
                    // multiply with probability of loss in left descendant
                    ln_prob_dupl_loss += log(extinction_probs[left_index]);
                }
                
            }
        }
        //
        
    }
    
    
    // create a local copy of the genes per branch
    const std::set<const TopologyNode*> &initial_genes = genes_per_branch_recent[individual_node.getIndex()];
    std::set<const TopologyNode*> remaining_genes = initial_genes;
    
    // get all duplication events among the genes
    std::vector<double> duplication_times;
    
    
    std::map<double, const TopologyNode *> dupl_times_2_nodes;
    for ( std::set<const TopologyNode*>::iterator it = remaining_genes.begin(); it != remaining_genes.end(); ++it)
    {
        const TopologyNode *ind = (*it);
        if ( ind->isRoot() == false )
        {
            const TopologyNode &parent = ind->getParent();
            double parent_age = parent.getAge();
            dupl_times_2_nodes[ parent_age ] = &parent;
        }
    }
    
    double current_time = individual_age;
    while ( current_time < parent_individual_age && dupl_times_2_nodes.size() > 0 )
    {
        
        const TopologyNode *parent = dupl_times_2_nodes.begin()->second;
        double parent_age = parent->getAge();
        current_time = parent_age;
        
        if ( parent_age < parent_individual_age )
        { //Duplication in the individual tree branch
            
            // get the left and right child of the parent
            const TopologyNode *left = &parent->getChild( 0 );
            const TopologyNode *right = &parent->getChild( 1 );
            if ( remaining_genes.find( left ) == remaining_genes.end() || remaining_genes.find( right ) == remaining_genes.end() )
            {
                // one of the children/genes does not belong to this individual tree branch
                return RbConstants::Double::neginf;
            }
            
            // We remove the coalescent event and the coalesced lineages
            dupl_times_2_nodes.erase( dupl_times_2_nodes.begin() );
            remaining_genes.erase( remaining_genes.find( left ) );
            remaining_genes.erase( remaining_genes.find( right ) );
            
            // We insert the parent in the vector of lineages in this branch
            remaining_genes.insert( parent );
            if ( parent->isRoot() == false )
            {
                const TopologyNode *grand_parent = &parent->getParent();
                dupl_times_2_nodes[ grand_parent->getAge() ] = grand_parent;
            }
            
            if ( fabs( parent_age - individual_age) > 1E-10 )
            {
                duplication_times.push_back( parent_age );
            }

            
            
        } //End if duplication in the individual tree branch
        else
        { //No more duplication in this individual tree branch
            
            // jump out of the while loop
            //                currentTime = speciesAge;
            break;
        }
        
        
    } // end of while loop
    
    if ( initial_genes.size() > 1 )
    {
        ln_prob_dupl_loss += computeLnDuplicationLossProbability(initial_genes.size(), duplication_times, individual_age, parent_individual_age, individual_node.getIndex(), individual_node.isRoot() == false);
    }
    
    
    // merge the two sets of individuals that go into the next species
    if ( individual_node.isRoot() == false )
    {
        std::set<const TopologyNode *> &genes = genes_per_branch_recent[ individual_node.getParent().getIndex() ];
        genes.insert( remaining_genes.begin(), remaining_genes.end());
    }

    
    std::set<const TopologyNode *> &outgoing_genes = genes_per_branch_ancient[ individual_node.getIndex() ];
    outgoing_genes.clear();
    outgoing_genes.insert( remaining_genes.begin(), remaining_genes.end());
    
    return ln_prob_dupl_loss;
}


void DuplicationLossProcess::redrawValue( void )
{
    
    simulateTree();
    
}


void DuplicationLossProcess::resetTipAllocations( void )
{
    
    const Tree &ind_tree = individual_tree->getValue();
    const std::vector< TopologyNode* > &individual_tree_nodes = ind_tree.getNodes();
    
    // first let's create a map from individual names to the nodes of the individual tree
    std::map<std::string, TopologyNode * > individual_names_2_individual_nodes;
    for (std::vector< TopologyNode *>::const_iterator it = individual_tree_nodes.begin(); it != individual_tree_nodes.end(); ++it)
    {
        TopologyNode *this_node = (*it);
        if ( this_node->isTip() )
        {
            const std::string &name = this_node->getName();
            individual_names_2_individual_nodes[name] = *it;
        }
    }
    
    // create a map for the individuals to branches
    genes_per_branch_recent  = std::vector< std::set< const TopologyNode* > >(ind_tree.getNumberOfNodes(), std::set< const TopologyNode* >() );
    genes_per_branch_ancient = std::vector< std::set< const TopologyNode* > >(ind_tree.getNumberOfNodes(), std::set< const TopologyNode* >() );
    for (size_t i=0; i<num_taxa; ++i)
    {
        const TopologyNode &n = value->getNode( i );
        const std::string &individual_name = n.getSpeciesName();
        TopologyNode *individual_node = individual_names_2_individual_nodes[individual_name];
        genes_per_branch_recent[ individual_node->getIndex() ].insert( &n );
    }
    
    
}


/**
 * Set the current value.
 */
void DuplicationLossProcess::setValue(Tree *v, bool f )
{
    
    // delegate to super class
    TypedDistribution<Tree>::setValue(v, f);
    
    resetTipAllocations();
    
}



void DuplicationLossProcess::setDuplicationRate(TypedDagNode<RbVector<double> >* d)
{
    
    removeParameter( duplication_rate );
    removeParameter( duplication_rates );

    duplication_rate  = NULL;
    duplication_rates = d;
    
    addParameter( duplication_rates );
    
}


void DuplicationLossProcess::setDuplicationRate(TypedDagNode<double>* d)
{
    
    removeParameter( duplication_rate );
    removeParameter( duplication_rates );
    
    duplication_rate  = d;
    duplication_rates = NULL;
    
    addParameter( duplication_rate );
}


void DuplicationLossProcess::setLossRate(TypedDagNode<RbVector<double> >* l)
{
    
    removeParameter( loss_rate );
    removeParameter( loss_rates );
    
    loss_rate  = NULL;
    loss_rates = l;
    
    addParameter( loss_rates );
}


void DuplicationLossProcess::setLossRate(TypedDagNode<double>* l)
{
    
    removeParameter( loss_rate );
    removeParameter( loss_rates );
    
    loss_rate = l;
    loss_rates = NULL;
    
    addParameter( loss_rate );

}


void DuplicationLossProcess::setGeneSamplingProbability(TypedDagNode<RbVector<double> >* s)
{
    
    removeParameter( gene_sampling_probability );
    
    gene_sampling_probability = s;
    
    addParameter( gene_sampling_probability );
    
}




void DuplicationLossProcess::simulateTree( void )
{
    
    // TODO: Need to implement a correct way to simulate the gene-tree under the process @Dominik
    
    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    const Tree &sp = individual_tree->getValue();
    const std::vector< TopologyNode* > &species_tree_nodes = sp.getNodes();
    // first let's create a map from species names to the nodes of the species tree
    std::map<std::string, TopologyNode * > species_names_2_nodes;
    for (std::vector< TopologyNode *>::const_iterator it = species_tree_nodes.begin(); it != species_tree_nodes.end(); ++it)
    {
        if ( (*it)->isTip() == true )
        {
            const std::string &name = (*it)->getName();
            species_names_2_nodes[name] = *it;
        }
    }
    
    
    std::map< const TopologyNode *, std::vector< TopologyNode* > > individuals_per_branch;
    
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
            throw RbException() << "Could not match a taxon with name" << species_name << " to any of the tips in the species tree.";
        }
        
        n->setAge( 0.0 );
        std::vector< TopologyNode * > &nodes_at_node = individuals_per_branch[ species_node ];
        nodes_at_node.push_back( n );
    }
    
    
    /*   for (std::map<std::string,std::string>::iterator it = gene2species.begin(); it != gene2species.end(); ++it)
     {
     const std::string &tipName = it->first;
     TopologyNode *n = new TopologyNode( tipName );
     const std::string &speciesName = it->second;
     TopologyNode *speciesNode = speciesNames2Nodes[speciesName];
     std::vector< TopologyNode * > &nodesAtNode = individualsPerBranch[ speciesNode ];
     nodesAtNode.push_back( n );
     }
     */
    
    std::map<TopologyNode *, double> nodes_2_ages;
    TopologyNode *root = NULL;
    // we loop over the nodes of the species tree in phylogenetic traversal
    for (std::vector<TopologyNode *>::const_iterator it = species_tree_nodes.begin(); it != species_tree_nodes.end(); ++it)
    {
        TopologyNode *sp_node = *it;
        const TopologyNode *sp_parent_node = NULL;
        double branch_length = RbConstants::Double::inf;
        if ( sp_node->isRoot() == false )
        {
            sp_parent_node = &sp_node->getParent();
            branch_length = sp_parent_node->getAge() - sp_node->getAge();
        }
        
        std::vector<TopologyNode*> initial_individuals_at_branch = individuals_per_branch[sp_node];
        double branch_ne = 1.0;
        
        double theta = 1.0 / branch_ne;
        
        double prev_coalescent_time = 0.0;
        
        size_t j = initial_individuals_at_branch.size();
        double n_pairs = j * (j-1) / 2.0;
        double lambda = n_pairs * theta;
        double u = RbStatistics::Exponential::rv( lambda, *rng);
        double next_coalescent_time = prev_coalescent_time + u;
        
        while ( next_coalescent_time < branch_length && j > 1 )
        {
            // randomly coalesce two lineages
            size_t index = static_cast<size_t>( floor(rng->uniform01()*initial_individuals_at_branch.size()) );
            TopologyNode *left = initial_individuals_at_branch[index];
            initial_individuals_at_branch.erase( initial_individuals_at_branch.begin() + index);
            
            index = static_cast<size_t>( floor(rng->uniform01()*initial_individuals_at_branch.size()) );
            TopologyNode *right = initial_individuals_at_branch[index];
            initial_individuals_at_branch.erase( initial_individuals_at_branch.begin() + index);
            
            TopologyNode *new_parent = new TopologyNode();
            new_parent->addChild(left);
            left->setParent(new_parent);
            new_parent->addChild(right);
            right->setParent(new_parent);
            
            root = new_parent;
            
            if ( root == NULL )
            {
                std::cerr << "Oh, the root is NULL :(" << std::endl;
            }
            
            initial_individuals_at_branch.push_back( new_parent );
            
            nodes_2_ages[new_parent] = next_coalescent_time + sp_node->getAge();
            
            
            prev_coalescent_time = next_coalescent_time;
            j--;
            n_pairs = j * (j-1) / 2.0;
            lambda = n_pairs * theta ;
            u = RbStatistics::Exponential::rv( lambda, *rng);
            next_coalescent_time = prev_coalescent_time + u;
        }
        
        if ( sp_parent_node != NULL )
        {
            std::vector<TopologyNode *> &incoming_lineages = individuals_per_branch[sp_parent_node];
            incoming_lineages.insert(incoming_lineages.end(), initial_individuals_at_branch.begin(), initial_individuals_at_branch.end());
        }
        
        
    }
    
    // the time tree object (topology + times)
    Tree *psi = new Tree();
    
    // internally we treat unrooted topologies the same as rooted
    psi->setRooted( true );
    
    // initialize the topology by setting the root
    psi->setRoot(root, true);
    
    for ( std::map<TopologyNode*, double>::iterator it = nodes_2_ages.begin(); it != nodes_2_ages.end(); ++it)
    {
        TopologyNode *node = it->first;
        node->setAge( it->second );
    }
    
    // finally store the new value
    value = psi;
    
    resetTipAllocations();
}


/** Swap a parameter of the distribution */
void DuplicationLossProcess::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if ( oldP == individual_tree )
    {
        individual_tree = static_cast<const TypedDagNode< Tree >* >( newP );
    }
    
    if ( oldP == duplication_rate )
    {
        duplication_rate = static_cast<const TypedDagNode< double >* >( newP );
    }
    
    if ( oldP == duplication_rates )
    {
        duplication_rates = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    
    if ( oldP == loss_rate )
    {
        loss_rate = static_cast<const TypedDagNode< double >* >( newP );
    }
    
    if ( oldP == loss_rates )
    {
        loss_rates = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    
    if ( oldP == gene_sampling_probability )
    {
        gene_sampling_probability = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    
}
