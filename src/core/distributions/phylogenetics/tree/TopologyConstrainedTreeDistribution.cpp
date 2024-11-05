#include <cstddef>
#include <algorithm>
#include <iosfwd>
#include <map>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "AbstractRootedTreeDistribution.h"
#include "UniformTopologyBranchLengthDistribution.h"
#include "Clade.h"
#include "TopologyConstrainedTreeDistribution.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "DistributionExponential.h"
#include "RbException.h"
#include "StochasticNode.h"
#include "Taxon.h"
#include "TopologyNode.h"
#include "RbBitSet.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "StringUtilities.h"
#include "TimeInterval.h"
#include "Tree.h"
#include "TreeChangeEventHandler.h"
#include "TreeChangeEventMessage.h"
#include "TreeUtilities.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { template <class valueType> class RbOrderedSet; }

using namespace RevBayesCore;


/**
 * Constructor.
 *
 * The constructor connects the parameters of the birth-death process (DAG structure)
 * and initializes the probability density by computing the combinatorial constant of the tree structure.
 *
 * \param[in]    c         Clade constraints.
 */
TopologyConstrainedTreeDistribution::TopologyConstrainedTreeDistribution(TypedDistribution<Tree>* base_dist,
                                                                          const std::vector<Clade> &c,
                                                                          Tree *t,
                                                                          std::int64_t age_check_precision) : TypedDistribution<Tree>( NULL ),
//    active_backbone_clades( base_dist->getValue().getNumberOfInteriorNodes(), RbBitSet() ),
    active_clades( base_dist->getValue().getNumberOfInteriorNodes(), RbBitSet() ),
    backbone_topology(NULL),
    backbone_topologies(NULL),
    base_distribution( base_dist ),
    dirty_nodes( base_dist->getValue().getNumberOfNodes(), true ),
    monophyly_constraints( c ),
    num_backbones( 0 ),
    use_multiple_backbones( false ),
    starting_tree( t ),
    rooting_known( false ),
    is_rooted( true )
{
    
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    
    // add the parameters of the base distribution
    const std::vector<const DagNode*>& pars = base_distribution->getParameters();
    for (std::vector<const DagNode*>::const_iterator it = pars.begin(); it != pars.end(); ++it)
    {
        this->addParameter( *it );
    }
    
    value = &base_distribution->getValue();
    
    // Are there any fossils in the starting tree?
    if (starting_tree != NULL)
    {
        std::vector<bool> fossils;
        for (size_t i = 0; i < t->getNumberOfTips(); ++i)
        {
            TopologyNode* node = &t->getNode(i);
            fossils.push_back( node->isFossil() );
        }
        
        bool no_fossil = std::none_of(fossils.begin(), fossils.end(), [](bool v) { return v; });
        
        if (!no_fossil)
        {
            delete value;
            
            AbstractRootedTreeDistribution* tree_base_distribution = dynamic_cast<AbstractRootedTreeDistribution*>( base_distribution );
            std::vector<Taxon> taxa = tree_base_distribution->getTaxa();
            
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
    }
    
    initializeBitSets();
    redrawValue( SimulationCondition::MCMC );
}


/**
 * Copy Constructor.
 *
 * The constructor connects the parameters of the birth-death process (DAG structure)
 * and initializes the probability density by computing the combinatorial constant of the tree structure.
 *
 * \param[in]    c         Clade constraints.
 */
TopologyConstrainedTreeDistribution::TopologyConstrainedTreeDistribution(const TopologyConstrainedTreeDistribution &d) : TypedDistribution<Tree>( d ),
    active_backbone_clades( d.active_backbone_clades ),
    active_clades( d.active_clades ),
    backbone_constraints( d.backbone_constraints ),
    backbone_mask( d.backbone_mask ),
    backbone_topology( d.backbone_topology ),
    backbone_topologies( d.backbone_topologies ),
    base_distribution( d.base_distribution->clone() ),
    dirty_nodes( d.dirty_nodes ),
    monophyly_constraints( d.monophyly_constraints ),
    stored_backbone_clades( d.stored_backbone_clades ),
    stored_clades( d.stored_clades ),
    num_backbones( d.num_backbones ),
    use_multiple_backbones( d.use_multiple_backbones ),
    starting_tree( (d.starting_tree==NULL ? NULL : d.starting_tree->clone()) ),
    rooting_known( d.rooting_known ),
    is_rooted( d.is_rooted )
{
    // the copy constructor of the TypedDistribution creates a new copy of the value
    // however, here we want to hold exactly the same value as the base-distribution
    // thus, we delete the newly created value
    value->getTreeChangeEventHandler().removeListener( this );
    delete value;
    
    // and then set it to the value of the base distribution
    value = &base_distribution->getValue();
    
    value->getTreeChangeEventHandler().addListener( this );
    
    
    // add the parameters of the base distribution
    const std::vector<const DagNode*>& pars = base_distribution->getParameters();
    for (std::vector<const DagNode*>::const_iterator it = pars.begin(); it != pars.end(); ++it)
    {
        this->addParameter( *it );
    }
    
}



TopologyConstrainedTreeDistribution::~TopologyConstrainedTreeDistribution()
{
    
    delete base_distribution;
    
    //value->getTreeChangeEventHandler().removeListener( this );
    
    // DO NOT DELETE THE VALUE
    // the base distribution is the actual owner of the value!!!
    // we simply avoid the deletion of the value by setting its pointer to NULL
    // our base class, the TypedDistribution thinks that it owns the value and thus deletes it
    value = NULL;
    
    delete starting_tree;
}



TopologyConstrainedTreeDistribution& TopologyConstrainedTreeDistribution::operator=(const TopologyConstrainedTreeDistribution &d)
{
    
    if ( this != &d )
    {
        TypedDistribution<Tree>::operator=( d );
        
        delete base_distribution;
        delete starting_tree;
        
        active_backbone_clades          = d.active_backbone_clades;
        active_clades                   = d.active_clades;
        backbone_constraints            = d.backbone_constraints;
        backbone_mask                   = d.backbone_mask;
        backbone_topology               = d.backbone_topology;
        backbone_topologies             = d.backbone_topologies;
        base_distribution               = d.base_distribution->clone();
        dirty_nodes                     = d.dirty_nodes;
        monophyly_constraints           = d.monophyly_constraints;
        stored_backbone_clades          = d.stored_backbone_clades;
        stored_clades                   = d.stored_clades;
        num_backbones                   = d.num_backbones;
        use_multiple_backbones          = d.use_multiple_backbones;
        starting_tree                   = (d.starting_tree == NULL ? NULL : d.starting_tree->clone());
        rooting_known                   = d.rooting_known;
        is_rooted                       = d.is_rooted;

        // add the parameters of the base distribution
        const std::vector<const DagNode*>& pars = base_distribution->getParameters();
        for (std::vector<const DagNode*>::const_iterator it = pars.begin(); it != pars.end(); ++it)
        {
            this->addParameter( *it );
        }
        
    }
    
    return *this;
}


TopologyConstrainedTreeDistribution* TopologyConstrainedTreeDistribution::clone( void ) const
{
    
    return new TopologyConstrainedTreeDistribution( *this );
}


/**
 * Compute the log-transformed probability of the current value under the current parameter values.
 *
 */
double TopologyConstrainedTreeDistribution::computeLnProbability( void )
{
    recursivelyUpdateClades( value->getRoot() );
    
    // first check if the current tree matches the clade constraints
    if ( matchesConstraints() == false )
    {
        return RbConstants::Double::neginf;
    }
    
    if ( matchesBackbone() == false )
    {
        return RbConstants::Double::neginf;
    }
    
    double lnProb = base_distribution->computeLnProbability();
    
    return lnProb;
}


void TopologyConstrainedTreeDistribution::initializeBitSets(void)
{
    // fill the monophyly constraints bitsets
    for (size_t i = 0; i < monophyly_constraints.size(); i++)
    {
        // clade constraint has only one match
        if (monophyly_constraints[i].hasOptionalConstraints() == false)
        {
            RbBitSet b( value->getNumberOfTips() );
            for (size_t j = 0; j < monophyly_constraints[i].size(); j++)
            {
                const std::map<std::string, size_t> &taxon_map = value->getTaxonBitSetMap();
                const std::string &name = monophyly_constraints[i].getTaxonName(j);
                std::map<std::string, size_t>::const_iterator it = taxon_map.find( name );
                if ( it == taxon_map.end() )
                {
                    throw RbException("Could not find taxon with name '" + name + "'.");
                }
                size_t k = it->second;
                
                b.set(k);
            }
            monophyly_constraints[i].setBitRepresentation( b );
        }
        // clade constraint allows optional matches
        else
        {
            std::vector<Clade> optional_constraints = monophyly_constraints[i].getOptionalConstraints();
            for (size_t j = 0; j < optional_constraints.size(); j++)
            {
                RbBitSet b( value->getNumberOfTips() );
                for (size_t k = 0; k < optional_constraints[j].size(); k++)
                {
                    const std::map<std::string, size_t> &taxon_map = value->getTaxonBitSetMap();
                    const std::string &name = optional_constraints[j].getTaxonName(k);
                    std::map<std::string, size_t>::const_iterator it = taxon_map.find( name );
                    if ( it == taxon_map.end() )
                    {
                        throw RbException("Could not find taxon with name '" + name + "'.");
                    }
                    size_t s = it->second;
                    
                    b.set(s);
                }
                optional_constraints[j].setBitRepresentation( b );
            }
            monophyly_constraints[i].setOptionalConstraints( optional_constraints );
        }
        
    }
    
    // reset the backbone constraints and mask
    backbone_constraints.clear();
    backbone_mask.clear();
    backbone_constraints.resize(num_backbones);
    backbone_mask.resize( num_backbones );
    
    // add the backbone constraints
    if ( backbone_topologies != NULL && use_multiple_backbones )
    {
        for (size_t i = 0; i < num_backbones; i++)
        {
            backbone_mask[i] = RbBitSet( value->getNumberOfTips() );
            backbone_mask[i] |= recursivelyAddBackboneConstraints( backbone_topologies->getValue()[i].getRoot(), i );
        }
    }
    else if ( backbone_topology != NULL && !use_multiple_backbones )
    {
        backbone_mask[0] = RbBitSet( value->getNumberOfTips() );
        backbone_mask[0] |= recursivelyAddBackboneConstraints( backbone_topology->getValue().getRoot(), 0 );
    }
    
}


void TopologyConstrainedTreeDistribution::fireTreeChangeEvent(const TopologyNode &n, const unsigned& m)
{
    if (m == TreeChangeEventMessage::DEFAULT || m == TreeChangeEventMessage::TOPOLOGY)
    {
        
        recursivelyFlagNodesDirty(n);
    }
}


/**
 * Touch the current value and reset some internal flags.
 * If the root age variable has been restored, then we need to change the root age of the tree too.
 */
void TopologyConstrainedTreeDistribution::getAffected(RbOrderedSet<DagNode *> &affected, const DagNode *affecter)
{
    
    // delegate to the base distribution
    base_distribution->getAffected(affected, affecter);
}


/**
 * We check here if all the constraints are satisfied.
 * These are hard constraints, that is, the clades must be monophyletic.
 *
 * \return     True if the constraints are matched, false otherwise.
 */
bool TopologyConstrainedTreeDistribution::matchesBackbone( void )
{
    
    // ensure that each backbone constraint is found in the corresponding active_backbone_clades
    for (size_t i = 0; i < num_backbones; i++)
    {
        bool is_negative_constraint = false;
        if (backbone_topology != NULL)
        {
            is_negative_constraint = backbone_topology->getValue().isNegativeConstraint();
        }
        else if (backbone_topologies != NULL)
        {
            is_negative_constraint = ( backbone_topologies->getValue() )[i].isNegativeConstraint();
        }
        
        std::vector<bool> negative_constraint_found( backbone_constraints[i].size(), false );
        for (size_t j = 0; j < backbone_constraints[i].size(); j++)
        {
            std::vector<RbBitSet>::iterator it = std::find(active_backbone_clades[i].begin(), active_backbone_clades[i].end(), backbone_constraints[i][j] );
            
            // the search fails if the positive/negative backbone constraint is not satisfied
            if (it == active_backbone_clades[i].end() && !is_negative_constraint )
            {
                // match fails if positive constraint is not found
                return false;
            }
            else if (it != active_backbone_clades[i].end() && is_negative_constraint )
            {
                // match fails if negative constraint is found
                negative_constraint_found[j] = true;
            }
        }
        
        // match fails if all negative backbone clades are found
        bool negative_constraint_failure = true;
        for (size_t j = 0; j < negative_constraint_found.size(); j++)
        {
            if (negative_constraint_found[j] == false)
            {
                negative_constraint_failure = false;
            }
        }
        if (negative_constraint_failure)
        {
            return false;
        }
    }
    
    // if no search has failed, then the match succeeds
    return true;
}


/**
 * We check here if all the monophyly constraints are satisfied.
 *
 * \return     True if the constraints are matched, false otherwise.
 */
bool TopologyConstrainedTreeDistribution::matchesConstraints( void )
{
    for (size_t i = 0; i < monophyly_constraints.size(); i++)
    {
        
        std::vector<Clade> constraints;
        if ( monophyly_constraints[i].hasOptionalConstraints() == true )
        {
            constraints = monophyly_constraints[i].getOptionalConstraints();
        }
        else
        {
            constraints.push_back(monophyly_constraints[i]);
        }
        
        std::vector<bool> constraint_satisfied( constraints.size(), false );
        for (size_t j = 0; j < constraints.size(); j++)
        {
            
            std::vector<RbBitSet>::iterator it = std::find(active_clades.begin(), active_clades.end(), constraints[j].getBitRepresentation() );
            
            if (it != active_clades.end() && constraints[j].isNegativeConstraint() == false )
            {
                constraint_satisfied[j] = true;
            }
            else if (it == active_clades.end() && constraints[j].isNegativeConstraint() )
            {
                constraint_satisfied[j] = true;
            }
        }
        
        // match fails if no optional positive or negative constraints satisfied
        bool any_satisfied = false;
        for (size_t j = 0; j < constraint_satisfied.size(); j++)
        {
            if ( constraint_satisfied[j] == true )
            {
                any_satisfied = true;
                break;
            }
        }
        
        if ( any_satisfied == false )
        {
            return false;
        }
    }
    
    return true;
}


void TopologyConstrainedTreeDistribution::recursivelyFlagNodesDirty(const TopologyNode& n)
{
    
    dirty_nodes[ n.getIndex() ] = true;
    
    if ( n.isRoot() == false )
    {
        recursivelyFlagNodesDirty(n.getParent());
    }
    
}


RbBitSet TopologyConstrainedTreeDistribution::recursivelyAddBackboneConstraints( const TopologyNode& node, size_t backbone_idx )
{
    RbBitSet tmp( value->getNumberOfTips() );
    
    if ( node.isTip() )
    {
        const std::map<std::string, size_t>& taxon_map = value->getTaxonBitSetMap();
        const std::string& name = node.getName();
        std::map<std::string, size_t>::const_iterator it = taxon_map.find(name);
        if (it == taxon_map.end()) {
            
            throw RbException("Taxon named " + it->first + " not found in tree's taxon map!");
        }
        tmp.set( it->second );
    }
    else
    {
        // get the child names
        for (size_t i = 0; i < node.getNumberOfChildren(); i++)
        {
            tmp |= recursivelyAddBackboneConstraints( node.getChild(i), backbone_idx );
        }
        
        if ( node.isRoot() == false )
        {
            backbone_constraints[backbone_idx].push_back(tmp);
        }
    }
    
    return tmp;
}


RbBitSet TopologyConstrainedTreeDistribution::recursivelyUpdateClades( const TopologyNode& node )
{
    if ( node.isTip() )
    {
        RbBitSet tmp = RbBitSet( value->getNumberOfTips() );
        const std::map<std::string, size_t>& taxon_map = value->getTaxonBitSetMap();
        const std::string& name = node.getName();
        std::map<std::string, size_t>::const_iterator it = taxon_map.find(name);
        tmp.set( it->second );
        return tmp;
    }
    else if ( node.isRoot() )
    {
        if ( dirty_nodes[node.getIndex()] == true )
        {
            for (size_t i = 0; i < node.getNumberOfChildren(); i++)
            {
                recursivelyUpdateClades( node.getChild(i) );
            }
            
            dirty_nodes[node.getIndex()] = false;
        }
        
        return RbBitSet( value->getNumberOfTips() ).set();
    }
    else
    {
        if ( dirty_nodes[node.getIndex()] == true )
        {
            RbBitSet tmp = RbBitSet( value->getNumberOfTips() );
            for (size_t i = 0; i < node.getNumberOfChildren(); i++)
            {
                tmp |= recursivelyUpdateClades( node.getChild(i) );
            }
            
            // update the clade
            size_t idx = node.getIndex() - value->getNumberOfTips();
            active_clades[idx] = tmp;
            
            for (size_t i = 0; i < num_backbones; i++)
            {
                active_backbone_clades[i][idx] = tmp & backbone_mask[i];
            }
            
            
            dirty_nodes[node.getIndex()] = false;
        }
        
        return active_clades[node.getIndex() - value->getNumberOfTips()];
    }
}


/**
 * Redraw the current value. We delegate this to the simulate method.
 */
void TopologyConstrainedTreeDistribution::redrawValue( SimulationCondition c )
{
    
    Tree* new_value = NULL;
    bool alwaysReturn = (c == SimulationCondition::MCMC);
    
    if ( starting_tree == NULL )
    {
        if ( rooting_known == false )
        {
//            base_distribution->redrawValue();
//            is_rooted = base_distribution->getValue().isRooted();
            is_rooted = true;
            rooting_known = true;
            value = NULL;
        }
            
        if ( is_rooted == true )
        {
            new_value = simulateRootedTree(alwaysReturn);
        }
        else
        {
            new_value = simulateUnrootedTree();
        }
        // base_distribution->redrawValue();
    }
    else
    {
        new_value = starting_tree->clone();
    }
    
    if ( value != NULL )
    {
        value->getTreeChangeEventHandler().removeListener( this );
    }
    new_value->getTreeChangeEventHandler().addListener( this );
    
    // if we don't own the tree, then we just replace the current pointer with the pointer
    // to the new value of the base distribution
    value = new_value;
    base_distribution->setValue( value );
    
    // recompute the active clades
    dirty_nodes = std::vector<bool>( value->getNumberOfNodes(), true );
    active_clades = std::vector<RbBitSet>(value->getNumberOfInteriorNodes(), RbBitSet());

    recursivelyUpdateClades( value->getRoot() );
    
    stored_clades          = active_clades;
    stored_backbone_clades = active_backbone_clades;
}

void TopologyConstrainedTreeDistribution::redrawValue( void )
{
    // if no condition is specified, assume the most restrictive
    redrawValue(SimulationCondition::VALIDATION);
}


void TopologyConstrainedTreeDistribution::setBackbone(const TypedDagNode<Tree> *backbone_one, const TypedDagNode<RbVector<Tree> > *backbone_many)
{
    if (backbone_one == NULL && backbone_many == NULL)
    {
        ; // do nothing
    }
    else if (backbone_one != NULL && backbone_many != NULL)
    {
        ; // do nothing
    }
    else
    {
        
        
        // clear old parameter
        if (backbone_topology != NULL)
        {
            this->removeParameter( backbone_topology );
            backbone_topology = NULL;
        }
        else
        {
            this->removeParameter( backbone_topologies );
            backbone_topologies = NULL;
        }
        
        // set new parameter
        if (backbone_one != NULL)
        {
            backbone_topology = backbone_one;
            num_backbones = 1;
            use_multiple_backbones = false;
            this->addParameter( backbone_one );
        }
        else
        {
            backbone_topologies = backbone_many;
            num_backbones = backbone_topologies->getValue().size();
            use_multiple_backbones = true;
            this->addParameter( backbone_many );
        }
        
        for (size_t i = 0; i < num_backbones; i++)
        {
            std::vector<RbBitSet>v( base_distribution->getValue().getNumberOfInteriorNodes(), RbBitSet() );
            active_backbone_clades.push_back(v);
        }
        backbone_mask = std::vector<RbBitSet>( num_backbones, RbBitSet(base_distribution->getValue().getNumberOfInteriorNodes()) );
        
        
        initializeBitSets();
        
        // redraw the current value
        if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
        {
            this->redrawValue();
        }
        
    }
}

void checkCladesConsistent(const std::vector<Clade>& clades)
{
    // complain if we have conflicts
    for(int i=0;i<clades.size();i++)
    {
        auto& c1 = clades[i];
        if ( c1.isNegativeConstraint() ) continue;

        for(int j=0;j<clades.size();j++)
        {
            if (i == j) continue;

            auto& c2 = clades[j];
            if ( c2.isNegativeConstraint() ) continue;

            // c1 conflicts with c2 then say which clades conflict
            if (c1.conflicts(c2))
            {
                RbException e;
                e<<"Clades conflict:\n"
                 <<"  Clade1 = "<<c1<<"\n"
                 <<"  Clade2 = "<<c2<<"\n";

                // Here we print the intersection.
                // But we should probably print the smallest of:
                // (c1 AND c2), (c1 - c2), (c2 - c1).
                e<<"  Overlap = ";
                auto both = c1.intersection(c2);
                for(auto& taxon: both)
                    e<<taxon.getName()<<" ";

                throw e;
            }

            // If c2 is within c1 then it shouldn't be older
            if (c1.isNestedWithin(c2))
            {
                if (c2.getAge() > c1.getAge())
                    throw RbException()<<"Clade ages conflict:\n"
                                       <<"Clade 1 = "<<c1<<" Age = "<<c1.getAge()<<"\n"
                                       <<"Clade 2 = "<<c2<<" Age = "<<c2.getAge();
            }
        }
    }
}

/**
 *
 */
Tree* TopologyConstrainedTreeDistribution::simulateRootedTree( bool alwaysReturn )
{
        
    // the time tree object (topology & times)
    Tree *psi = new Tree();

    // internally we treat unrooted topologies the same as rooted
    psi->setRooted( true );

    AbstractRootedTreeDistribution* tree_base_distribution = dynamic_cast<AbstractRootedTreeDistribution*>( base_distribution );
    if ( tree_base_distribution == NULL )
    {
        throw RbException("dnConstrainedTopology cannot simulate from the base distribution. Use the 'initialTree' argument to provide a custom starting tree.");
    }
    size_t num_taxa = tree_base_distribution->getNumberOfTaxa();
    const std::vector<Taxon> &taxa = tree_base_distribution->getTaxa();

    // add a clounter variable of how many missing taxa we have already added
    size_t n_added_missing_taxa = 0;

    // create the tip nodes
    std::vector<TopologyNode*> nodes;
    for (size_t i=0; i<num_taxa; ++i)
    {
        // create the i-th taxon
        TopologyNode* node = new TopologyNode( taxa[i], i );

        // set the age of this tip node
        node->setAge( taxa[i].getAge() );

        // add the new node to the list
        nodes.push_back( node );
    }

    double ra = tree_base_distribution->getRootAge();
    double max_age = tree_base_distribution->getOriginAge();

    // we need a sorted vector of constraints, starting with the smallest
    std::vector<Clade> sorted_clades;

    for (size_t i=0; i<monophyly_constraints.size(); ++i)
    {
        Clade& monophyly_constraint = monophyly_constraints[i];
        if ( monophyly_constraint.getAge() > max_age )
        {
            throw RbException("Cannot simulate tree: clade constraints are older than the origin age.");
        }

        // set the ages of each of the taxa in the constraint
//        monophyly_constraint.setAges( taxa );

        // populate sorted clades vector
        if ( monophyly_constraint.size() > 1 && monophyly_constraint.size() < num_taxa )
        {
            if ( monophyly_constraint.hasOptionalConstraints() == true )
            {
                std::vector<Clade> optional_constraints = monophyly_constraint.getOptionalConstraints();
                size_t idx = (size_t)( GLOBAL_RNG->uniform01() * optional_constraints.size() );
                sorted_clades.push_back( optional_constraints[idx] );
            }
            else
            {
                sorted_clades.push_back( monophyly_constraint );
            }
        }

    }

    // create a clade that contains all species
    Clade all_species = Clade(taxa);
    all_species.setAge( ra );
    sorted_clades.push_back(all_species);

    // complain if we have conflicts
    checkCladesConsistent(sorted_clades);

    size_t num_clades = sorted_clades.size();
    std::sort(sorted_clades.begin(), sorted_clades.end(), cladeSmaller);

    /*
     * Walk clade constraints from tips to root.
     * For each clade, make a tree with other clades ("virtual taxa") or tip taxa as the leaves.
     * To find virtual taxa, walk back from the current taxon to the start of the sorted_clades list:
     *   If we find a nested clade, add it as a child.
     *   Remove its tip taxa from the unclaimed_taxa set to avoid adding grandchildren as children.
     */

    /*
     * An alternative idea would be to make a tree with all leaves as children of the root.
     * Then for each clade we try find the MRCA, and see if we can create a new child that
     * contains some of the children of the MCRA.  It should contain all the clade taxa and
     * no others.  If we cannot do that, then the clade is not consistent with the tree.
     *
     * For constraints with rooted splits that separate TAXA1 from TAXA2+root, we can use the
     * BUILD algorithm to find a consistent tree or report failure.
     *
     * However, perhaps a better idea would be to use MCMC to find a solution to the constraints.
     * There isn't a general algorithm for determining if tree constraints are consistent,
     * especially when we include NOT constraints and OR constraints.
     */
    
    std::vector<Clade> virtual_taxa;
    int i = -1;
    for (std::vector<Clade>::iterator it = sorted_clades.begin(); it != sorted_clades.end(); it++)
    {
        ++i;
        const Clade &c = *it;

        // ignore negative clade constraints during simulation
        if ( c.isNegativeConstraint() ) continue;

        std::vector<Taxon> unclaimed_taxa = c.getTaxa();
        std::vector<Clade> children;

        int j = i;
        std::vector<Clade>::reverse_iterator jt(it);
        for (; jt != sorted_clades.rend(); jt++)
        {
            j--;
            const Clade &c_nested = *jt;

            // ignore negative clade constraints during simulation
            if ( c_nested.isNegativeConstraint() ) continue;

            const std::vector<Taxon>& taxa_nested = c_nested.getTaxa();

            bool found_all = true;
            bool found_some = false;
            for (auto& taxon_nested: taxa_nested)
            {
                std::vector<Taxon>::iterator kt = std::find(unclaimed_taxa.begin(), unclaimed_taxa.end(), taxon_nested);
                if ( kt != unclaimed_taxa.end() )
                {
                    unclaimed_taxa.erase( kt );
                    found_some = true;
                }
                else
                {
                    found_all = false;
                }
            }

            if ( found_all == true )
            {
                //                c.addTaxon( virtual_taxa[j] );
                //                taxa.push_back( virtual_taxa[j] );
                children.push_back( virtual_taxa[j] );
            }

            // We check for conflicts before we try to construct the tree.
            // So any overlapping clades should be nested.
            assert(found_all == found_some);
        }

        std::vector<TopologyNode*> nodes_in_clade;

        for (size_t k=0; k<unclaimed_taxa.size(); ++k)
        {
            const Taxon& taxon = unclaimed_taxa[k];

            Clade tmp_clade = Clade( taxon );
            tmp_clade.setAge( taxon.getAge() );
            children.push_back( tmp_clade );
        }

        for (size_t k=0; k<children.size(); ++k)
        {
            const Clade& clade = children[k];

            for (size_t j = 0; j < nodes.size(); ++j)
            {
                if (nodes[j]->getClade() == clade)
                {
                    nodes_in_clade.push_back( nodes[j] );
                    nodes.erase( nodes.begin()+j );
                    break;
                }
            }
        }


        // here we need to start adding our "fake" tips
        for ( int index_missing_species = 0; index_missing_species < c.getNumberMissingTaxa(); ++index_missing_species)
        {
            ++n_added_missing_taxa;
            TopologyNode* missing_node = new TopologyNode("Missing_Taxon_" + ( c.getCladeName() != "" ? c.getCladeName() + "_" : "") + StringUtilities::to_string(n_added_missing_taxa) );
            missing_node->setAge( 0.0 );
            nodes_in_clade.push_back( missing_node );
        }

        double clade_age = c.getAge();

        double max_node_age = 0;
        for (size_t k=0; k<nodes_in_clade.size(); ++k)
        {
            const TopologyNode* node = nodes_in_clade[k];

            if ( node->getAge() > max_node_age )
            {
                max_node_age = node->getAge();
            }
        }

        if ( clade_age <= max_node_age )
        {
            // Get the rng
//            RandomNumberGenerator* rng = GLOBAL_RNG;

//            clade_age = rng->uniform01() * ( max_age - max_node_age ) + max_node_age;
            clade_age = tree_base_distribution->simulateCladeAge(nodes_in_clade.size(), max_age, 0, max_node_age, alwaysReturn);
        }

        tree_base_distribution->simulateClade(nodes_in_clade, clade_age, 0.0, alwaysReturn);
        nodes.push_back( nodes_in_clade[0] );

        std::vector<Taxon> v_taxa;
        nodes_in_clade[0]->getTaxa(v_taxa);
        Clade new_clade = Clade(v_taxa);
        new_clade.setAge( nodes_in_clade[0]->getAge() );
        virtual_taxa.push_back( new_clade );
    }

    TopologyNode *root = nodes[0];

    // initialize the topology by setting the root
    psi->setRoot(root, true);

    return psi;
}


/**
 *
 */
Tree* TopologyConstrainedTreeDistribution::simulateUnrootedTree( void )
{
    
    // the time tree object (topology & times)
    Tree *psi = new Tree();
    
    // internally we treat unrooted topologies the same as rooted
    psi->setRooted( false );
    
    UniformTopologyBranchLengthDistribution* tree_base_distribution = dynamic_cast<UniformTopologyBranchLengthDistribution*>( base_distribution );
    if ( tree_base_distribution == NULL )
    {
        throw RbException("dnConstrainedTopology cannot simulate from the base distribution. Use the 'initialTree' argument to provide a custom starting tree.");
    }
    const std::vector<Taxon> &taxa = tree_base_distribution->getTaxa();
    size_t num_taxa = taxa.size();
    
    // create the tip nodes
    std::vector<TopologyNode*> nodes;
    for (size_t i=0; i<num_taxa; ++i)
    {
        
        // create the i-th taxon
        TopologyNode* node = new TopologyNode( taxa[i], i );
        
        // add the new node to the list
        nodes.push_back( node );
        
    }
    
    if ( backbone_topology != NULL )
    {
        psi = backbone_topology->getValue().clone();
        std::vector<TopologyNode*> inserted_nodes = psi->getNodes();
        
        for (size_t i=0; i<num_taxa; ++i)
        {
            
            bool contains = false;
            for ( size_t j=0; j<inserted_nodes.size(); ++j )
            {
                if ( inserted_nodes[j]->getName() == taxa[i].getName() )
                {
                    contains = true;
                    break;
                }
            }
            
            if ( contains == false )
            {
                // randomly pick a branch to add
                size_t index = size_t(GLOBAL_RNG->uniform01() * inserted_nodes.size());
                while ( inserted_nodes[index]->isRoot() ) {
                    index = size_t(GLOBAL_RNG->uniform01() * inserted_nodes.size());
                }
                
                TopologyNode* new_node = new TopologyNode(taxa[i]);
                TopologyNode* tmp_node = new TopologyNode();
                
                tmp_node->addChild( new_node );
                new_node->setParent( tmp_node );
                
                TopologyNode* child = inserted_nodes[index];
                TopologyNode* parent = &child->getParent();
                
                if ( parent->getNumberOfChildren() > 3 )
                {
                    throw RbException("Wrong number of children.");
                }
                
                parent->removeChild( child );
                parent->addChild( tmp_node );
                tmp_node->setParent( parent );
                
                child->setParent( tmp_node );
                tmp_node->addChild( child );
                
                tmp_node->setBranchLength( child->getBranchLength() * 0.5 );
                child->setBranchLength( child->getBranchLength() * 0.5 );
                new_node->setBranchLength( RbStatistics::Exponential::rv(10.0, *GLOBAL_RNG) );
                
                inserted_nodes.push_back( tmp_node );
                inserted_nodes.push_back( new_node );
                
                if ( tmp_node->getNumberOfChildren() != 2 )
                {
                    throw RbException("Wrong number of children.");
                }
                if ( parent->getNumberOfChildren() > 3 )
                {
                    throw RbException("Wrong number of children.");
                }
            }
            
        }
        
        // initialize the topology by setting the root
        psi->setRoot(&psi->getRoot(), true);
        
        return psi;
    }
    
    // we need a sorted vector of constraints, starting with the smallest
    std::vector<Clade> sorted_clades;
    
    for (size_t i = 0; i < monophyly_constraints.size(); ++i)
    {
        
        // set ages for optional constraints
        if ( monophyly_constraints[i].hasOptionalConstraints() == true )
        {
            std::vector<Clade> optional_constraints = monophyly_constraints[i].getOptionalConstraints();
            monophyly_constraints[i].setOptionalConstraints( optional_constraints );
        }
        // populate sorted clades vector
        if ( monophyly_constraints[i].size() > 1 && monophyly_constraints[i].size() < num_taxa )
        {
        
            if ( monophyly_constraints[i].hasOptionalConstraints() == true )
            {
                std::vector<Clade> optional_constraints = monophyly_constraints[i].getOptionalConstraints();
                size_t idx = (size_t)( GLOBAL_RNG->uniform01() * optional_constraints.size() );
                sorted_clades.push_back( optional_constraints[idx] );
            }
            else
            {
                sorted_clades.push_back( monophyly_constraints[i] );
            }
        }
        
    }
    
    
    // create a clade that contains all species
    Clade all_species = Clade(taxa);
    sorted_clades.push_back(all_species);

    
    std::vector<Clade> virtual_taxa;
    int i = -1;
    for (std::vector<Clade>::iterator it = sorted_clades.begin(); it != sorted_clades.end(); it++)
    {
        // ignore negative clade constraints during simulation
        if ( it->isNegativeConstraint() == true )
        {
            continue;
        }
        
        ++i;
        const Clade &c = *it;
        std::vector<Taxon> taxa = c.getTaxa();
        std::vector<Clade> clades;
        
        int j = i;
        std::vector<Clade>::reverse_iterator jt(it);
        for (; jt != sorted_clades.rend(); jt++)
        {
            // ignore negative clade constraints during simulation
            if ( jt->isNegativeConstraint() == true )
            {
                continue;
            }
            
            j--;
            const Clade &c_nested = *jt;
            std::vector<Taxon> taxa_nested = c_nested.getTaxa();
            
            bool found_all = true;
            bool found_some = false;
            for (size_t k = 0; k < taxa_nested.size(); ++k)
            {
                std::vector<Taxon>::iterator kt = std::find(taxa.begin(), taxa.end(), taxa_nested[k]);
                if ( kt != taxa.end() )
                {
                    taxa.erase( kt );
                    found_some = true;
                }
                else
                {
                    found_all = false;
                }
                
            }
            
            if ( found_all == true )
            {
                clades.push_back( virtual_taxa[j] );
            }
            
            if ( found_all == false && found_some == true )
            {
                throw RbException("Cannot simulate tree: conflicting monophyletic clade constraints. Check that all clade constraints are properly nested.");
            }
            
        }
        
        
        std::vector<TopologyNode*> nodes_in_clade;
        
        
        for (size_t k = 0; k < taxa.size(); ++k)
        {
            Clade tmp_clade = Clade( taxa[k] );
            clades.push_back( tmp_clade );
        }
        
        for (size_t k = 0; k < clades.size(); ++k)
        {
            for (size_t j = 0; j < nodes.size(); ++j)
            {
                if (nodes[j]->getClade() == clades[k])
                {
                    nodes_in_clade.push_back( nodes[j] );
                    nodes.erase( nodes.begin()+j );
                    break;
                }
                
            }
            
        }
        
        tree_base_distribution->simulateClade(nodes_in_clade);
        nodes.push_back( nodes_in_clade[0] );
        
        std::vector<Taxon> v_taxa;
        nodes_in_clade[0]->getTaxa(v_taxa);
        Clade new_clade = Clade(v_taxa);
        virtual_taxa.push_back( new_clade );
    }
    
    TopologyNode *root = nodes[0];
    
    // initialize the topology by setting the root
    psi->setRoot(root, true);
    
    return psi;
}


/**
 * Set the DAG node.
 */
void TopologyConstrainedTreeDistribution::setStochasticNode( StochasticNode<Tree> *n )
{
    
    // delegate to base class first
    TypedDistribution<Tree>::setStochasticNode( n );
    
    if ( base_distribution != NULL )
    {
        base_distribution->setStochasticNode( n );
    }
    
}


/**
 * Set the current value.
 */
void TopologyConstrainedTreeDistribution::setValue(Tree *v, bool f )
{
    value->getTreeChangeEventHandler().removeListener( this );
    
    // we set our value to the same value as the base distribution
    // but first we need to make sure that our base class doesn't delete the value
    value = NULL;
    
    // and the we can set it for both ourselves and the base distribution
    TypedDistribution<Tree>::setValue(v, f);
    base_distribution->setValue(v, f);
    
    value->getTreeChangeEventHandler().addListener( this );
    
    initializeBitSets();
    
    // recompute the active clades
    dirty_nodes = std::vector<bool>( value->getNumberOfNodes(), true );
    
    recursivelyUpdateClades( value->getRoot() );
    
    stored_clades          = active_clades;
    stored_backbone_clades = active_backbone_clades;
}


/**
 * Swap the parameters held by this distribution.
 *
 *
 * \param[in]    oldP      Pointer to the old parameter.
 * \param[in]    newP      Pointer to the new parameter.
 */
void TopologyConstrainedTreeDistribution::swapParameterInternal( const DagNode *oldP, const DagNode *newP )
{
    
    if ( oldP == backbone_topologies )
    {
        backbone_topologies = static_cast<const TypedDagNode<RbVector<Tree> >* >( newP );
    }
    else if ( oldP == backbone_topology )
    {
        backbone_topology = static_cast<const TypedDagNode<Tree>* >( newP );
    }
    else
    {
        base_distribution->swapParameter(oldP,newP);
    }
    
}


/**
 * Touch the current value and reset some internal flags.
 * If the root age variable has been restored, then we need to change the root age of the tree too.
 */
void TopologyConstrainedTreeDistribution::touchSpecialization(const DagNode *affecter, bool touchAll)
{
    stored_clades = active_clades;
    stored_backbone_clades = active_backbone_clades;
    
    // if the root age wasn't the affecter, we'll set it in the base distribution here
    base_distribution->touch(affecter, touchAll);
}

void TopologyConstrainedTreeDistribution::keepSpecialization(const DagNode *affecter)
{
    stored_clades = active_clades;
    stored_backbone_clades = active_backbone_clades;
    
    base_distribution->keep(affecter);
}

void TopologyConstrainedTreeDistribution::restoreSpecialization(const DagNode *restorer)
{
    active_clades = stored_clades;
    active_backbone_clades = stored_backbone_clades;
    
    base_distribution->restore(restorer);
    
}
