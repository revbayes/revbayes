#include <stddef.h>
#include <algorithm>
#include <iosfwd>
#include <map>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "AbstractRootedTreeDistribution.h"
#include "Clade.h"
#include "TopologyConstrainedTreeDistribution.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
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
#include "TypedDagNode.h"
#include "TypedDistribution.h"

using std::vector;
using std::set;

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
TopologyConstrainedTreeDistribution::TopologyConstrainedTreeDistribution(TypedDistribution<Tree>* base_dist, const std::vector<Clade> &c) : TypedDistribution<Tree>( NULL ),
//    active_backbone_clades( base_dist->getValue().getNumberOfInteriorNodes(), RbBitSet() ),
    active_clades( base_dist->getValue().getNumberOfInteriorNodes(), RbBitSet() ),
    backbone_topology(NULL),
    backbone_topologies(NULL),
    base_distribution( base_dist ),
    dirty_nodes( base_dist->getValue().getNumberOfNodes(), true ),
    monophyly_constraints( c ),
    num_backbones( 0 ),
    use_multiple_backbones( false )

{
    AbstractRootedTreeDistribution* tree_base_distribution = dynamic_cast<AbstractRootedTreeDistribution*>(base_distribution);
    if (tree_base_distribution == NULL)
    {
        throw(RbException("Can only constrain tree distributions of type AbstractRootedTreeDistribution"));
    }
    
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
    
    initializeBitSets();
    redrawValue();
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
    use_multiple_backbones( d.use_multiple_backbones )
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
    
}



TopologyConstrainedTreeDistribution& TopologyConstrainedTreeDistribution::operator=(const TopologyConstrainedTreeDistribution &d)
{
    
    if ( this != &d )
    {
        TypedDistribution<Tree>::operator=( d );
        
        delete base_distribution;
        
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
        if (monophyly_constraints[i].isOptionalMatch() == false)
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
void TopologyConstrainedTreeDistribution::getAffected(RbOrderedSet<RevBayesCore::DagNode *> &affected, RevBayesCore::DagNode *affecter)
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
    
    //    std::cout << base_distribution->getValue() << "\n";
    // ensure that each backbone constraint is found in the corresponding active_backbone_clades
    for (size_t i = 0; i < num_backbones; i++)
    {
        bool is_negative_constraint = false;
        if (backbone_topology != NULL) {
            is_negative_constraint = backbone_topology->getValue().isNegativeConstraint();
        }
        else if (backbone_topologies != NULL) {
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
        for (size_t j = 0; j < negative_constraint_found.size(); j++) {
            if (negative_constraint_found[j] == false)
                negative_constraint_failure = false;
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
        if ( monophyly_constraints[i].isOptionalMatch() == true )
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

//    bool TopologyConstrainedTreeDistribution::matchesConstraints( void )
//    {
//        for (size_t i = 0; i < monophyly_constraints.size(); i++)
//        {
//            std::vector<RbBitSet>::iterator it = std::find(active_clades.begin(), active_clades.end(), monophyly_constraints[i].getBitRepresentation() );
//
//            if (it == active_clades.end() && !monophyly_constraints[i].isNegativeConstraint() )
//            {
//                // match fails if positive constraint is not found
//                return false;
//            }
//            else if (it != active_clades.end() && monophyly_constraints[i].isNegativeConstraint() )
//            {
//                // match fails if negative constraint is found
//                return false;
//            }
//        }
//
//        return true;
//    }




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
        
        return RbBitSet( value->getNumberOfTips(), true );
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
void TopologyConstrainedTreeDistribution::redrawValue( void )
{
    
    Tree* new_value = simulateTree();
    // base_distribution->redrawValue();
    
    value->getTreeChangeEventHandler().removeListener( this );
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




void TopologyConstrainedTreeDistribution::setBackbone(const TypedDagNode<Tree> *backbone_one, const TypedDagNode<RbVector<Tree> > *backbone_many)
{
    if (backbone_one == NULL && backbone_many == NULL) {
        ; // do nothing
    } else if (backbone_one != NULL && backbone_many != NULL) {
        ; // do nothing
    } else {
        
        
        // clear old parameter
        if (backbone_topology != NULL) {
            this->removeParameter( backbone_topology );
            backbone_topology = NULL;
        } else {
            this->removeParameter( backbone_topologies );
            backbone_topologies = NULL;
        }
        
        // set new parameter
        if (backbone_one != NULL) {
            backbone_topology = backbone_one;
            num_backbones = 1;
            use_multiple_backbones = false;
            this->addParameter( backbone_one );
        } else {
            backbone_topologies = backbone_many;
            num_backbones = backbone_topologies->getValue().size();
            use_multiple_backbones = true;
            this->addParameter( backbone_many );
        }
        
        for (size_t i = 0; i < num_backbones; i++) {
            std::vector<RbBitSet>v( base_distribution->getValue().getNumberOfInteriorNodes(), RbBitSet() );
            active_backbone_clades.push_back(v);
        }
        backbone_mask = std::vector<RbBitSet>( num_backbones, base_distribution->getValue().getNumberOfInteriorNodes() );
        
        
        initializeBitSets();
        
        // redraw the current value
        if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
        {
            this->redrawValue();
        }
        
    }
}

void set_ages_for_constraint_top(Clade& clade, const vector<Taxon>& taxa)
{
    // set the ages of each of the taxa in the constraint
    for (size_t j = 0; j < clade.size(); ++j)
    {
        bool found = false;
        for (auto& taxon: taxa)
        {
            if ( taxon.getName() == clade.getTaxonName(j) )
            {
                clade.setTaxonAge(j, taxon.getAge());
                found = true;
                break;
            }
        }
        if (not found)
            throw RbException("set_ages_for_constraint: can't find taxon " + clade.getTaxonName(j) + " in full taxon set!");
    }
}

void set_ages_for_constraint(Clade& clade, const vector<Taxon>& taxa)
{
    // set the ages of each of the taxa in the constraint
    set_ages_for_constraint_top( clade, taxa );

    // set ages for optional constraints
    std::vector<Clade> optional_constraints = clade.getOptionalConstraints();
    for (auto& optional_constraint: optional_constraints)
        set_ages_for_constraint_top( optional_constraint, taxa );

    clade.setOptionalConstraints( optional_constraints );
}

bool clade_nested_within(const Clade& clade1, const Clade& clade2)
{
    set<Taxon> taxa1;
    for(auto& taxon: clade1.getTaxa())
        taxa1.insert(taxon);

    set<Taxon> taxa2;
    for(auto& taxon: clade2.getTaxa())
        taxa2.insert(taxon);

    return std::includes(clade2.begin(), clade2.end(), clade1.begin(), clade1.end());
}

bool clades_overlap(const Clade& clade1, const Clade& clade2)
{
    set<Taxon> taxa1;
    for(auto& taxon: clade1.getTaxa())
        taxa1.insert(taxon);

    set<Taxon> taxa2;
    for(auto& taxon: clade2.getTaxa())
        taxa2.insert(taxon);

    auto i = taxa1.begin();
    auto j = taxa2.begin();
    while (i != taxa1.end() && j != taxa2.end())
    {
      if (*i == *j)
        return true;
      else if (*i < *j)
        ++i;
      else
        ++j;
    }
    return false;
}

bool clades_conflict(const Clade& clade1, const Clade& clade2)
{
    return clades_overlap(clade1, clade2) and
        (not clade_nested_within(clade1, clade2)) and
        (not clade_nested_within(clade2, clade2));
}

/**
 *
 */
Tree* TopologyConstrainedTreeDistribution::simulateTree( void )
{
    
    // the time tree object (topology & times)
    Tree *psi = new Tree();
    
    // internally we treat unrooted topologies the same as rooted
    psi->setRooted( true );
    
    AbstractRootedTreeDistribution* tree_base_distribution = dynamic_cast<AbstractRootedTreeDistribution*>( base_distribution );
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
    
    for (auto& monophyly_constraint: monophyly_constraints)
    {
        if ( monophyly_constraint.getAge() > max_age )
        {
            throw RbException("Cannot simulate tree: clade constraints are older than the origin age.");
        }

        // set the ages of each of the taxa in the constraint
        set_ages_for_constraint( monophyly_constraint, taxa );

        // populate sorted clades vector
        if ( monophyly_constraint.size() > 1 && monophyly_constraint.size() < num_taxa )
        {
            if ( monophyly_constraint.isOptionalMatch() == true )
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
    
//    for(std::vector<Clade>::iterator it = sorted_clades.begin(); it != sorted_clades.end(); it++)
//    {
//        std::cout << it->getAge() << std::endl;
//    }
    
    // DO WE NEED TO SORT THE TAXA?
    // try this crummy bubble sort
    size_t num_clades = sorted_clades.size();
    for (int i = 0; i < num_clades - 1; i++) {
        for(int j = 0; j < num_clades - i - 1; j++){
            if (sorted_clades[j].getAge() > sorted_clades[j+1].getAge()) {
                std::swap(sorted_clades[j], sorted_clades[j+1]);
            }
        }
    }
    
//    for(std::vector<Clade>::iterator it = sorted_clades.begin(); it != sorted_clades.end(); it++)
//    {
//        std::cout << it->getAge() << std::endl;
//    }

    
    std::vector<Clade> virtual_taxa;
    int i = -1;
    for (std::vector<Clade>::iterator it = sorted_clades.begin(); it != sorted_clades.end(); it++)
    {
        // ignore negative clade constraints during simulation
        if ( it->isNegativeConstraint() == true )
        {
            continue;
        }
        
//        std::cout << it->getAge() << std::endl;
        
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
                //                c.addTaxon( virtual_taxa[j] );
                //                taxa.push_back( virtual_taxa[j] );
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
            tmp_clade.setAge( taxa[k].getAge() );
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
        

        // here we need to start adding our "fake" tips
        for ( int index_missing_species = 0; index_missing_species < c.getNumberMissingTaxa(); ++index_missing_species)
        {
            ++n_added_missing_taxa;
            TopologyNode* missing_node = new TopologyNode("Missing_Taxon_" + StringUtilities::to_string(n_added_missing_taxa) );
            missing_node->setAge( 0.0 );
            nodes_in_clade.push_back( missing_node );
        }
        
        double clade_age = c.getAge();
        
        double max_node_age = 0;
        for (size_t j = 0; j < nodes_in_clade.size(); ++j)
        {
            if ( nodes_in_clade[j]->getAge() > max_node_age )
            {
                max_node_age = nodes_in_clade[j]->getAge();
            }
        }
        
        if ( clade_age <= max_node_age )
        {
            // Get the rng
//            RandomNumberGenerator* rng = GLOBAL_RNG;
            
//            clade_age = rng->uniform01() * ( max_age - max_node_age ) + max_node_age;
            clade_age = tree_base_distribution->simulateCladeAge(nodes_in_clade.size(), max_age, 0, max_node_age);
        }
        
        tree_base_distribution->simulateClade(nodes_in_clade, clade_age, 0.0);
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
void TopologyConstrainedTreeDistribution::touchSpecialization(DagNode *affecter, bool touchAll)
{
    stored_clades = active_clades;
    stored_backbone_clades = active_backbone_clades;
    
    // if the root age wasn't the affecter, we'll set it in the base distribution here
    base_distribution->touch(affecter, touchAll);
}

void TopologyConstrainedTreeDistribution::keepSpecialization(DagNode *affecter)
{
    stored_clades = active_clades;
    stored_backbone_clades = active_backbone_clades;
    
    base_distribution->keep(affecter);
}

void TopologyConstrainedTreeDistribution::restoreSpecialization(DagNode *restorer)
{
    active_clades = stored_clades;
    active_backbone_clades = stored_backbone_clades;
    
    base_distribution->restore(restorer);
    
}
