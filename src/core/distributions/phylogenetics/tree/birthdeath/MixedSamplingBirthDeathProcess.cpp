////
////  MixedSamplingBirthDeathProcess.cpp
////  revbayes-proj
////
////  Created by Michael Landis on 11/16/17.
////  Copyright © 2017 Michael Landis. All rights reserved.
////
//
//#include "MixedSamplingBirthDeathProcess.h"
//#include "AbstractRootedTreeDistribution.h"
//#include "Clade.h"
//#include "RandomNumberFactory.h"
//#include "RandomNumberGenerator.h"
//#include "RbConstants.h"
//#include "RbException.h"
//#include "StochasticNode.h"
//#include "Taxon.h"
//#include "TopologyNode.h"
//#include "TreeUtilities.h"
//
//#include <algorithm>
//#include <cmath>
//
//using namespace RevBayesCore;
//
//
///**
// * Constructor.
// *
// * The constructor connects the parameters of the birth-death process (DAG structure)
// * and initializes the probability density by computing the combinatorial constant of the tree structure.
// *
// * \param[in]    c         Clade constraints.
// */
//MixedSamplingBirthDeathProcess::MixedSamplingBirthDeathProcess(TypedDistribution<Tree>* base_dist_1, TypedDistribution<Tree>* base_dist_2, const TypedDagNode<Tree>* bb) : TypedDistribution<Tree>( NULL ),
////    active_backbone_clades( base_dist->getValue().getNumberOfInteriorNodes(), RbBitSet() ),
//active_clades( base_dist_1->getValue().getNumberOfInteriorNodes() + base_dist_2->getValue().getNumberOfInteriorNodes() + 1, RbBitSet() ),
//backbone_topology( bb ),
//base_distribution_1( base_dist_1 ),
//base_distribution_2( base_dist_2 ),
//dirty_nodes( base_dist_1->getValue().getNumberOfNodes() + base_dist_2->getValue().getNumberOfNodes() + 1, true ),
//num_backbones( 0 )
//{
//    AbstractRootedTreeDistribution* tree_base_distribution_1 = dynamic_cast<AbstractRootedTreeDistribution*>(base_distribution_1);
//    if (tree_base_distribution_1 == NULL)
//    {
//        throw(RbException("Can only constrain tree distributions of type AbstractRootedTreeDistribution"));
//    }
//    AbstractRootedTreeDistribution* tree_base_distribution_2 = dynamic_cast<AbstractRootedTreeDistribution*>(base_distribution_2);
//    if (tree_base_distribution_2 == NULL)
//    {
//        throw(RbException("Can only constrain tree distributions of type AbstractRootedTreeDistribution"));
//    }
//    
//    
//    // add the parameters to our set (in the base class)
//    // in that way other class can easily access the set of our parameters
//    // this will also ensure that the parameters are not getting deleted before we do
//    
//    // add the parameters of the distribution
//    const std::vector<const DagNode*>& pars_1 = base_distribution_1->getParameters();
//    for (std::vector<const DagNode*>::const_iterator it = pars_1.begin(); it != pars_1.end(); ++it)
//    {
//        this->addParameter( *it );
//    }
//    const std::vector<const DagNode*>& pars_2 = base_distribution_2->getParameters();
//    for (std::vector<const DagNode*>::const_iterator it = pars_2.begin(); it != pars_2.end(); ++it)
//    {
//        this->addParameter( *it );
//    }
//    
//    addParameter( bb );
//    
//    // assign tree to value
////    value = &base_distribution_1->getValue();
//    
//    initializeBitSets();
//    redrawValue();
//}
//
//
///**
// * Copy Constructor.
// *
// * The constructor connects the parameters of the birth-death process (DAG structure)
// * and initializes the probability density by computing the combinatorial constant of the tree structure.
// *
// * \param[in]    c         Clade constraints.
// */
//MixedSamplingBirthDeathProcess::MixedSamplingBirthDeathProcess(const MixedSamplingBirthDeathProcess &d) : TypedDistribution<Tree>( d ),
//active_backbone_clades( d.active_backbone_clades ),
//active_clades( d.active_clades ),
//backbone_constraints( d.backbone_constraints ),
//backbone_mask( d.backbone_mask ),
//backbone_topology( d.backbone_topology ),
//base_distribution_1( d.base_distribution_1->clone() ),
//base_distribution_2( d.base_distribution_2->clone() ),
//dirty_nodes( d.dirty_nodes ),
//stored_backbone_clades( d.stored_backbone_clades ),
//stored_clades( d.stored_clades ),
//num_backbones( d.num_backbones )
//{
//    // the copy constructor of the TypedDistribution creates a new copy of the value
//    // however, here we want to hold exactly the same value as the base-distribution
//    // thus, we delete the newly created value
//    value->getTreeChangeEventHandler().removeListener( this );
//    delete value;
//    
//    // and then set it to the value of the base distribution
//    value = &base_distribution_1->getValue();
//    
//    value->getTreeChangeEventHandler().addListener( this );
//}
//
//
//
//MixedSamplingBirthDeathProcess::~MixedSamplingBirthDeathProcess()
//{
//    
//    delete base_distribution_1;
//    delete base_distribution_2;
//    
//    //value->getTreeChangeEventHandler().removeListener( this );
//    
//    // DO NOT DELETE THE VALUE
//    // the base distribution is the actual owner of the value!!!
//    // we simply avoid the deletion of the value by setting its pointer to NULL
//    // our base class, the TypedDistribution thinks that it owns the value and thus deletes it
//    value = NULL;
//    
//}
//
//
//MixedSamplingBirthDeathProcess* MixedSamplingBirthDeathProcess::clone( void ) const
//{
//    
//    return new MixedSamplingBirthDeathProcess( *this );
//}
//
//
///**
// * Compute the log-transformed probability of the current value under the current parameter values.
// *
// */
//double MixedSamplingBirthDeathProcess::computeLnProbability( void )
//{
//    recursivelyUpdateClades( value->getRoot() );
//    
//    // first check if the current tree matches the clade constraints
////    if ( !matchesConstraints() )
////    {
////        return RbConstants::Double::neginf;
////    }
//    
//    if ( !matchesBackbone() )
//    {
//        return RbConstants::Double::neginf;
//    }
//    
//    double lnProb = base_distribution_1->computeLnProbability() + base_distribution_2->computeLnProbability();
//    
//    return lnProb;
//}
//
//
//void MixedSamplingBirthDeathProcess::initializeBitSets(void)
//{
//    // reset the backbone constraints and mask
//    backbone_constraints.clear();
//    backbone_mask.clear();
//    backbone_constraints.resize(num_backbones);
//    backbone_mask.resize( num_backbones );
//    
//    // add the backbone constraints
//    backbone_mask[0] = RbBitSet( value->getNumberOfTips() );
//    backbone_mask[0] |= recursivelyAddBackboneConstraints( backbone_topology->getValue().getRoot(), 0 );
//}
//
//
//void MixedSamplingBirthDeathProcess::fireTreeChangeEvent(const TopologyNode &n, const unsigned& m)
//{
//    if (m == TreeChangeEventMessage::DEFAULT || m == TreeChangeEventMessage::TOPOLOGY)
//    {
//        
//        recursivelyFlagNodesDirty(n);
//    }
//}
//
//
///**
// * We check here if all the constraints are satisfied.
// * These are hard constraints, that is, the clades must be monophyletic.
// *
// * \return     True if the constraints are matched, false otherwise.
// */
//bool MixedSamplingBirthDeathProcess::matchesBackbone( void )
//{
//    
//    //    std::cout << base_distribution->getValue() << "\n";
//    // ensure that each backbone constraint is found in the corresponding active_backbone_clades
//    for (size_t i = 0; i < num_backbones; i++)
//    {
//        bool is_negative_constraint = false;
//        
//        std::vector<bool> negative_constraint_found( backbone_constraints[i].size(), false );
//        for (size_t j = 0; j < backbone_constraints[i].size(); j++)
//        {
//            std::vector<RbBitSet>::iterator it = std::find(active_backbone_clades[i].begin(), active_backbone_clades[i].end(), backbone_constraints[i][j] );
//            
//            // the search fails if the positive/negative backbone constraint is not satisfied
//            if (it == active_backbone_clades[i].end() && !is_negative_constraint )
//            {
//                // match fails if positive constraint is not found
//                return false;
//            }
//            else if (it != active_backbone_clades[i].end() && is_negative_constraint )
//            {
//                // match fails if negative constraint is found
//                negative_constraint_found[j] = true;
//            }
//        }
//        
//        // match fails if all negative backbone clades are found
//        bool negative_constraint_failure = true;
//        for (size_t j = 0; j < negative_constraint_found.size(); j++) {
//            if (negative_constraint_found[j] == false)
//                negative_constraint_failure = false;
//        }
//        if (negative_constraint_failure)
//        {
//            return false;
//        }
//    }
//    
//    // if no search has failed, then the match succeeds
//    return true;
//}
//
//
///**
// * We check here if all the monophyly constraints are satisfied.
// *
// * \return     True if the constraints are matched, false otherwise.
// */
////bool MixedSamplingBirthDeathProcess::matchesConstraints( void )
////{
////    for (size_t i = 0; i < monophyly_constraints.size(); i++)
////    {
////        
////        std::vector<Clade> constraints;
////        if (monophyly_constraints[i].isOptionalMatch())
////        {
////            constraints = monophyly_constraints[i].getOptionalConstraints();
////        }
////        else
////        {
////            constraints.push_back(monophyly_constraints[i]);
////        }
////        
////        std::vector<bool> constraint_satisfied( constraints.size(), false );
////        for (size_t j = 0; j < constraints.size(); j++) {
////            
////            std::vector<RbBitSet>::iterator it = std::find(active_clades.begin(), active_clades.end(), constraints[j].getBitRepresentation() );
////            
////            if (it != active_clades.end() && !constraints[j].isNegativeConstraint() )
////            {
////                constraint_satisfied[j] = true;
////            }
////            else if (it == active_clades.end() && constraints[j].isNegativeConstraint() )
////            {
////                constraint_satisfied[j] = true;
////            }
////        }
////        
////        // match fails if no optional positive or negative constraints satisfied
////        bool any_satisfied = false;
////        for (size_t j = 0; j < constraint_satisfied.size(); j++)
////        {
////            if (constraint_satisfied[j])
////            {
////                any_satisfied = true;
////                break;
////            }
////        }
////        if (!any_satisfied)
////            return false;
////        
////    }
////    
////    return true;
////}
//
//void MixedSamplingBirthDeathProcess::recursivelyFlagNodesDirty(const TopologyNode& n)
//{
//    
//    
//    dirty_nodes[ n.getIndex() ] = true;
//    
//    if ( n.isRoot() )
//        return;
//    
//    recursivelyFlagNodesDirty(n.getParent());
//    
//}
//
//
//RbBitSet MixedSamplingBirthDeathProcess::recursivelyAddBackboneConstraints( const TopologyNode& node, size_t backbone_idx )
//{
//    RbBitSet tmp( value->getNumberOfTips() );
//    
//    if ( node.isTip() )
//    {
//        const std::map<std::string, size_t>& taxon_map = value->getTaxonBitSetMap();
//        const std::string& name = node.getName();
//        std::map<std::string, size_t>::const_iterator it = taxon_map.find(name);
//        if (it == taxon_map.end()) {
//            
//            throw RbException("Taxon named " + it->first + " not found in tree's taxon map!");
//        }
//        tmp.set( it->second );
//    }
//    else
//    {
//        // get the child names
//        for (size_t i = 0; i < node.getNumberOfChildren(); i++)
//        {
//            tmp |= recursivelyAddBackboneConstraints( node.getChild(i), backbone_idx );
//        }
//        
//        if ( node.isRoot() == false )
//        {
//            backbone_constraints[backbone_idx].push_back(tmp);
//        }
//    }
//    
//    return tmp;
//}
//
//
//RbBitSet MixedSamplingBirthDeathProcess::recursivelyUpdateClades( const TopologyNode& node )
//{
//    if ( node.isTip() )
//    {
//        RbBitSet tmp = RbBitSet( value->getNumberOfTips() );
//        const std::map<std::string, size_t>& taxon_map = value->getTaxonBitSetMap();
//        const std::string& name = node.getName();
//        std::map<std::string, size_t>::const_iterator it = taxon_map.find(name);
//        tmp.set( it->second );
//        return tmp;
//    }
//    else if ( node.isRoot() )
//    {
//        if ( dirty_nodes[node.getIndex()] == true )
//        {
//            for (size_t i = 0; i < node.getNumberOfChildren(); i++)
//            {
//                recursivelyUpdateClades( node.getChild(i) );
//            }
//            
//            dirty_nodes[node.getIndex()] = false;
//        }
//        
//        return RbBitSet( value->getNumberOfTips(), true );
//    }
//    else
//    {
//        if ( dirty_nodes[node.getIndex()] == true )
//        {
//            RbBitSet tmp = RbBitSet( value->getNumberOfTips() );
//            for (size_t i = 0; i < node.getNumberOfChildren(); i++)
//            {
//                tmp |= recursivelyUpdateClades( node.getChild(i) );
//            }
//            
//            // update the clade
//            size_t idx = node.getIndex() - value->getNumberOfTips();
//            active_clades[idx] = tmp;
//            
//            for (size_t i = 0; i < num_backbones; i++) {
//                active_backbone_clades[i][idx] = tmp & backbone_mask[i];
//            }
//            
//            
//            dirty_nodes[node.getIndex()] = false;
//        }
//        
//        return active_clades[node.getIndex() - value->getNumberOfTips()];
//    }
//}
//
//
///**
// * Redraw the current value. We delegate this to the simulate method.
// */
//void MixedSamplingBirthDeathProcess::redrawValue( void )
//{
//    
//    Tree* new_value = simulateTree();
//    // base_distribution->redrawValue();
//    
//    value->getTreeChangeEventHandler().removeListener( this );
//    new_value->getTreeChangeEventHandler().addListener( this );
//    
//    // if we don't own the tree, then we just replace the current pointer with the pointer
//    // to the new value of the base distribution
//    value = new_value;
//    base_distribution->setValue( value );
//    
//    // recompute the active clades
//    dirty_nodes = std::vector<bool>( value->getNumberOfNodes(), true );
//    
//    recursivelyUpdateClades( value->getRoot() );
//    
//    stored_clades          = active_clades;
//    stored_backbone_clades = active_backbone_clades;
//}
//
//
//
///*
//void MixedSamplingBirthDeathProcess::setBackbone(const TypedDagNode<Tree> *backbone_one), const TypedDagNode<RbVector<Tree> > *backbone_many)
//{
//    if (backbone_one == NULL && backbone_many == NULL) {
//        ; // do nothing
//    } else if (backbone_one != NULL && backbone_many != NULL) {
//        ; // do nothing
//    } else {
//        
//        // clear old parameter
//        if (backbone_topology != NULL) {
//            this->removeParameter( backbone_topology );
//            backbone_topology = NULL;
//        } else {
//            this->removeParameter( backbone_topologies );
//            backbone_topologies = NULL;
//        }
//        
//        // set new parameter
//        if (backbone_one != NULL) {
//            backbone_topology = backbone_one;
//            num_backbones = 1;
//            use_multiple_backbones = false;
//            this->addParameter( backbone_one );
//        } else {
//            backbone_topologies = backbone_many;
//            num_backbones = backbone_topologies->getValue().size();
//            use_multiple_backbones = true;
//            this->addParameter( backbone_many );
//        }
//        
//        for (size_t i = 0; i < num_backbones; i++) {
//            std::vector<RbBitSet>v( base_distribution->getValue().getNumberOfInteriorNodes(), RbBitSet() );
//            active_backbone_clades.push_back(v);
//        }
//        backbone_mask = std::vector<RbBitSet>( num_backbones, base_distribution->getValue().getNumberOfInteriorNodes() );
//        
//        
//        initializeBitSets();
//        
//        // redraw the current value
//        if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
//        {
//            this->redrawValue();
//        }
//        
//    }
//}
//*/
// 
// 
///**
// *
// */
//Tree* MixedSamplingBirthDeathProcess::simulateTree( void )
//{
//    
//    // the time tree object (topology & times)
//    Tree *psi = new Tree();
//    
//    // internally we treat unrooted topologies the same as rooted
//    psi->setRooted( true );
//    
//    AbstractRootedTreeDistribution* tree_base_distribution_1 = dynamic_cast<AbstractRootedTreeDistribution*>( base_distribution_1 );
//    tree_base_distribution_1->redrawValue();
//    
//    AbstractRootedTreeDistribution* tree_base_distribution_2 = dynamic_cast<AbstractRootedTreeDistribution*>( base_distribution_2 );
//    tree_base_distribution_2->redrawValue();
//    
//    size_t num_taxa_1 = tree_base_distribution_1->getNumberOfTaxa();
//    const std::vector<Taxon> &taxa_1 = tree_base_distribution_1->getTaxa();
//    
//    size_t num_taxa_2 = tree_base_distribution_2->getNumberOfTaxa();
//    const std::vector<Taxon> &taxa_2 = tree_base_distribution_2->getTaxa();
//    
//    // create the tip nodes
//    std::vector<TopologyNode*> nodes;
//    for (size_t i=0; i < num_taxa_1; ++i)
//    {
//        // create the i-th taxon
//        TopologyNode* node = new TopologyNode( taxa_1[i], i );
//        
//        // set the age of this tip node
//        node->setAge( taxa_1[i].getAge() );
//        
//        // add the new node to the list
//        nodes.push_back( node );
//        
//    }
//    
//    
//    double ra = tree_base_distribution_1->getRootAge();
//    double max_age = tree_base_distribution_1->getOriginAge();
//    
//    
//    // placeholder to compile
//    TopologyNode* root = &tree_base_distribution_1->getValue().getRoot();
//    
//    // initialize the topology by setting the root
//    psi->setRoot(root, true);
//    
//    return psi;
//}
//
//
///**
// * Set the DAG node.
// */
//void MixedSamplingBirthDeathProcess::setStochasticNode( StochasticNode<Tree> *n )
//{
//    
//    // delegate to base class first
//    TypedDistribution<Tree>::setStochasticNode( n );
//    
//    if ( base_distribution != NULL )
//    {
//        base_distribution->setStochasticNode( n );
//    }
//    
//}
//
//
///**
// * Set the current value.
// */
//void MixedSamplingBirthDeathProcess::setValue(Tree *v, bool f )
//{
//    value->getTreeChangeEventHandler().removeListener( this );
//    
//    // we set our value to the same value as the base distribution
//    // but first we need to make sure that our base class doesn't delete the value
//    value = NULL;
//    
//    // and the we can set it for both ourselves and the base distribution
//    TypedDistribution<Tree>::setValue(v, f);
//    base_distribution->setValue(v, f);
//    
//    value->getTreeChangeEventHandler().addListener( this );
//    
//    initializeBitSets();
//    
//    // recompute the active clades
//    dirty_nodes = std::vector<bool>( value->getNumberOfNodes(), true );
//    
//    recursivelyUpdateClades( value->getRoot() );
//    
//    stored_clades          = active_clades;
//    stored_backbone_clades = active_backbone_clades;
//}
//
//
///**
// * Swap the parameters held by this distribution.
// *
// *
// * \param[in]    oldP      Pointer to the old parameter.
// * \param[in]    newP      Pointer to the new parameter.
// */
//void MixedSamplingBirthDeathProcess::swapParameterInternal( const DagNode *oldP, const DagNode *newP )
//{
//    
//    if ( oldP == backbone_topology )
//    {
//        backbone_topology = static_cast<const TypedDagNode<Tree>* >( newP );
//    }
//    else
//    {
//        base_distribution_1->swapParameter(oldP,newP);
//        base_distribution_2->swapParameter(oldP,newP);
//    }
//    
//}
//
///**
// * Touch the current value and reset some internal flags.
// * If the root age variable has been restored, then we need to change the root age of the tree too.
// */
//void MixedSamplingBirthDeathProcess::touchSpecialization(DagNode *affecter, bool touchAll)
//{
//    stored_clades = active_clades;
//    stored_backbone_clades = active_backbone_clades;
//    
//    // if the root age wasn't the affecter, we'll set it in the base distribution here
//    base_distribution->touch(affecter, touchAll);
//}
//
//void MixedSamplingBirthDeathProcess::keepSpecialization(DagNode *affecter)
//{
//    stored_clades = active_clades;
//    stored_backbone_clades = active_backbone_clades;
//    
//    base_distribution->keep(affecter);
//}
//
//void MixedSamplingBirthDeathProcess::restoreSpecialization(DagNode *restorer)
//{
//    active_clades = stored_clades;
//    active_backbone_clades = stored_backbone_clades;
//    
//    base_distribution->restore(restorer);
//    
//}
