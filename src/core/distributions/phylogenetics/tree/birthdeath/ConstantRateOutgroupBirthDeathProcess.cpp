//
//  ConstantRateOutgroupBirthDeathProcess.cpp
//  revbayes-proj
//
//  Created by Michael Landis on 5/5/16.
//  Copyright © 2016 Michael Landis. All rights reserved.
//

#include "ConstantRateOutgroupBirthDeathProcess.h"

#include <cmath>
#include <algorithm>
#include <iostream>
#include <map>
#include <string>
#include <utility>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbMathCombinatorialFunctions.h"
#include "AbstractRootedTreeDistribution.h"
#include "TimeInterval.h"
#include "Tree.h"
#include "TreeChangeEventHandler.h"
#include "TreeChangeEventMessage.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }


using namespace RevBayesCore;


/**
 * Constructor.
 * We delegate most parameters to the base class and initialize the members.
 *
 * \param[in]    ra				Age or the root (=time of the process).
 * \param[in]    s              Speciation rate.
 * \param[in]    e              Extinction rate.
 * \param[in]    ri             Sampling probability at present time.
 * \param[in]    cdt            Condition of the process (none/survival/#Taxa).
 * \param[in]    tn             Taxa.
 */
ConstantRateOutgroupBirthDeathProcess::ConstantRateOutgroupBirthDeathProcess(    const TypedDagNode<double> *ra,
                                                                                 const TypedDagNode<double> *s,
                                                                                 const TypedDagNode<double> *e,
                                                                                 const TypedDagNode<double> *ri,
                                                                                 const TypedDagNode<double> *ro,
                                                                                 const std::string& ssi,
                                                                                 const std::string& sso,
                                                                                 const std::string& cdt,
                                                                                 const std::vector<Taxon> &tn,
                                                                                 const std::vector<Taxon> &tno) : AbstractBirthDeathProcess( ra, cdt, tn, false, NULL ),
lambda( s ),
mu( e ),
rho_ingroup( ri ),
rho_outgroup( ro ),
taxa_outgroup( tno ),
sampling_strategy_outgroup( sso ),
sampling_strategy_ingroup( ssi ),
active_clades( tn.size() - 1, RbBitSet() ),
dirty_nodes( 2 * tn.size() - 1, true )

{
    addParameter( lambda );
    addParameter( mu );
    addParameter( rho_ingroup );
    addParameter( rho_outgroup );
    
    // get taxon ingroup
    for (size_t i = 0; i < taxa.size(); i++) {
        std::vector<Taxon>::iterator it = std::find( taxa_outgroup.begin(), taxa_outgroup.end(), this->taxa[i] );
        if (it == taxa_outgroup.end()) {
            taxa_ingroup.push_back( this->taxa[i] );
        }
    }
    clade_ingroup = Clade( taxa_ingroup );
    
    simulateOutgroupTree();
    
    initializeBitSets();
    
    // need a second round of simulations that satisfy the clade constraints
    // ... for now assume the topology is provided ...
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
ConstantRateOutgroupBirthDeathProcess* ConstantRateOutgroupBirthDeathProcess::clone( void ) const
{
    
    return new ConstantRateOutgroupBirthDeathProcess( *this );
}



bool ConstantRateOutgroupBirthDeathProcess::isLnProbabilityNonZero( void )
{
    // first check if the current tree matches the clade constraints
    return true;
    //recursivelyUpdateClades( value->getRoot() );
    //return( matchesConstraints() );
    
}

/**
 * Compute the log-transformed probability of the current value under the current parameter values.
 *
 * \return    The log-probability density.
 */
double ConstantRateOutgroupBirthDeathProcess::computeLnProbabilityTimes( void ) const
{

    // proceed with probabilities
    
    double lnProbTimes = 0.0;

    // the process starts at the root age with one lineage
    TopologyNode* root = &value->getRoot();
    double process_time = root->getAge();
    
    // variable declarations and initialization
    double birth_rate = lambda->getValue();
    double death_rate = mu->getValue();
    
    double rho_in = 1.0;
    if ( sampling_strategy_ingroup == "uniform" ) {
        rho_in = rho_ingroup->getValue();
    }
    double rho_out = 1.0;
    if ( sampling_strategy_outgroup == "uniform" ) {
        rho_out = rho_outgroup->getValue();
    }
    
    // get node ages for ingroup
    TopologyNode* ingroup_mrca = findMrcaNode( clade_ingroup );
    std::vector<double> ingroup_ages = std::vector<double>();
    recursivelyCollectCladeAges(ingroup_mrca, ingroup_ages);
    std::sort(ingroup_ages.begin(), ingroup_ages.end());
    
    // remaining node ages are from the outgroup
    std::vector<double> outgroup_ages = std::vector<double>();
    recursivelyCollectCladeAges(root, outgroup_ages, ingroup_mrca);
    std::sort(outgroup_ages.begin(), outgroup_ages.end());
    
    // ingroup divergence time probabilities
    for (size_t i=0; i<ingroup_ages.size(); i++)
    {
        double t = ingroup_ages[i];
        lnProbTimes += log(birth_rate * pOne(t, birth_rate, death_rate, rho_in));
    }
    
    // outgroup divergence time probabilities
    lnProbTimes += log(pOne(process_time, birth_rate, death_rate, rho_out));
    for (size_t i=0; i<outgroup_ages.size(); i++)
    {
        double t = outgroup_ages[i];
        lnProbTimes += log(birth_rate * pOne(t, birth_rate, death_rate, rho_out));
    }
    
    // sampling strategies for ingroup/outgroup
    if (sampling_strategy_ingroup == "diversified") {
        lnProbTimes += computeLnProbabilityDiversifiedSampling( ingroup_ages, ingroup_mrca->getAge(), rho_in, taxa_ingroup.size() );
    }
    if (sampling_strategy_outgroup == "diversified") {
        lnProbTimes += computeLnProbabilityDiversifiedSampling( outgroup_ages, root->getAge(), rho_out, taxa_outgroup.size() );
    }

    // condition on one ingroup and one outgroup lineage surviving, under their respective sampling probs
    if (this->condition == "survival")
    {
        double ln_factor_to_remove = 2.0 * -log( pSurvival(0, root->getAge()) );
        double ln_factor_to_add = -log( pSurvival(0, root->getAge(), rho_out) * pSurvival(0, root->getAge(), rho_in ) );
        lnProbTimes += ln_factor_to_add - ln_factor_to_remove;
    }
    
    
    return lnProbTimes;
    
}

double ConstantRateOutgroupBirthDeathProcess::computeLnProbabilityDiversifiedSampling(std::vector<double> ages, double presentTime, double samplingProb, size_t num_taxa) const
{
    
    if (ages.size() < 1) {
        std::cout << "MJL: Warning! something weird w/ OutgroupBDP ages.size() for diversified sampling\n";
        return 0.0;
    }
    
    // We use equation (5) of Hoehna et al. "Inferring Speciation and Extinction Rates under Different Sampling Schemes"
    double last_event = ages[ages.size()-1];

    double p_0_T = 1.0 - pSurvival(0,presentTime,1.0) * exp( rateIntegral(0,presentTime) );
    double p_0_t = (1.0 - pSurvival(last_event,presentTime,1.0) * exp( rateIntegral(last_event,presentTime) ));
    double F_t = p_0_t / p_0_T;

    // get an estimate of the actual number of taxa
    double m = round(num_taxa / samplingProb);
    int m_int = m;
    if (m_int == 0) m_int = 1;
    double lnP = (m-num_taxa) * log(F_t) + log(RbMath::choose(m_int,num_taxa));
    
    return lnP;
}

void ConstantRateOutgroupBirthDeathProcess::initializeBitSets(void)
{
    // fill the monophyly constraints bitsets
    RbBitSet b( value->getNumberOfTips() );
    const std::map<std::string, size_t> &taxon_map = value->getTaxonBitSetMap();
    for (size_t j = 0; j < clade_ingroup.size(); j++)
    {
        const std::string &name = clade_ingroup.getTaxonName(j);
        std::map<std::string, size_t>::const_iterator it = taxon_map.find( name );
        if ( it == taxon_map.end() )
        {
            throw RbException("Could not find taxon with name '" + name + "'.");
        }
        size_t k = it->second;
        
        b.set(k);
    }
    clade_ingroup.setBitRepresentation( b );
    
}

TopologyNode* ConstantRateOutgroupBirthDeathProcess::findMrcaNode( Clade c ) const
{
    // get an arbitrary ingroup taxon
    TopologyNode* node = &value->getTipNodeWithName( c.getTaxonName(0) );
    
    // traverse tip-to-root path until the clade's bit representation is found
    bool found = false;
    while ( !node->isRoot() && !found )
    {
        // get active_clades index for node
        size_t idx = node->getIndex() - value->getNumberOfTips();
        
        if (active_clades[idx] == clade_ingroup.getBitRepresentation())
        {
            found = true;
        }
        
        if (!found)
        {
            node = &node->getParent();
        }
    };
    
    return node;
}

void ConstantRateOutgroupBirthDeathProcess::fireTreeChangeEvent(const TopologyNode &n, const unsigned& m)
{
    if (m == TreeChangeEventMessage::DEFAULT || m == TreeChangeEventMessage::TOPOLOGY)
    {
        
        recursivelyFlagNodesDirty(n);
    }
}


bool ConstantRateOutgroupBirthDeathProcess::matchesConstraints( void )
{
   
    bool constraint_satisfied = false;
    std::vector<RbBitSet>::iterator it = std::find(active_clades.begin(), active_clades.end(), clade_ingroup.getBitRepresentation() );
    
    if (it != active_clades.end())
    {
        constraint_satisfied = true;
    }
    
    return constraint_satisfied;
}

/**
 * Compute the probability of survival if the process starts with one species at time start and ends at time end.
 *
 *
 * \param[in]    start      Start time of the process.
 * \param[in]    end        End/stopping time of the process.
 *
 * \return Speciation rate at time t.
 */
double ConstantRateOutgroupBirthDeathProcess::pSurvival(double start, double end, double sampling_prob) const
{
    // variable declarations and initialization
    double birth_rate = lambda->getValue();
    double death_rate = mu->getValue();
//    double sampling_prob = rho_ingroup->getValue();
    
    double p0 = pZero(end, birth_rate, death_rate, sampling_prob);
    
    return 1.0 - p0;
}

double ConstantRateOutgroupBirthDeathProcess::pSurvival(double start, double end) const
{
    // variable declarations and initialization
    double birth_rate = lambda->getValue();
    double death_rate = mu->getValue();
    double sampling_prob = rho_outgroup->getValue();
    
    double p0 = pZero(end, birth_rate, death_rate, sampling_prob);
    
    return 1.0 - p0;
}


double ConstantRateOutgroupBirthDeathProcess::rateIntegral(double t_low, double t_high) const
{
    
    double rate = (lambda->getValue() - mu->getValue()) * (t_low - t_high);
    
    return rate;
}


void ConstantRateOutgroupBirthDeathProcess::recursivelyCollectCladeAges( TopologyNode* node,  std::vector<double>& ages, TopologyNode* barrier_node ) const
{
//    std::cout << node << " " << barrier_node << "\n";
    if (node == barrier_node) return;
    
    if (node->isTip()) return;
    
    ages.push_back( node->getAge() );
//    std::cout << node << " " << node->getAge() << " " << node->getName() << "\n";
    
    const std::vector<TopologyNode*> children = node->getChildren();
    for (size_t i = 0; i < children.size(); i++) {
        recursivelyCollectCladeAges( children[i], ages, barrier_node );
    }
//    std::cout << node->getAge() << " " << node->getName() << "\n";
    
}


/**
 * Simulate new speciation times.
 */
double ConstantRateOutgroupBirthDeathProcess::simulateDivergenceTime(double origin, double present) const
{
    
    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    // get the parameters
    double age = origin - present;
    double b = lambda->getValue();
    double d = mu->getValue();
    double r = rho_ingroup->getValue();
    
    // get a random draw
    double u = rng->uniform01();
    
    
    // compute the time for this draw
    // see Hartmann et al. 2010 and Stadler 2011
    double t = 0.0;
    if ( b > d )
    {
        t = ( log( ( (b-d) / (1 - (u)*(1-((b-d)*exp((d-b)*age))/(r*b+(b*(1-r)-d)*exp((d-b)*age) ) ) ) - (b*(1-r)-d) ) / (r * b) ) )  /  (b-d);
    }
    else
    {
        t = ( log( ( (b-d) / (1 - (u)*(1-(b-d)/(r*b*exp((b-d)*age)+(b*(1-r)-d) ) ) ) - (b*(1-r)-d) ) / (r * b) ) )  /  (b-d);
    }
    
    return present + t;
}



/**
 *
 */
void ConstantRateOutgroupBirthDeathProcess::simulateOutgroupTree( void )
{
    
    // the time tree object (topology & times)
    Tree *psi = new Tree();
    
    // internally we treat unrooted topologies the same as rooted
    psi->setRooted( true );
    
    // create the tip nodes
    std::vector<TopologyNode*> nodes;
    std::vector<TopologyNode*> nodes_ingroup;
    std::vector<TopologyNode*> nodes_outgroup;
    for (size_t i=0; i<taxa_ingroup.size(); ++i)
    {
        // create the i-th taxon
        TopologyNode* node = new TopologyNode( taxa_ingroup[i], i );
        
        // set the age of this tip node
        node->setAge( taxa_ingroup[i].getAge() );
        
        // add the new node to the list
        nodes_ingroup.push_back( node );
        nodes.push_back( node );
        
    }
    for (size_t i=0; i<taxa_outgroup.size(); ++i)
    {
        // create the i-th taxon
        TopologyNode* node = new TopologyNode( taxa_outgroup[i], i );
        
        // set the age of this tip node
        node->setAge( taxa_outgroup[i].getAge() );
        
        // add the new node to the list
        nodes_outgroup.push_back( node );
        nodes.push_back( node );
        
    }
    
    double ra = getRootAge();
    double max_age = getOriginAge();
    
    double max_node_age = 0;
    for (size_t j = 0; j < nodes.size(); ++j)
    {
        if ( nodes[j]->getAge() > max_node_age )
        {
            max_node_age = nodes[j]->getAge();
        }
    }
    if ( ra <= max_node_age )
    {
        if (ra > 0.0)
        {
            throw(RbException("Root age younger than oldest taxon age"));
        }
        
        // Get the rng
        RandomNumberGenerator* rng = GLOBAL_RNG;
        
        ra = rng->uniform01() * ( max_age - max_node_age ) + max_node_age;
    }
    
    simulateClade(nodes, ra, max_age, true);
    
    TopologyNode *root = nodes[0];
    
    // initialize the topology by setting the root
    psi->setRoot(root, true);
    
    // finally store the new value
    delete value;
    value = psi;
    
}


/**
 * Swap the parameters held by this distribution.
 *
 *
 * \param[in]    oldP      Pointer to the old parameter.
 * \param[in]    newP      Pointer to the new parameter.
 */
void ConstantRateOutgroupBirthDeathProcess::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == lambda)
    {
        lambda = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == mu)
    {
        mu = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == rho_ingroup)
    {
        rho_ingroup = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == rho_outgroup)
    {
        rho_outgroup = static_cast<const TypedDagNode<double>* >( newP );
    }
    else
    {
        // delegate the super-class
        AbstractBirthDeathProcess::swapParameterInternal(oldP, newP);
    }
    
}



double ConstantRateOutgroupBirthDeathProcess::pZero(double t, double b, double d, double r) const
{
    double v = 1.0 - r*(b-d)/(r*b+(b*(1-r)-d)*exp(-(b-d)*t));
    return v;
}


double ConstantRateOutgroupBirthDeathProcess::pOne(double t, double b, double d, double r) const
{
    double v = r*pow(b-d,2) * exp(-(b-d)*t) / pow(r*b+(b*(1-r)-d)*exp(-(b-d)*t),2);
    return v;
}



void ConstantRateOutgroupBirthDeathProcess::recursivelyFlagNodesDirty(const TopologyNode& n)
{
    dirty_nodes[ n.getIndex() ] = true;
    
    if ( n.isRoot() )
        return;
    
    recursivelyFlagNodesDirty(n.getParent());
    
}

RbBitSet ConstantRateOutgroupBirthDeathProcess::recursivelyUpdateClades( const TopologyNode& node )
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
            dirty_nodes[node.getIndex()] = false;
        }
        
        return active_clades[node.getIndex() - value->getNumberOfTips()];
    }
}


/**
 * Set the current value.
 */
void ConstantRateOutgroupBirthDeathProcess::setValue(Tree *v, bool f )
{
    value->getTreeChangeEventHandler().removeListener( this );
    
    // we set our value to the same value as the base distribution
    // but first we need to make sure that our base class doesn't delete the value
    value = v;
    
    // and the we can set it for both ourselves and the base distribution
//    TypedDistribution<Tree>::setValue(v, f);
//    base_distribution->setValue(v, f);
    
    value->getTreeChangeEventHandler().addListener( this );
    
    initializeBitSets();
    
    // recompute the active clades
    dirty_nodes = std::vector<bool>( value->getNumberOfNodes(), true );
    
    recursivelyUpdateClades( value->getRoot() );
    
    stored_clades          = active_clades;
}



void ConstantRateOutgroupBirthDeathProcess::touchSpecialization(const DagNode *affecter, bool touchAll)
{
    AbstractRootedTreeDistribution::touchSpecialization(affecter, touchAll);
    stored_clades = active_clades;
}

void ConstantRateOutgroupBirthDeathProcess::keepSpecialization(const DagNode *affecter)
{
    stored_clades = active_clades;
}

void ConstantRateOutgroupBirthDeathProcess::restoreSpecialization(const DagNode *restorer)
{
    active_clades = stored_clades;
}
