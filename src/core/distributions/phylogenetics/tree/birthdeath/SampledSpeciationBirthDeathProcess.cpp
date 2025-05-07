#include <algorithm>
#include <cmath>
#include <set>
#include <cstddef>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include "DistributionExponential.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "StochasticNode.h"
#include "TopologyNode.h"
#include "SampledSpeciationBirthDeathProcess.h"
#include "AbstractCharacterHistoryBirthDeathProcess.h"
#include "BranchHistory.h"
#include "BranchHistoryDiscrete.h"
#include "CharacterEvent.h"
#include "CharacterEventCompare.h"
#include "CharacterEventDiscrete.h"
#include "CharacterHistoryDiscrete.h"
#include "RbException.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "Taxon.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { template <class valueType> class RbOrderedSet; }

using namespace RevBayesCore;

SampledSpeciationBirthDeathProcess::SampledSpeciationBirthDeathProcess(const TypedDagNode<double> *a,
                                                                       const TypedDagNode<double> *s,
                                                                       const TypedDagNode<double> *e,
                                                                       const TypedDagNode< double > *r,
                                                                       const std::vector<Taxon> &n) :
AbstractCharacterHistoryBirthDeathProcess(),
root_age( a ),
speciation( s ),
extinction( e ),
rho( r ),
branch_histories( NULL, 1, 1, true ),
taxa( n ),
active_likelihood( std::vector<size_t>(2*n.size()-1, 0) ),
storedLikelihood( std::vector<std::vector<double> >(2*n.size()-1, std::vector<double>(2, 0.0))),
changed_nodes( std::vector<bool>(2*n.size()-1, false) ),
dirty_nodes( std::vector<bool>(2*n.size()-1, true) ),
scalingFactors( std::vector<std::vector<double> >(2*n.size()-1, std::vector<double>(2,0.0) ) ),
totalScaling( 0.0 )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( root_age );
    addParameter( speciation );
    addParameter( extinction );
    addParameter( rho );
    
    num_taxa = taxa.size();
    
    // the combinatorial factor for the probability of a labelled history is
    // 2^{n-1} / ( n! * (n-1)! )
    // but since the probability of the divergence times contains the factor (n-1)! we simply store
    // 2^{n-1} / n!
    double lnFact = 0.0;
    for (size_t i = 2; i <= num_taxa; i++)
    {
        lnFact += std::log(i);
    }
    
    logTreeTopologyProb = (num_taxa - 1) * RbConstants::LN2 - lnFact ;
    
    simulateTree();
    
}

SampledSpeciationBirthDeathProcess::~SampledSpeciationBirthDeathProcess()
{
    
}


void SampledSpeciationBirthDeathProcess::assertParentChildEdges(TopologyNode* n)
{
    if (!n->isRoot())
    {
        TopologyNode* p = &n->getParent();
        bool is_child = false;
        for (size_t i = 0; i < p->getNumberOfChildren(); i++)
        {
            if (n == &p->getChild(i))
                is_child = true;
        }
        
        if (!is_child)
        {
            
            std::cout << "is_child mismatch\n";
        }
        
    }
    
    if (!n->isTip())
    {
        bool is_parent0 = false;
        TopologyNode* ch0 = &n->getChild(0);
        if (n == &ch0->getParent())
            is_parent0 = true;
        
        bool is_parent1 = false;
        TopologyNode* ch1 = &n->getChild(1);
        if (n == &ch1->getParent())
            is_parent1= true;
        
        if (!is_parent0)
        {
            
            std::cout << "is_parent0 mismatch\n";
        }
        if (!is_parent1)
        {
            
            std::cout << "is_parent1 mismatch\n";
        }
        
    }

    
    ;
}

void SampledSpeciationBirthDeathProcess::assignNodes( TopologyNode* n, size_t& tip_index, size_t& int_index)
{
    
    if (n->isTip())
    {
        n->setTaxon(taxa[tip_index]);
//        n->setIndex(tip_index++);
        tip_index++;
    }
    else
    {
        for (size_t i = 0; i < n->getNumberOfChildren(); i++)
        {
            assignNodes( &n->getChild(i), tip_index, int_index );
        }
//        n->setIndex(int_index++);
    }
}


/* Clone function */
SampledSpeciationBirthDeathProcess* SampledSpeciationBirthDeathProcess::clone( void ) const
{
    SampledSpeciationBirthDeathProcess* v = new SampledSpeciationBirthDeathProcess( *this );
    v->computeLnProbability();
    return v;
}


/* Compute probability */
double SampledSpeciationBirthDeathProcess::computeLnProbability( void )
{
    // for now
    totalScaling = 0;
    
    // Variable declarations and initialization
    double lnProb = 0.0;
    double age = root_age->getValue();
    
    // we need to check that the root age matches
    if ( age != value->getRoot().getAge() )
    {
        return RbConstants::Double::neginf;
    }
    
    
    // check that the ages are in correct chronological order
    // i.e., no child is older than its parent
    const std::vector<TopologyNode*>& nodes = value->getNodes();
    for (std::vector<TopologyNode*>::const_iterator it = nodes.begin(); it != nodes.end(); it++)
    {
        
        const TopologyNode &the_node = *(*it);
        if ( the_node.isRoot() == false )
        {
            
            if ( (the_node.getAge() - (*it)->getParent().getAge()) > 0 && the_node.isSampledAncestorTip() == false )
            {
                return RbConstants::Double::neginf;
            }
            else if ( (the_node.getAge() - (*it)->getParent().getAge()) > 1E-6 && the_node.isSampledAncestorTip() == true )
            {
                return RbConstants::Double::neginf;
            }
            
        }
        
    }
    
    // check that the sampled ancestor nodes have a zero branch length
    for (std::vector<TopologyNode*>::const_iterator it = nodes.begin(); it != nodes.end(); it++)
    {
        const TopologyNode &the_node = *(*it);
        if ( the_node.isSampledAncestorTip() == true )
        {
            
            if ( the_node.isFossil() == false )
            {
                return RbConstants::Double::neginf;
            }
            else if ( the_node.getBranchLength() > 1E-6 )
            {
                return RbConstants::Double::neginf;
            }
            
        }
    }
    
    // add the survival of a second species if we condition on the MRCA
    lnProb += computeRootLikelihood();
    
    return lnProb; // + logTreeTopologyProb;
}

double SampledSpeciationBirthDeathProcess::computeLineageUnsampledByPresentProbability(double t_low, double t_sample)
{
    
    // First, get the probability that a lineage at time t_low leaves at least one sampled descendants at time t_sample
    double birthRate = speciation->getValue();
    double deathRate = extinction->getValue();
    double samplingProb = rho->getValue();

    // Prob( N(T)>0 | N(t)=1 )
    double p0 = (samplingProb * (birthRate - deathRate)) /
                (samplingProb * birthRate + ( (1.0 - samplingProb) * birthRate - deathRate) * exp(-(birthRate - deathRate) * (t_sample - t_low)));
    
    // Second, get the complement of the event (no sampled descendants)
    // Prob( N(T)=0 | N(t)=1) = 1 - Prob( N(T)>0 | N(t)=1 )
    double p1 = 1.0 - p0;
    
    return p1;
}

void SampledSpeciationBirthDeathProcess::computeNodeProbability(const RevBayesCore::TopologyNode &node, size_t node_index)
{
//    if (false && node.isRoot())
//    {
//        return;
//    }
//    elsee
    if (!node.isTip())
    {
        // this is an internal node
        const TopologyNode &left = node.getChild(0);
        size_t left_index = left.getIndex();
        computeNodeProbability( left, left_index );
        const TopologyNode &right = node.getChild(1);
        size_t right_index = right.getIndex();
        computeNodeProbability( right, right_index );
    }
    
    // check for recomputation
    bool bypassDirtyFlag = !false;
    
    if ( dirty_nodes[node_index] || bypassDirtyFlag )
    {
        
        double lnProb = 0.0;
        
        // get value
        const BranchHistory& bh = branch_histories[ node_index ];
        const std::multiset<CharacterEvent*,CharacterEventCompare>& hist = bh.getHistory();
        
        // process parameters
        const double &birth       = speciation->getValue();
        const double &death       = extinction->getValue();
        const double &sample_prob = rho->getValue();

        // branch/tree variables
        double branch_length = 0.0;
        double prev_age = 0.0;
        if (node.isRoot()) {
            branch_length = node.getAge();
            prev_age = 2 * node.getAge();
        }
        else {
            branch_length = node.getBranchLength();
            prev_age = node.getParent().getAge();
        }
        
        
        double end_age       = node.getAge();
        double sample_age    = 0.0; // NB: assumes the process ends at the present, T==0
        double prev_time     = 0.0;
        
        // compute probability for the observed and sampled speciation events on the branch
        for (std::multiset<CharacterEvent*,CharacterEventCompare>::const_reverse_iterator it=hist.rbegin(); it!=hist.rend(); ++it)
        {
            CharacterEvent* event = *it;
            double curr_time = event->getAge(); // CHECK THIS AGE
            double time_interval = curr_time - prev_time;
            double curr_age = prev_age - time_interval;

            
   
            // compute probability one lineage goes extinct by the present
            if (!node.isRoot()) {
                
                // compute the probability that the next event was a birth event
                //            std::cout << "\t\t\tlnProb\t" << lnProb << "\n";
                double v = log(birth) - (birth + death) * time_interval;
                lnProb += v ; // log(birth) - (birth + death) * time_interval;

                double p = computeLineageUnsampledByPresentProbability(-curr_age, sample_age);
                lnProb += log(p);
            }
            else
            {
                double v = log(birth) - (birth + death) * time_interval;
                lnProb += v ; // log(birth) - (birth + death) * time_interval;
            }
            
            // for survive,extinct and extinct,survive
            lnProb += 1.00 * log(2);
            
            // advance time
            prev_time = curr_time;
            prev_age  = curr_age;
        }
        
        double time_interval = prev_age - end_age;
        double v = 0.0;
        if ( node.isTip() ) {
            // if node is a tip, no further events occurred
            v = log(sample_prob) - (birth + death) * time_interval;
        }
        else {
            // if node is not a tip, the next event is a speciation event
            v = log(birth) - (birth + death) * time_interval;
        }
        lnProb += v;
//        std::cout << "\t" << (node.isTip() ? "C" : "A") << "-event\n";
//        std::cout << "\t\ttime_int\t" << time_interval << "\n";
//        std::cout << "\t\tprev_age\t" << prev_age << "\n";
//        std::cout << "\t\tcurr_age\t" << end_age << "\n";
//        std::cout << "\t\t\tlnProb\t" << lnProb << "\n";
//        std::cout << "---\n";
    
    
        // store likelihood
        storedLikelihood[node_index][ active_likelihood[node_index] ] = lnProb;

        // mark as computed
        dirty_nodes[node_index] = false;
    }
}

double SampledSpeciationBirthDeathProcess::computeRootLikelihood( void )
{
    
    const TopologyNode &root = value->getRoot();
    
    // fill the likelihoods
    if (!true) {
        ; //computeNodeProbability(root, root.getIndex() );
    }
    else {
        const TopologyNode &left = root.getChild(0);
        size_t left_index = left.getIndex();
        computeNodeProbability( left, left_index );
        const TopologyNode &right = root.getChild(1);
        size_t right_index = right.getIndex();
        computeNodeProbability( right, right_index );
    }
    
    // sum lnProbs across all nodes
    double lnProb = 0.0;
    for (size_t i = 0; i < storedLikelihood.size(); ++i)
    {
        lnProb += storedLikelihood[i][active_likelihood[i]];
    }
    
    return lnProb;
}


void SampledSpeciationBirthDeathProcess::executeMethod(const std::string &n, const std::vector<const DagNode *> &args, RbVector<long> &rv) const
{
    
    if ( n == "numberEvents" )
    {
        size_t num_branches = branch_histories.getNumberBranches();
        rv.clear();
        rv.resize( num_branches );
        
        for (size_t i = 0; i < num_branches; ++i)
        {
            rv[i] = int(branch_histories[i].getNumberEvents());
        }
        
    }
    else
    {
        throw RbException("The heterogeneous rate birth-death process does not have a member method called '" + n + "'.");
    }
    
}

void SampledSpeciationBirthDeathProcess::executeMethod(const std::string &n, const std::vector<const DagNode *> &args, RbVector<RbVector<double> > &rv) const
{
    
    if ( n == "eventTimes" )
    {
        size_t num_branches = branch_histories.getNumberBranches();
        rv.clear();
        rv.resize( num_branches );
        size_t num_slots = 20;
        for (size_t i = 0; i < num_branches; ++i)
        {
            rv[i].resize(num_slots, 0.0);
            
            size_t j = 0;
            std::multiset<CharacterEvent*, CharacterEventCompare> h = branch_histories[i].getHistory();
            for ( std::multiset<CharacterEvent*, CharacterEventCompare>::reverse_iterator it = h.rbegin(); it != h.rend(); it++ )
            {
                rv[i][j++] = (*it)->getAge(); // CHECK THIS AGE
                if (j > num_slots) break;
            }
        }
        
    }
    else
    {
        throw RbException("The heterogeneous rate birth-death process does not have a member method called '" + n + "'.");
    }
    
}

void SampledSpeciationBirthDeathProcess::getAffected(RbOrderedSet<DagNode *> &affected, const DagNode *affecter)
{
    
    if ( affecter == root_age)
    {
        dag_node->initiateGetAffectedNodes( affected );
    }
    
}


/**
 * Get the character history object.
 */
CharacterHistoryDiscrete& SampledSpeciationBirthDeathProcess::getCharacterHistory( void )
{
    
    return branch_histories;
}

/**
 * Get the character history object.
 */
const CharacterHistoryDiscrete& SampledSpeciationBirthDeathProcess::getCharacterHistory( void ) const
{
    
    return branch_histories;
}


void SampledSpeciationBirthDeathProcess::getLineagesAtAge(TopologyNode* n, std::vector<TopologyNode*>& nodes, double t)
{
    
    // add suitable branches for hidden speciation events to vector
    if (n->isRoot())
    {
        ; // ignore root (no subtending branch)
    }
    else if (n->getAge() < t && n->getParent().getAge() > t)
    {
        nodes.push_back(n);
    }
    
    // end recursion at tips
    if (n->isTip())
        return;
    
    // recurse
    for (size_t i = 0; i < n->getNumberOfChildren(); i++)
        getLineagesAtAge(&n->getChild(i), nodes, t);
    
}


/**
 * Keep the current value and reset some internal flags. Nothing to do here.
 */
void SampledSpeciationBirthDeathProcess::keepSpecialization(const DagNode *affecter)
{
    
    if ( affecter == root_age )
    {
        dag_node->keepAffected();
    }
    
}


void SampledSpeciationBirthDeathProcess::redrawValue( void )
{
    simulateTree();
}


/**
 * Restore the current value and reset some internal flags.
 * If the root age variable has been restored, then we need to change the root age of the tree too.
 */
void SampledSpeciationBirthDeathProcess::restoreSpecialization(const DagNode *affecter)
{
    
    if ( affecter == root_age )
    {
        value->getNode( value->getRoot().getIndex() ).setAge( root_age->getValue() );
        dag_node->restoreAffected();
    }
    
}



void SampledSpeciationBirthDeathProcess::setValue(Tree *v, bool force)
{
    
    // delegate to the parent class
    TypedDistribution< Tree >::setValue(v, force);
    
    branch_histories.setTree( value );
    
    // the true density is unknown, so this is an approximate simulated value
//    simulateEventsForTreeAdHoc();
    
    // TODO: consider running a mini-RJMCMC for fixed process parameters to initialize values
}

/** Simulate events for the given tree */
void SampledSpeciationBirthDeathProcess::simulateEventsForTreeAdHoc( void )
{
    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    // reset histories
    branch_histories.setTree( value );
    
    for (size_t i = 0; i < value->getNumberOfNodes(); i++)
    {
        TopologyNode& node = value->getNode(i);
        if (!node.isRoot())
        {
            double startAge = node.getParent().getAge();
            double endAge = node.getAge();
            double currAge = startAge;
            while (currAge > endAge)
            {
                double t = RbStatistics::Exponential::rv(speciation->getValue(), *rng);
                currAge -= t;
                
                if (currAge > endAge)
                {
                    double p = exp(computeLineageUnsampledByPresentProbability(-currAge, 0.0));
                    double u = rng->uniform01();
                    if (u < p)
                    {
                        double pos = (startAge - currAge) / (startAge - endAge);
                        CharacterEvent* evt = new CharacterEventDiscrete(0, 0, pos);
                        
                        branch_histories.addEvent(evt, i);
                    }
                }
            }
        }
    }
}

void SampledSpeciationBirthDeathProcess::simulateTree( void )
{
    double rootAge = root_age->getValue();
    bool failed = true;
    size_t max_tries = 10000;
    for (size_t sim_idx = 0; sim_idx < max_tries; sim_idx++)
    {
        // Create the time tree object (topology + times)
        Tree *psi = new Tree();
        
        // Root the topology by setting the appropriate flag
        psi->setRooted( true );
        
        // Create the root node and a vector of nodes
        TopologyNode* root = new TopologyNode();
        psi->setRoot(root, true);
        root->setAge(rootAge);
        std::set<TopologyNode* > nodes;
        std::vector<double> unsampledLineageAges;
        nodes.insert(root);
        
        // Simulate a tree
        simulateEvent(root, nodes, unsampledLineageAges, 0.0, rootAge);
        psi->setRoot(root, true);
        
//        std::cout << psi->getNumberOfTips() << "\n";
        // redo if the wrong number of taxa were sampled
        if (num_taxa != psi->getNumberOfTips() || root->getNumberOfChildren() != 2)
        {
            delete psi;
        }
        else
        {
            failed = false;
            
            
            // Initialize the topology by setting the root
            root->setParent(NULL);
            
            // get random taxon label indices
            std::vector<size_t> taxon_idx;
            for (size_t i = 0; i < taxa.size(); i++)
                taxon_idx.push_back(i);

            deprecated::random_shuffle(taxon_idx.begin(), taxon_idx.end());
            
            // Set names for terminal taxa
            for (size_t i = 0; i < taxa.size(); i++)
            {
                if (psi->getNode(i).isTip())
                    psi->getNode(i).setTaxon( taxa[ taxon_idx[i] ] );
                else
                    throw RbException("Tip nodes not assigned correct indices");
            }
            
            // Set the character histories
            branch_histories.setTree( psi );
            simulateUnsampledLineages(psi, unsampledLineageAges);
//            simulateUnsampledRootLineages( psi );
            
            // Error checking index problems
            for (size_t i = 0; i < psi->getNumberOfNodes(); i++)
                assertParentChildEdges( &psi->getNode(i) );

            // update value
            value = psi;
            
            // done!
            break;
        }
    }
    
    // If this is reached, the simulation failed
    if (failed) {
        simulateEventsForTreeAdHoc();
        throw RbException("The speciation sampled birth-death process failed to simulate a starting tree after 10000 tries.");
    }
}

size_t SampledSpeciationBirthDeathProcess::simulateEvent( TopologyNode* n, std::set<TopologyNode*>& nodes, std::vector<double>& unsampledLineageAges, double time, double maxTime)
{

    double age = maxTime - time;
    if (time > maxTime)
        age = 0.0;
    n->setAge(age);
    
    if (time > maxTime)
    {
        // tip nodes have no children
        n->removeAllChildren();
        
        // return tip count only if sampled
        if (GLOBAL_RNG->uniform01() < rho->getValue())
        {
//            std::cout << "s-tip node\t0.0\t" << n << "\n";
            return 1;
        }
        else
        {
//            std::cout << "u-tip node\t0.0\t" << n << "\n";
            return 0;
        }
    }
    else
    {
        // get BD process parameters
        double b = speciation->getValue();
        double d = extinction->getValue();
        
        // is the next event a birth?
        double u = GLOBAL_RNG->uniform01();
        bool is_birth = ( u < (b/(b+d)) || n->isRoot());
        
        // recurse
        size_t left_sampled  = 0;
        size_t right_sampled = 0;
        TopologyNode* left   = new TopologyNode();
        TopologyNode* right  = new TopologyNode();
        if (is_birth)
        {
            n->addChild(left);
            left->setParent(n);
            nodes.insert(left);
            double t_left = time + RbStatistics::Exponential::rv(b+d, *GLOBAL_RNG);
//            std::cout << "left\t" << t_left << "\t" << left << "\t" << n << "\n";
            left_sampled += simulateEvent(left, nodes, unsampledLineageAges, t_left, maxTime);
          
            n->addChild(right);
            right->setParent(n);
            nodes.insert(right);
            double t_right = time + RbStatistics::Exponential::rv(b+d, *GLOBAL_RNG);
//            std::cout << "right\t" << t_right << "\t" << right << "\t" << n << "\n";
            right_sampled += simulateEvent(right, nodes, unsampledLineageAges, t_right, maxTime);
        }
        else
        {
//            std::cout << "extinct\t" << n << "\n";
            return 0;
        }
        
        // recursion ends, do clean-up
//        std::cout << "CLEANUP\t" << n << "\t" << maxTime - n->getAge() << "\t" << left_sampled << "\t" << right_sampled << "\n";
//        std::cout << "PRECLEANUP\t" << n << "\t" << n->getNumberOfChildren() << "\n";
        if (left_sampled == 0 && right_sampled == 0)
        {
            // completely unsampled clade
            n->removeAllChildren();
            nodes.erase(left);
            nodes.erase(right);
        }
        else if ( n->isRoot() && (left_sampled == 0 || right_sampled == 0) )
        {
            // bad root state
            n->removeAllChildren();
            nodes.erase(left);
            nodes.erase(right);
        }
        else if (left_sampled == 0)
        {
            // remove left node
            n->removeChild(left);
            nodes.erase(left);
//            delete left;
            
            // convert node to unsampled birth event
            unsampledLineageAges.push_back(time);
            
            // patch parent to right node
            TopologyNode* p = &n->getParent();
            right = &n->getChild(0);
            n->setParent(NULL);
            n->removeChild(right); // right
            p->removeChild(n);
            p->addChild(right);
            right->setParent(p);

            // delete node
            nodes.erase(n);
//            delete n;

        }
        else if (right_sampled == 0)
        {
            // remove right node
            n->removeChild(right);
            nodes.erase(right);
//            delete right;
            
            // convert node to unsampled birth event
            unsampledLineageAges.push_back(time);
            
            // patch parent to left node
            TopologyNode* p = &n->getParent();
            left = &n->getChild(0);
            n->setParent(NULL);
            n->removeChild(left); // left
            p->removeChild(n);
            p->addChild(left);
            left->setParent(p);
            
            // delete node
            nodes.erase(n);
//            delete n;
        }
    
        return left_sampled + right_sampled;
    }
}

void SampledSpeciationBirthDeathProcess::simulateUnsampledLineages(Tree* t, std::vector<double> ages)
{
    TopologyNode* root = &t->getRoot();
    for (size_t i = 0; i < ages.size(); i++)
    {
        // find lineages alive at that age
        std::vector<TopologyNode*> nodes;
        getLineagesAtAge(root, nodes, ages[i]);
        
        // sample a random lineage
        size_t u = (size_t)(GLOBAL_RNG->uniform01() * nodes.size());
        size_t node_index = nodes[u]->getIndex();
        double time = nodes[u]->getAge() + nodes[u]->getBranchLength() - ages[i];
//        std::cout << nodes[u]->getIndex() << " " <<  nodes[u]->getBranchLength() << " " << nodes[u]->getAge() << " " << time << "\n";
//        std::cout << nodes[u]->getParent().getIndex() << " " << nodes[u]->getParent().getAge() <<  " -> " << nodes[u]->getIndex() << " " << nodes[u]->getAge() <<"\n";
//        std::cout << "\n";
        
        CharacterEvent* evt = new CharacterEventDiscrete(0, 0, time);
        branch_histories.addEvent(evt, node_index);
    }
    
    
    // simulate events for root lineage
    double rootBranchEndTime = root_age->getValue();
    double b = speciation->getValue();
    double d = extinction->getValue();
    
    bool keep = false;
    
    // root branch has no events if the death rate is 0
    if (d == 0.0) keep = true;
    std::vector<double> rootEventTimes;
    
    while (!keep) {
        bool stop = false;
        rootEventTimes.clear();
        double currentTime = 0.0;
        while (!stop) {
            double dt = RbStatistics::Exponential::rv(b+d, *GLOBAL_RNG);
            
            // reject if death event
            if (GLOBAL_RNG->uniform01() < (d / (b+d)))
            {
                keep = false;
                stop = true;
            }
            
            // reject if birth event with sampled descendants
            else if (GLOBAL_RNG->uniform01() < (1.0 - computeLineageUnsampledByPresentProbability(currentTime+dt, 0.0)))
            {
                keep = false;
                stop = true;
            }
            
            // accept and finish if branch length exceeded
            else if (currentTime+dt >= rootBranchEndTime)
            {
                keep = true;
                stop = true;
            }
            
            // accept and continue if birth event within branch length
            else if (currentTime+dt < rootBranchEndTime)
            {
                currentTime += dt;
                rootEventTimes.push_back(currentTime);
            }
        }
    };
     
    for (size_t i = 0; i < rootEventTimes.size(); i++)
    {
        CharacterEvent* evt = new CharacterEventDiscrete(0, 0, rootEventTimes[i]);
        branch_histories.addEvent(evt, root->getIndex());
    }
    
}



/** Swap a parameter of the distribution */
void SampledSpeciationBirthDeathProcess::swapParameterInternal( const DagNode *oldP, const DagNode *newP )
{
    if (oldP == root_age)
    {
        root_age = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == speciation)
    {
        speciation = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == extinction)
    {
        extinction = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == rho)
    {
        rho = static_cast<const TypedDagNode<double>* >( newP );
    }
    
}

/**
 * Touch the current value and reset some internal flags.
 * If the root age variable has been restored, then we need to change the root age of the tree too.
 */
void SampledSpeciationBirthDeathProcess::touchSpecialization(const DagNode *affecter, bool touchAll)
{
    
    if ( affecter == root_age )
    {
        value->getNode( value->getRoot().getIndex() ).setAge( root_age->getValue() );
        dag_node->touchAffected();
    }
    
}
