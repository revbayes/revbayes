#include <string>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "AbstractHomologousDiscreteCharacterData.h"
#include "GeneralizedLineageHeterogeneousBirthDeathSamplingProcess.h"
#include "RlAbstractHomologousDiscreteCharacterData.h"
#include "RlString.h"
#include "HomologousDiscreteCharacterData.h"
#include "DiscreteCharacterState.h"
#include "DiscreteTaxonData.h"
#include "NaturalNumbersState.h"
#include "TopologyNode.h"
#include "Taxon.h"
#include "Tree.h"
#include "RlUserInterface.h"


using namespace RevBayesCore;

/**
 * Constructor.
*/
GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::GeneralizedLineageHeterogeneousBirthDeathSamplingProcess(
	const std::vector<Taxon>&                                       taxa_,
	const TypedDagNode<double>*                                     age_,
	const std::string&                                              condition_type_,
	const TypedDagNode<Simplex >*                                   root_frequency_,
	const size_t                                                    num_states_,
	bool                                                            use_origin_,
	bool                                                            zero_indexed_,
	size_t                                                          n_proc_,
	double                                                          abs_tol_,
    double                                                          rel_tol_
) : TypedDistribution<Tree>( new TreeDiscreteCharacterData() ),
	n_proc(n_proc_),
	taxa(taxa_),
	age(age_),
	condition_type(condition_type_),
	num_states(num_states_),
	use_origin(use_origin_),
	zero_indexed(zero_indexed_),
	root_frequency(root_frequency_),
	dirty_nodes( std::vector<bool>(2 * taxa.size() - 1, true) )
{
	try {
		// create the pointer
		tp_ptr = Plugin::loader().createTensorPhyloLik();

		// set the random number generator
		RandomNumberGenerator* rng = GLOBAL_RNG;
		unsigned int seed = rng->getSeed();
		tp_ptr->setSeed( size_t(seed) );

	} catch (...) {
		throw RbException("TensorPhylo is not loaded (use loadPlugin(...)).");
	}
	
	// turn on/off debug
	tp_ptr->setConditionalProbCompatibilityMode(false);
	tp_ptr->setNumberOfThreads(n_proc);
	tp_ptr->setAbsoluteTolerance(abs_tol_);
	tp_ptr->setRelativeTolerance(rel_tol_);
	tp_ptr->setLikelihoodApproximator(TensorPhylo::Interface::approximatorVersion_t::SEQUENTIAL_BRANCHWISE);
//	tp_ptr->setLikelihoodApproximator(TensorPhylo::Interface::approximatorVersion_t::PARALLEL_BRANCHWISE);
//	tp_ptr->setDebugMode(TensorPhylo::Interface::debugMode_t::DBG_FILE, "test.out");

	// add the parameters
	addParameter(age);
	addParameter(root_frequency);

	// add this object to the tree change event handler
    value->getTreeChangeEventHandler().addListener( this );

    // simulate the tree
    simulateTree();

    // resize vectors for this number of nodes
    resizeVectors(value->getNumberOfNodes());

    // update the kernel
    updateRootFrequency(true);

    // turn on/off debug
	// tp_ptr->setApplyTreeLikCorrection(false);

    // set the condition type
    tp_ptr->setConditionalProbabilityType(TensorPhylo::Interface::TIME);
    if ( condition_type != "time" )
    {
        if ( use_origin )
        {
        	if ( condition_type == "survival" or condition_type == "sampledExtant" )
        	{
        		tp_ptr->setConditionalProbabilityType(TensorPhylo::Interface::STEM_SURVIVAL);
        	}
        	else if ( condition_type == "sampled" )
        	{
        		tp_ptr->setConditionalProbabilityType(TensorPhylo::Interface::STEM_ONE_SAMPLE);
        	}
        	else if ( condition_type == "tree" )
        	{
        		tp_ptr->setConditionalProbabilityType(TensorPhylo::Interface::STEM_TWO_SAMPLES);
        	}
        	else if ( condition_type == "treeExtant" )
        	{
        		tp_ptr->setConditionalProbabilityType(TensorPhylo::Interface::STEM_TWO_EXT_SAMPLES);
        	}
        	else
        	{
        		throw RbException("Invalid condition selected.");
        	}
        }
        else
        {
        	if ( condition_type == "survival" or condition_type == "sampledExtant" )
        	{
        		tp_ptr->setConditionalProbabilityType(TensorPhylo::Interface::ROOT_MRCA);
        	}
        	else if ( condition_type == "sampled" )
        	{
        		tp_ptr->setConditionalProbabilityType(TensorPhylo::Interface::ROOT_SAMPLING);
        	}
        	else if ( condition_type == "sampledMRCA" )
        	{
        		tp_ptr->setConditionalProbabilityType(TensorPhylo::Interface::ROOT_SAMPLING_AND_MRCA);
        	}
        	else
        	{
        		throw RbException("Invalid condition selected.");
        	}
        }
    }

}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
GeneralizedLineageHeterogeneousBirthDeathSamplingProcess* GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::clone( void ) const
{
	GeneralizedLineageHeterogeneousBirthDeathSamplingProcess* tmp = new GeneralizedLineageHeterogeneousBirthDeathSamplingProcess( *this );
    tmp->getValue().getTreeChangeEventHandler().addListener(tmp);
    return tmp;
}

/**
 * Destructor. Because we added ourselves as a reference to tau when we added a listener to its
 * TreeChangeEventHandler, we need to remove ourselves as a reference and possibly delete tau
 * when we die. All other parameters are handled by others.
 */
GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::~GeneralizedLineageHeterogeneousBirthDeathSamplingProcess( void )
{
    // We don't delete the params, because they might be used somewhere else too. The model needs to do that!

    // remove myself from the tree listeners
    value->getTreeChangeEventHandler().removeListener( this );
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::buildSerialSampledRandomBinaryTree(Tree *psi, std::vector<TopologyNode*> &nodes, const std::vector<double> &ages)
{

	// Get the rng
	RandomNumberGenerator* rng = GLOBAL_RNG;

	// make a vector of active nodes, and a vector of serial nodes
	std::vector<TopologyNode*> active_nodes;
	std::vector<TopologyNode*> extinct_nodes;

	for(int i = 0; i < nodes.size(); ++i)
	{
		if ( nodes.at(i)->getAge() == 0.0 )
		{
			active_nodes.push_back( nodes.at(i) );
		}
		else
		{
			extinct_nodes.push_back( nodes.at(i) );
		}
	}

	// loop backward through ages
	size_t num_ages = ages.size();
	double current_time = 0.0;

	for(int i = num_ages - 1; i >= 0; i--)
	{
		// get the current time
		double current_time = ages[i];

		// check if any extinct nods become active
		size_t num_extinct = extinct_nodes.size();
		for(int j = num_extinct - 1; j >= 0; --j)
		{
			if ( extinct_nodes.at(j)->getAge() < current_time )
			{
				// add the extinct node to the active nodes list, remove it from the extinct nodes list
				active_nodes.push_back( extinct_nodes.at(j) );
				extinct_nodes.erase( extinct_nodes.begin() + long(j) );
			}
		}

		// randomly draw one child (arbitrarily called left) node from the list of active nodes
		size_t left = static_cast<size_t>( floor( rng->uniform01() * active_nodes.size() ) );
		TopologyNode* leftChild = active_nodes.at(left);

		// remove the randomly drawn node from the list
		active_nodes.erase( active_nodes.begin() + long(left) );

		// randomly draw one child (arbitrarily called left) node from the list of active nodes
		size_t right = static_cast<size_t>( floor( rng->uniform01() * active_nodes.size() ) );
		TopologyNode* rightChild = active_nodes.at(right);

		// remove the randomly drawn node from the list
		active_nodes.erase( active_nodes.begin() + long(right) );

		// add the parent
		size_t num_taxa = taxa.size();
		TopologyNode* parent = new TopologyNode(i + num_taxa);
		parent->addChild(leftChild);
		parent->addChild(rightChild);
		leftChild->setParent(parent);
		rightChild->setParent(parent);
		parent->setAge( current_time );
		active_nodes.push_back(parent);
		nodes.push_back(parent);

	}

}

double GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::computeLnProbability(void)
{

	// make sure we need to recompute
	if ( probability_dirty == false )
	{
		return current_ln_prob;
	}
	else
	{
		// make sure all parameters are up-to-date
		// on the tensorphylo side
		prepareParameters(false);
	}

	// store the old likelihood
    old_ln_prob = current_ln_prob;

    // calculate a likelihood!
    // tp_ptr->writeStateToFile("params.dat");
	// tp_ptr->loadStateFromFile("state.dat");
    current_ln_prob = tp_ptr->computeLogLikelihood();

    // we are now able to reset tp approximator
    tp_can_reset = true;

    // flag the likelihood as up-to-date
    probability_dirty = false;

    // NOTE: this likelihood differs from the one computed in other revbayes
    // distributions because we include the probability density of the root node.
    // for lineage-homogeneous (state-independent) models, that probability density
    // drops out if we condition on survival, in which case this likelihood
    // will be equivalent to a likelihood from other revbayes birth-death models

    return current_ln_prob;
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::dumpModel(std::string file_name)
{
	tp_ptr->writeStateToFile(file_name);
}

/**
 * Resize various vectors depending on the current number of nodes.
 */
void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::resizeVectors(size_t num_nodes)
{
    dirty_nodes = std::vector<bool>(num_nodes, true);
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::fireTreeChangeEvent(const TopologyNode &n, const unsigned& m)
{

	// mark tree as dirty
	tree_dirty = true;

	// mark nodes dirty
    recursivelyFlagNodeDirty( n );

	// update the newick string to tensorphylo
//	updateTree(); // TODO: move me somewhere later

}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::recursivelyFlagNodeDirty( const RevBayesCore::TopologyNode &n ) {

    // we need to flag this node and all ancestral nodes for recomputation
    size_t index = n.getIndex();

    // if this node is already dirty, the also all the ancestral nodes must have been flagged as dirty
    if ( dirty_nodes[index] == false )
    {

    	// the root doesn't have an ancestor
        if ( n.isRoot() == false )
        {
            recursivelyFlagNodeDirty( n.getParent() );
        }

        // set the flag
        dirty_nodes[index] = true;

    }

}


const RevBayesCore::AbstractHomologousDiscreteCharacterData& GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::getCharacterData() const
{
    return static_cast<TreeDiscreteCharacterData*>(this->value)->getCharacterData();
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::redrawValue(void)
{
	// just simulate with a uniform tree and all missing data
    simulateTree();
    
    // simulate tree under SSE model (implementation in progress)
    // simulateSSETree();
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::setLambda(const TypedDagNode< RbVector<double> >* param)
{
	if ( lambda_const != NULL || lambda_var != NULL )
	{
		throw RbException("Tried to set lambda twice.");
	}

	// set the value
	lambda_const = param;
	lambda_times = NULL;

	// include the parameter
	addParameter(lambda_const);

	// flag for update
	lambda_dirty = true;
	probability_dirty = true;
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::setLambda(const TypedDagNode< RbVector< RbVector<double> > >* param, const TypedDagNode< RbVector<double> >* times)
{
	if ( lambda_const != NULL || lambda_var != NULL )
	{
		throw RbException("Tried to set lambda twice.");
	}

	// check for ascending times
	checkTimesAreAscending(times->getValue());

	// set the value
	lambda_var   = param;
	lambda_times = times;

	// include the parameter
	addParameter(lambda_var);
	addParameter(lambda_times);

	// flag for update
	lambda_dirty = true;
	probability_dirty = true;

}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::setMu(const TypedDagNode< RbVector<double> >* param)
{
	if ( mu_const != NULL || mu_var != NULL )
	{
		throw RbException("Tried to set mu twice.");
	}

	// set the value
	mu_const = param;
	mu_times = NULL;

	// include the parameter
	addParameter(mu_const);

	// flag for update
	mu_dirty = true;
	probability_dirty = true;
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::setMu(const TypedDagNode< RbVector< RbVector<double> > >* param, const TypedDagNode< RbVector<double> >* times)
{
	if ( mu_const != NULL || mu_var != NULL )
	{
		throw RbException("Tried to set mu twice.");
	}

	// check for ascending times
	checkTimesAreAscending(times->getValue());

	// set the value
	mu_var   = param;
	mu_times = times;

	// include the parameter
	addParameter(mu_var);
	addParameter(mu_times);

	// flag for update
	mu_dirty = true;
	probability_dirty = true;
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::setPhi(const TypedDagNode< RbVector<double> >* param)
{
	if ( phi_const != NULL || phi_var != NULL )
	{
		throw RbException("Tried to set phi twice.");
	}

	// set the value
	phi_const = param;
	phi_times = NULL;

	// include the parameter
	addParameter(phi_const);

	// flag for update
	phi_dirty = true;
	probability_dirty = true;
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::setPhi(const TypedDagNode< RbVector< RbVector<double> > >* param, const TypedDagNode< RbVector<double> >* times)
{
	if ( phi_const != NULL || phi_var != NULL )
	{
		throw RbException("Tried to set phi twice.");
	}

	// check for ascending times
	checkTimesAreAscending(times->getValue());

	// set the value
	phi_var   = param;
	phi_times = times;

	// include the parameter
	addParameter(phi_var);
	addParameter(phi_times);

	// flag for update
	phi_dirty = true;
	probability_dirty = true;
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::setDelta(const TypedDagNode< RbVector<double> >* param)
{
	if ( delta_const != NULL || delta_var != NULL )
	{
		throw RbException("Tried to set delta twice.");
	}

	// set the value
	delta_const = param;
	delta_times = NULL;

	// include the parameter
	addParameter(delta_const);

	// flag for update
	delta_dirty = true;
	probability_dirty = true;
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::setDelta(const TypedDagNode< RbVector< RbVector<double> > >* param, const TypedDagNode< RbVector<double> >* times)
{
	if ( delta_const != NULL || delta_var != NULL )
	{
		throw RbException("Tried to set delta twice.");
	}

	// check for ascending times
	checkTimesAreAscending(times->getValue());

	// set the value
	delta_var   = param;
	delta_times = times;

	// include the parameter
	addParameter(delta_var);
	addParameter(delta_times);

	// flag for update
	delta_dirty = true;
	probability_dirty = true;
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::setUpsilon(const TypedDagNode< RbVector< RbVector<double> > >* param, const TypedDagNode< RbVector<double> >* times)
{
	if ( upsilon != NULL )
	{
		throw RbException("Tried to set upsilon twice.");
	}

	// set the value
	upsilon       = param;
	upsilon_times = times;

	// include the parameter
	addParameter(upsilon);
	addParameter(upsilon_times);

	// flag for update
	upsilon_dirty = true;
	probability_dirty = true;
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::setGamma(const TypedDagNode< RbVector< RbVector<double> > >* param, const TypedDagNode< RbVector<double> >* times)
{
	if ( gamma != NULL )
	{
		throw RbException("Tried to set gamma twice.");
	}

	// set the value
	gamma       = param;
	gamma_times = times;

	// include the parameter
	addParameter(gamma);
	addParameter(gamma_times);

	// flag for update
	gamma_dirty = true;
	probability_dirty = true;
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::setRho(const TypedDagNode< double >* param)
{
	if ( rho_simple != NULL || rho != NULL )
	{
		throw RbException("Tried to set rho twice.");
	}

	// set the value
	rho_simple = param;
	rho_times  = NULL;

	// include the parameter
	addParameter(rho_simple);

	// flag for update
	rho_dirty = true;
	probability_dirty = true;
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::setRho(const TypedDagNode< RbVector< RbVector<double> > >* param, const TypedDagNode< RbVector<double> >* times)
{
	if ( rho_simple != NULL || rho != NULL )
	{
		throw RbException("Tried to set rho twice.");
	}

	// set the value
	rho       = param;
	rho_times = times;

	// include the parameter
	addParameter(rho);
	addParameter(rho_times);

	// flag for update
	rho_dirty = true;
	probability_dirty = true;
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::setXi(const TypedDagNode< RbVector< RbVector<double> > >* param, const TypedDagNode< RbVector<double> >* times)
{
	if ( xi != NULL )
	{
		throw RbException("Tried to set xi twice.");
	}

	// set the value
	xi       = param;
	xi_times = times;

	// include the parameter
	addParameter(xi);
	addParameter(xi_times);

	// flag for update
	xi_dirty = true;
	probability_dirty = true;
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::setEta(const TypedDagNode< double >* param)
{
	if ( eta_simple != NULL || eta_const != NULL || eta_var != NULL )
	{
		throw RbException("Tried to set eta twice.");
	}

	// set the value
	eta_simple = param;
	eta_times  = NULL;

	// include the parameter
	addParameter(eta_simple);

	// flag for update
	eta_dirty = true;
	probability_dirty = true;
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::setEta(const TypedDagNode< RateGenerator >* param)
{
	if ( eta_simple != NULL || eta_const != NULL || eta_var != NULL )
	{
		throw RbException("Tried to set eta twice.");
	}

	// set the value
	eta_const = param;
	eta_times = NULL;

	// include the parameter
	addParameter(eta_const);

	// flag for update
	eta_dirty = true;
	probability_dirty = true;
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::setEta(const TypedDagNode< RbVector< RateGenerator > >* param, const TypedDagNode< RbVector<double> >* times)
{
	if ( eta_simple != NULL ||  eta_const != NULL || eta_var != NULL )
	{
		throw RbException("Tried to set eta twice.");
	}

	// check for ascending times
	checkTimesAreAscending(times->getValue());

	// set the value
	eta_var   = param;
	eta_times = times;

	// include the parameter
	addParameter(eta_var);
	addParameter(eta_times);

	// flag for update
	eta_dirty = true;
	probability_dirty = true;
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::setOmega(const TypedDagNode< CladogeneticProbabilityMatrix >* param)
{
	if ( omega_const != NULL || omega_var != NULL )
	{
		throw RbException("Tried to set omega twice.");
	}

	// set the value
	omega_const = param;
	omega_times = NULL;

	// include the parameter
	addParameter(omega_const);

	// flag for update
	omega_dirty = true;
	probability_dirty = true;
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::setOmega(const TypedDagNode< RbVector< CladogeneticProbabilityMatrix > >* param, const TypedDagNode< RbVector<double> >* times)
{
	if ( omega_const != NULL || omega_var != NULL )
	{
		throw RbException("Tried to set omega twice.");
	}

	// check for ascending times
	checkTimesAreAscending(times->getValue());

	// set the value
	omega_var   = param;
	omega_times = times;

	// include the parameter
	addParameter(omega_var);
	addParameter(omega_times);

	// flag for update
	omega_dirty = true;
	probability_dirty = true;
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::setZeta(const TypedDagNode< RbVector< MatrixReal > >* param)
{
	if ( zeta != NULL )
	{
		throw RbException("Tried to set zeta twice.");
	}

	// set the value
	zeta = param;

	// include the parameter
	addParameter(zeta);

	// flag for update
	zeta_dirty = true;
	probability_dirty = true;
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::setValue(Tree *v, bool f)
{

    // Suppress outdegree-1 internal nodes (= sampled ancestors)
    v->suppressOutdegreeOneNodes(true);

	// Check that tree is binary. This may still not be the case if there are multifurcations.
    if (v->isBinary() == false)
    {
        throw RbException("The process is only implemented for binary trees.");
    }

    // check the taxon labels in the current tree and new tree
    std::vector<std::string> param_taxa;
    size_t num_taxa = taxa.size();
    for(size_t i = 0; i < num_taxa; ++i)
    {
    	param_taxa.push_back( taxa[i].getName() );
    }

//    // get the taxa names from the input tree
//    std::vector<Taxon> input_taxa_obj = v->getTaxa();
//    std::vector<std::string> input_taxa;
//    for (size_t i = 0; i < input_taxa_obj.size(); ++i)
//	{
//    	input_taxa.push_back( input_taxa_obj[i].getSpeciesName() );
//    }
//
//    // get the taxa names from the current tree
//    std::vector<Taxon> current_taxa_obj = value->getTaxa();
//    std::vector<std::string> current_taxa;
//    for (size_t i = 0; i < current_taxa_obj.size(); ++i)
//	{
//    	current_taxa.push_back( current_taxa_obj[i].getSpeciesName() );
//    }

    std::vector<std::string> input_taxa   = v->getSpeciesNames();
    std::vector<std::string> current_taxa = value->getSpeciesNames();

    // check that the number of taxa match
    if ( input_taxa.size() != num_taxa )
    {
    	throw RbException("Number of taxa in the new tree does not match number of taxa in the old tree.");
    }

    // sort the labels
    std::sort(param_taxa.begin(),   param_taxa.end());
    std::sort(input_taxa.begin(),   input_taxa.end());
    std::sort(current_taxa.begin(), current_taxa.end());
    for(size_t i = 0; i < num_taxa; ++i)
    {
    	if ( current_taxa[i] != param_taxa[i] )
    	{
    		throw RbException("Something went wrong. The taxon labels in the current tree do not match the taxon labels in the distribution object.");
    	}
    	if ( input_taxa[i] != current_taxa[i] )
    	{
    		throw RbException("Taxon labels in the new tree do not match with taxon labels in the original tree.");
    	}
    }

    // swap the listener
    value->getTreeChangeEventHandler().removeListener( this );

    // set the tree
    static_cast<TreeDiscreteCharacterData *>(this->value)->setTree( *v );

    // clear memory
    delete v;

    // set the taxon data
    for(size_t i = 0; i < num_taxa; ++i)
    {
    	// get the taxon data and node in the new tree
    	Taxon &taxon = taxa[i];
    	TopologyNode& node = this->value->getTipNodeWithName( taxon.getSpeciesName() );

    	// set the taxon data for the node
    	node.setTaxon( taxon );
    }

    // add this object to the tree change event handler
    value->getTreeChangeEventHandler().addListener( this );

    // update the kernel
    updateTree(true);

    // mark the likelihood dirty
    probability_dirty = true;

}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::drawStochasticCharacterMap(std::vector<std::string>& character_histories, bool use_simmap_default)
{
	// draw the stochastic map
	TensorPhylo::Interface::mapHistories_t history = tp_ptr->drawHistory();

	// translate the map to a vector of strings
	size_t node_index = 0;
	for(TensorPhylo::Interface::mapHistories_t::iterator it = history.begin(); it != history.end(); ++it)
	{
		// get the history for the branch
		std::vector< std::pair<double, size_t> > this_history = it->second;
        size_t n_events = this_history.size();

		// create the string
        std::string simmap_string = "{";
        if (use_simmap_default == true)
        {
            std::vector< std::pair<double, size_t> >::reverse_iterator first_event = this_history.rend();
            first_event--; // apparently, this isn't for the reverse_iterator; uncomment to verify for yourself
            int event_idx = 0;
            for(std::vector< std::pair<double, size_t> >::reverse_iterator jt = this_history.rbegin(); jt != this_history.rend(); ++jt)
            {
                simmap_string = simmap_string + StringUtilities::toString(jt->second) + "," + StringUtilities::toString(jt->first);
                if ( jt != first_event && event_idx != (n_events - 1))
                {
                    simmap_string = simmap_string + ":";
                }
                event_idx++;
            }
        }
        else
        {
            std::vector< std::pair<double, size_t> >::iterator first_event = this_history.begin();
            for(std::vector< std::pair<double, size_t> >::iterator jt = this_history.begin(); jt != this_history.end(); ++jt)
            {
                if ( jt != first_event)
                {
                    simmap_string = simmap_string + ":";
                }
                simmap_string = simmap_string + StringUtilities::toString(jt->second) + "," + StringUtilities::toString(jt->first);
            }
        }
        simmap_string = simmap_string + "}";

        // add the string to the character histories
        character_histories[it->first - 1] = simmap_string;
	}
    
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::drawStochasticCharacterMap(std::vector<std::string>& character_histories, std::vector<double>& branch_lambda, std::vector<double>& branch_mu, std::vector<double>& branch_phi, std::vector<double>& branch_delta, std::vector<long>& num_events)
{
	// draw the stochastic map
	TensorPhylo::Interface::mapHistories_t history = tp_ptr->drawHistoryAndComputeRates(branch_lambda, branch_mu, branch_phi, branch_delta, num_events);

	// translate the map to a vector of strings
	for(TensorPhylo::Interface::mapHistories_t::iterator it = history.begin(); it != history.end(); ++it)
	{
		// get the history for the branch
		std::vector< std::pair<double, size_t> > this_history = it->second;

		// create the string
        std::string simmap_string = "{";
        std::vector< std::pair<double, size_t> >::reverse_iterator first_event = this_history.rend();
        //first_event--; // apparently, this isn't for the reverse_iterator; uncomment to verify for yourself
        for(std::vector< std::pair<double, size_t> >::reverse_iterator jt = this_history.rbegin(); jt != this_history.rend(); ++jt)
        {
            simmap_string = simmap_string + StringUtilities::toString(jt->second) + "," + StringUtilities::toString(jt->first);
            if ( jt != first_event )
            {
                simmap_string = simmap_string + ":";
            }
        }
        simmap_string = simmap_string + "}";

        // add the string to the character histories
        character_histories[it->first - 1] = simmap_string;

	}

}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::drawJointConditionalAncestralStates(std::vector<size_t>& startStates, std::vector<size_t>& endStates)
{
	// draw the ancestral states
	TensorPhylo::Interface::mapHistories_t history = tp_ptr->drawAncestralStates();

	// translate the ancestral states to revbayes format
	for(TensorPhylo::Interface::mapHistories_t::iterator it = history.begin(); it != history.end(); ++it)
	{
		// get the history for the branch
		std::vector< std::pair<double, size_t> > this_history = it->second;

        // add the states
		startStates[it->first - 1] = this_history[0].second;
		endStates[it->first - 1]   = this_history[1].second;
	}

}


void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::getAffected(RbOrderedSet<DagNode *>& affected, DagNode* affecter)
{
    if ( affecter == age )
    {
        dag_node->initiateGetAffectedNodes( affected );
    }
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::initializeEmptyCharData(void)
{
	// get the tip labels from the tree
    std::vector<std::string> tips = value->getTipNames();

    // create a new data object
    HomologousDiscreteCharacterData<NaturalNumbersState> *tip_data = new HomologousDiscreteCharacterData<NaturalNumbersState>();
    for (size_t i = 0; i < tips.size(); i++)
    {
        DiscreteTaxonData<NaturalNumbersState> this_tip_data = DiscreteTaxonData<NaturalNumbersState>(tips[i]);
        NaturalNumbersState state = NaturalNumbersState(0, num_states);
        state.setState("?");
        this_tip_data.addCharacter(state);
        tip_data->addTaxonData(this_tip_data);
    }

    // store the new value for the discrete data
    static_cast<TreeDiscreteCharacterData*>(this->value)->setCharacterData(tip_data);

    // update the kernel
    updateData(true);

    // mark the likelihood dirty
    probability_dirty = true;

}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::keepSpecialization(const DagNode* affecter)
{

    if ( affecter == age )
    {
        dag_node->keepAffected();
    }

    // clear all flags
    for (std::vector<bool>::iterator it = this->dirty_nodes.begin(); it != this->dirty_nodes.end(); ++it)
    {
        (*it) = false;
    }

    probability_dirty    = false;
    tree_dirty           = false;
    root_frequency_dirty = false;
    lambda_dirty         = false;
    mu_dirty             = false;
    phi_dirty            = false;
    delta_dirty          = false;
    upsilon_dirty        = false;
    gamma_dirty          = false;
    rho_dirty            = false;
    xi_dirty             = false;
    eta_dirty            = false;
    omega_dirty          = false;
    zeta_dirty           = false;
    old_ln_prob          = current_ln_prob; // make sure we don't accidently restore the outdated likelihood

    // tell tp to keep the approximator
    if (tp_can_reset) {
    	tp_ptr->keepSchedulerAndApproximator();
    	tp_can_reset = false;
    }

}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::prepareParameters(bool force)
{

	// make sure all the parameters are up-to-date
	updateTree(force);
	updateData(force);
	updateRootFrequency(force);
	updateLambda(force);
	updateMu(force);
	updatePhi(force);
	updateDelta(force);
	updateUpsilon(force);
	updateGamma(force);
	updateRho(force);
	updateXi(force);
	updateEta(force);
	updateOmega(force);
	updateZeta(force);

}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::restoreSpecialization(const DagNode *restorer)
{

    // reset the flags
    for (std::vector<bool>::iterator it = dirty_nodes.begin(); it != dirty_nodes.end(); ++it)
    {
        (*it) = false;
    }

    if ( restorer == age )
    {
        if ( use_origin == false)
        {
            value->getRoot().setAge( age->getValue() );
        }

        if ( dag_node != NULL )
        {
            dag_node->touchAffected();
        }

//        // reset the tree
//        if (tp_can_reset) {
//        	updateTree(true);
//        }

    }
    
    // MJL added 230307
    if ( restorer == this->dag_node )
    {

//        // reset the tree
//        if (tp_can_reset) {
//        	updateTree(true);
//        }

        // make sure we update the likelihood
//        probability_dirty = true;

    }

    if ( restorer != this->dag_node )
    {
    	if ( restorer == root_frequency )
    	{
    		root_frequency_dirty = true;
    	}
    	else if ( restorer == lambda_const || restorer == lambda_var || restorer == lambda_times  )
    	{
    		lambda_dirty = true;
    	}
    	else if ( restorer == mu_const || restorer == mu_var || restorer == mu_times  )
    	{
    		mu_dirty = true;
    	}
    	else if ( restorer == phi_const || restorer == phi_var || restorer == phi_times  )
    	{
    		phi_dirty = true;
    	}
    	else if ( restorer == delta_const || restorer == delta_var || restorer == delta_times )
    	{
    		delta_dirty = true;
    	}
    	else if ( restorer == upsilon || restorer == upsilon_times  )
    	{
    		upsilon_dirty = true;
    	}
    	else if ( restorer == gamma || restorer == gamma_times )
    	{
    		gamma_dirty = true;
    	}
    	else if ( restorer == rho_simple || restorer == rho || restorer == rho_times  )
    	{
    		rho_dirty = true;
    	}
    	else if ( restorer == xi || restorer == xi_times  )
    	{
    		xi_dirty = true;
    	}
    	else if ( restorer == eta_simple || restorer == eta_const || restorer == eta_var || restorer == eta_times  )
    	{
    		eta_dirty = true;
    	}
    	else if ( restorer == omega_const || restorer == omega_var || restorer == omega_times  )
    	{
    		omega_dirty = true;
    	}
    	else if ( restorer == zeta  )
    	{
    		zeta_dirty = true;
    	}
    }

    // reset scheduler and approximator on tp side
    if ( tp_can_reset ) {
    	tp_ptr->resetSchedulerAndApproximator();
    	tp_can_reset = false;
    }

    // clean the likelihood
    probability_dirty = false;
    current_ln_prob   = old_ln_prob;

}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::simulateSSETree(void) {
    RBOUT("Warning: simulating tree under SSE model (not yet implemented).");
    
    // simulating a tree
    // std::string sim_newick_string = tp_ptr->simulateTree();
    // this->setValue( static_cast<Tree*>( tp_ptr->getTreeValue() ) );
    // this->setValue( static_cast<DiscreteCharacterData*>( tp_ptr->getDataValue() ) );
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::simulateTree(void)
{

	// Warning: simulating tree under uniform model.
	RBOUT("Warning: simulating tree under uniform model.");

    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;

    // Create the time tree object (topology + times)
    Tree *psi = new TreeDiscreteCharacterData();

    // Root the topology by setting the appropriate flag
    psi->setRooted( true );

    // make a vector of tip nodes
    std::vector<TopologyNode* > nodes;

    // set tip names
    size_t num_taxa = taxa.size();
    for (size_t i=0; i<num_taxa; ++i)
    {
        // get the node from the list
        TopologyNode* node = new TopologyNode(i);

        // set name and age
        const std::string& name = taxa[i].getName();
        node->setName(name);
        node->setSpeciesName(taxa[i].getSpeciesName());
        node->setAge(taxa[i].getAge());
        node->setTaxon(taxa[i]);

        // add to tips
        nodes.push_back(node);
    }

    // get times for simulation
    std::vector<double> ages = simulateCoalescentAges();

    // recursively build the tree
    buildSerialSampledRandomBinaryTree(psi, nodes, ages);

    // initialize the topology by setting the root
    TopologyNode* root = nodes[nodes.size() - 1];
    psi->setRoot(root, true);

    // store the new value for the tree
    value->getTreeChangeEventHandler().removeListener( this );
    value = psi;
    value->getTreeChangeEventHandler().addListener( this );

    // update the kernel
    updateTree(true);

    // set tip states to missing
    initializeEmptyCharData();

}

std::vector<double> GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::simulateCoalescentAges() const
{
    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;

    // create the container for the ages
    std::vector<double> ages;

    // for each tip, simulate an age between max age and tip age
    double max_age = age->getValue();
    size_t num_taxa = taxa.size();
    size_t num_ages;
    if ( use_origin == false )
    {
    	// root age provided
    	ages.push_back(max_age);
    	num_ages = num_taxa - 2;
    }
    else
    {
    	// stem age provided
    	num_ages = num_taxa - 1;
    }

    for(size_t i = 0; i < num_ages; ++i)
    {
    	// get the age of the tip
    	double a = taxa[i + 1].getAge();

    	// simulate the age of a node
    	double new_age = a + rng->uniform01() * (max_age - a);

    	//add the age to the vector of ages
    	ages.push_back(new_age);
    }

    // sort the ages (from youngest to oldest)
    std::sort(ages.begin(), ages.end(), std::greater<double>());

    // return the ages
    return ages;
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::touchSpecialization(const DagNode *affecter, bool touchAll)
{

    // we don't need to reset tp approximator yet
    tp_can_reset = false;

    if ( affecter == age )
    {
        if ( use_origin == false)
        {
            value->getRoot().setAge( age->getValue() );
        }

        if ( dag_node != NULL )
        {
            dag_node->touchAffected();
        }

        // update the tree
        tree_dirty = true;
//        updateTree();

        // make sure we update the likelihood
        probability_dirty = true;

    }
    
    // MJL added 230307
    if ( affecter == this->dag_node )
    {
        // update the tree
        tree_dirty = true;

        // make sure we update the likelihood
        probability_dirty = true;

    }

    if ( affecter != this->dag_node and affecter != age )
    {

    	// mark all nodes as dirty
        for (std::vector<bool>::iterator it = dirty_nodes.begin(); it != dirty_nodes.end(); ++it)
        {
            (*it) = true;
        }

        // mark the affecting parameter as dirty
        if ( affecter == root_frequency )
		{
			root_frequency_dirty = true;
			probability_dirty = true;
		}
        else if ( affecter == lambda_const || affecter == lambda_var || affecter == lambda_times  )
    	{
    		lambda_dirty = true;
    		probability_dirty = true;
    	}
    	else if ( affecter == mu_const || affecter == mu_var || affecter == mu_times  )
    	{
    		mu_dirty = true;
    		probability_dirty = true;
    	}
    	else if ( affecter == phi_const || affecter == phi_var || affecter == phi_times  )
    	{
    		phi_dirty = true;
    		probability_dirty = true;
    	}
    	else if ( affecter == delta_const || affecter == delta_var || affecter == delta_times )
    	{
    		delta_dirty = true;
    		probability_dirty = true;
    	}
    	else if ( affecter == upsilon || affecter == upsilon_times  )
    	{
    		upsilon_dirty = true;
    		probability_dirty = true;
    	}
    	else if ( affecter == gamma || affecter == gamma_times )
    	{
    		gamma_dirty = true;
    		probability_dirty = true;
    	}
    	else if ( affecter == rho_simple || affecter == rho || affecter == rho_times  )
    	{
    		rho_dirty = true;
    		probability_dirty = true;
    	}
    	else if ( affecter == xi || affecter == xi_times  )
    	{
    		xi_dirty = true;
    		probability_dirty = true;
    	}
    	else if ( affecter == eta_simple || affecter == eta_const || affecter == eta_var || affecter == eta_times  )
    	{
    		eta_dirty = true;
    		probability_dirty = true;
    	}
    	else if ( affecter == omega_const || affecter == omega_var || affecter == omega_times  )
    	{
    		omega_dirty = true;
    		probability_dirty = true;
    	}
    	else if ( affecter == zeta  )
    	{
    		zeta_dirty = true;
    		probability_dirty = true;
    	}
    }

}

std::vector<double> GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::RbToStd(const Simplex &obj)
{
	std::vector<double> std_object;
	for(size_t i = 0; i < obj.size(); ++i) {
		std_object.push_back(obj[i]);
	}
	return std_object;
}

std::vector<double> GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::RbToStd(const RbVector<double> &obj)
{
	std::vector<double> std_object;
	for(size_t i = 0; i < obj.size(); ++i) {
		std_object.push_back(obj[i]);
	}
	return std_object;
}

std::vector< std::vector<double> > GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::RbToStd(const RbVector< RbVector<double>> &obj)
{
	std::vector< std::vector<double> > std_object;
	for(size_t i = 0; i < obj.size(); ++i) {
		std::vector<double> sub;
		for(size_t j = 0; j < obj[i].size(); ++j) {
			sub.push_back(obj[i][j]);
		}
		std_object.push_back(sub);
	}
	return std_object;
}

std::vector< std::vector< std::vector<double> > > GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::RbToStd(const RbVector< RateGenerator > &obj)
{
	std::vector< std::vector< std::vector<double> > > std_object;
	for(size_t i = 0; i < obj.size(); ++i) {
		std::vector< std::vector<double> > sub;
		const RateGenerator &this_rate_generator = obj[i];
		for(size_t j = 0; j < num_states; ++j) {
			std::vector<double> sub_sub;
			for(size_t k = 0; k < num_states; ++k)
			{
				sub_sub.push_back( this_rate_generator.getRate(j, k, 0.0, 1.0) );
			}
			sub.push_back(sub_sub);
		}
		std_object.push_back(sub);
	}
	return std_object;
}

std::vector< std::vector< std::vector<double> > > GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::RbToStd(const RbVector< MatrixReal > &obj)
{
	std::vector< std::vector< std::vector<double> > > std_object;
	for(size_t i = 0; i < obj.size(); ++i) {
		std::vector< std::vector<double> > sub;
		const MatrixReal &this_matrix = obj[i];
		size_t nrow = this_matrix.getNumberOfRows();
		size_t ncol = this_matrix.getNumberOfColumns();
		for(size_t j = 0; j < nrow; ++j) {
			std::vector<double> sub_sub;
			for(size_t k = 0; k < ncol; ++k)
			{
				sub_sub.push_back( this_matrix[j][k] );
			}
			sub.push_back(sub_sub);
		}
		std_object.push_back(sub);
	}
	return std_object;
}

std::vector< std::map< std::vector<unsigned>, double > > GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::RbToStd(const RbVector< CladogeneticProbabilityMatrix > &obj)
{
	/* Original
	std::vector< std::map< std::vector<unsigned >, double >> std_object;
	for(size_t i = 0; i < obj.size(); ++i)
	{
		std_object.push_back( obj[i].getEventMap() );
	}*/
	// Copy and convert from 1-indexing to 0-indexing
	std::vector< std::map< std::vector<unsigned >, double >> std_object(obj.size());
	for(size_t i = 0; i < obj.size(); ++i)
	{
		for(std::map< std::vector<unsigned >, double >::const_iterator it=obj[i].getEventMap().begin(); it != obj[i].getEventMap().end(); ++it)
		{
			// Recover the map entries
			std::vector<unsigned> key = it->first;
			double val = it->second;
			// 1-indexing to 0-indexing
			if ( zero_indexed == false )
			{
				for(size_t iK=0; iK<key.size(); ++iK) key[iK] -= 1;
			}
			std_object[i][key] = val;
		}
	}
	return std_object;
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::checkTimesAreAscending(const RbVector<double> &obj) {

	bool isValid = true;
	for(size_t i = 0; i < obj.size(); ++i) {
		if ( std::fabs(obj[i] - 0.0) < std::numeric_limits<double>::epsilon() ) {
			isValid = false;
			break;
		}
		if ( i != obj.size() - 1 ) {
			if (obj[i] > obj[i + 1]) {
				isValid = false;
				break;
			}
		}
	}

	if (!isValid) {
		throw RbException("Problem with rate change times. Please make sure that times are in ascending order, and that there is no time that is zero. The first element of the time vector should be the start age of the youngest epoch, etc.");
	}

}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updateTree(bool force)
{
	if ( force or tree_dirty )
	{

		// MRM/BP 7/26/23: Bruno and I made a change here to collapse zero-length branches
		// into sampled ancestors before getting the newick string. This solves a problem
		// but may create a new one: might be expensive to traverse entire tree and collapse
		// ancestors every time we change the tree?
		// possible alternative: use a different algorithm to make the newick string that
		// correctly creates sampled ancestor nodes even when the branch has zero length

		// collapse sampled ancestors
		Tree* tmp_tree = this->getValue().clone();
		tmp_tree->collapseSampledAncestors();

		// get the newick string
		std::string var = tmp_tree->getNewickRepresentation();
//		std::string var = this->getValue().getNewickRepresentation();

		if ( use_origin )
		{
			// strip out trailing zeros
			char pattern = ':';
			while ( true )
			{
				// if we found a colon stop
				if ( var.back() == pattern )
				{
					break;
				}
				// otherwise, pop off the last character
				var.pop_back();
			}

			// now add the tail
			double origin_age  = age->getValue();
			double root_age    = this->getValue().getRoot().getAge();
			double tail_length = origin_age - root_age;
			var += std::to_string(tail_length);
		}
		else
		{
			// strip off the tail
			// strip out trailing zeros
            char pattern = ':';
			while ( true )
			{
				// if we found a colon stop
				if ( var.back() == pattern )
				{
					break;
				}
				// otherwise, pop off the last character
				var.pop_back();
			}

			// now add the trailing zeros
			var += "0.00000";
		}

		// make sure there's a closing semicolon
		var += ";";

		// set the tree
//		dirty_nodes = std::vector<bool>(dirty_nodes.size(), true);
		tp_ptr->setTree(var, dirty_nodes);
//		tp_ptr->setTree(var);

		// mark tree as clean
		tree_dirty = false;

	}
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updateData(bool force)
{
	if ( force )
	{
		// get the data
		const RevBayesCore::AbstractHomologousDiscreteCharacterData &dat = getCharacterData();

		size_t num_taxa = taxa.size();
		std::map< std::string, std::vector<double> > taxon_map;
		std::vector<std::string> taxon_names(num_taxa);
		for(size_t i = 0; i < num_taxa; ++i)
		{
			const AbstractDiscreteTaxonData& this_taxon_data = dat.getTaxonData(taxa[i].getName());

			// if the data are weighted, use the weights
			if ( this_taxon_data[0].isWeighted() )
			{
				taxon_map[taxa[i].getName()] = this_taxon_data[0].getWeights();
			}
			else
			{
				RbBitSet this_bit_set = this_taxon_data[0].getState();
				std::vector<double> this_taxon_states(num_states);
				for(size_t j = 0; j < num_states; ++j)
				{
					this_taxon_states[j] = this_bit_set[j] == true ? 1.0 : 0.0;
				}
				taxon_map[taxa[i].getName()] = this_taxon_states;
			}

			taxon_names[i] = taxa[i].getName();
		}

		// set the data
		tp_ptr->setData(taxon_names, taxon_map);

	}
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updateRootFrequency(bool force)
{
	if ( force or root_frequency_dirty )
	{
		// create empty vector
		std::vector<double> var;
		if ( root_frequency != NULL )
		{
			var = RbToStd( root_frequency->getValue() );
		}
		else
		{
			throw RbException("Root frequencies cannot be NULL.");
		}

		// set the parameters
		tp_ptr->setRootPrior(var);

	}

	// flag as clean
	root_frequency_dirty = false;

}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updateLambda(bool force)
{
	if ( force or lambda_dirty )
	{
		if ( lambda_const != NULL )
		{
			// create intermediate parameters
			RbVector< RbVector<double> > intermediate_parameters;
			intermediate_parameters.push_back( RbVector<double>( lambda_const->getValue() ) );

			// create empty vectors
			std::vector< std::vector<double> > params = RbToStd( intermediate_parameters );
			std::vector<double>                times;

			// set the parameters
			tp_ptr->setLambda(times, params);
		}
		else if ( lambda_var != NULL )
		{

			// check for ascending times
			checkTimesAreAscending(lambda_times->getValue());

			// create empty vectors
			std::vector< std::vector<double> > params = RbToStd( lambda_var->getValue() );
			std::vector<double>                times  = RbToStd( lambda_times->getValue() );

			// set the parameters
			tp_ptr->setLambda(times, params);

		}

	}

	// flag as clean
	lambda_dirty = false;

}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updateMu(bool force)
{
	if ( force or mu_dirty )
	{
		if ( mu_const != NULL )
		{
			// create intermediate parameters
			RbVector< RbVector<double> > intermediate_parameters;
			intermediate_parameters.push_back( RbVector<double>( mu_const->getValue() ) );

			// create empty vectors
			std::vector< std::vector<double> > params = RbToStd( intermediate_parameters );
			std::vector<double>                times;

			// set the parameters
			tp_ptr->setMu(times, params);
		}
		else if ( mu_var != NULL )
		{

			// check for ascending times
			checkTimesAreAscending(mu_times->getValue());

			// create empty vectors
			std::vector< std::vector<double> > params = RbToStd( mu_var->getValue() );
			std::vector<double>                times  = RbToStd( mu_times->getValue() );

			// set the parameters
			tp_ptr->setMu(times, params);
		}

	}

	// flag as clean
	mu_dirty = false;

}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updatePhi(bool force)
{
	if ( force or phi_dirty )
	{
		if ( phi_const != NULL )
		{
			// create intermediate parameters
			RbVector< RbVector<double> > intermediate_parameters;
			intermediate_parameters.push_back( RbVector<double>( phi_const->getValue() ) );

			// create empty vectors
			std::vector< std::vector<double> > params = RbToStd( intermediate_parameters );
			std::vector<double>                times;

			// set the parameters
			tp_ptr->setPhi(times, params);
		}
		else if ( phi_var != NULL )
		{

			// check for ascending times
			checkTimesAreAscending(phi_times->getValue());

			// create empty vectors
			std::vector< std::vector<double> > params = RbToStd( phi_var->getValue() );
			std::vector<double>                times  = RbToStd( phi_times->getValue() );

			// set the parameters
			tp_ptr->setPhi(times, params);
		}
	}

	// flag as clean
	phi_dirty = false;

}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updateDelta(bool force)
{
	if ( force or delta_dirty )
	{
		if ( delta_const != NULL )
		{
			// create intermediate parameters
			RbVector< RbVector<double> > intermediate_parameters;
			intermediate_parameters.push_back( RbVector<double>( delta_const->getValue() ) );

			// create empty vectors
			std::vector< std::vector<double> > params = RbToStd( intermediate_parameters );
			std::vector<double>                times;

			// set the parameters
			tp_ptr->setDelta(times, params);
		}
		else if ( delta_var != NULL )
		{

			// check for ascending times
			checkTimesAreAscending(delta_times->getValue());

			// create empty vectors
			std::vector< std::vector<double> > params = RbToStd( delta_var->getValue() );
			std::vector<double>                times  = RbToStd( delta_times->getValue() );

			// set the parameters
			tp_ptr->setDelta(times, params);
		}
	}

	// flag as clean
	delta_dirty = false;

}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updateUpsilon(bool force)
{
	if ( force or upsilon_dirty )
	{
		// create empty vectors
		std::vector< std::vector<double> > params;
		std::vector<double>                times;

		if ( upsilon != NULL && upsilon_times != NULL )
		{
			// convert to std
			params = RbToStd( upsilon->getValue() );
			times  = RbToStd( upsilon_times->getValue() );

			// set the parameters
			tp_ptr->setMassSpeciationEvents(times, params);
		}

	}

	// flag as clean
	upsilon_dirty = false;

}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updateGamma(bool force)
{
	if ( force or gamma_dirty )
	{
		// create empty vectors
		std::vector< std::vector<double> > params;
		std::vector<double>                times;

		if ( gamma != NULL && gamma_times != NULL )
		{
			// convert to std
			params = RbToStd( gamma->getValue() );
			times  = RbToStd( gamma_times->getValue() );

			// set the parameters
			tp_ptr->setMassExtinctionEvents(times, params);
		}

	}

	// flag as clean
	gamma_dirty = false;

}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updateRho(bool force)
{
	if ( force or rho_dirty )
	{
		if ( rho_simple != NULL )
		{
			// create empty vectors
			std::vector< std::vector<double> > params;
			std::vector<double>                times(1, 0.0);

			// replicate the value of rho for each state
			params.push_back( std::vector<double>( num_states, rho_simple->getValue() ) );

			// set the parameters
			tp_ptr->setMassSamplingEvents(times, params);
		}
		else if ( rho != NULL )
		{
			// create empty vectors
			std::vector< std::vector<double> > params;
			std::vector<double>                times;

			if ( rho != NULL && rho_times != NULL )
			{
				// convert to std
				params = RbToStd( rho->getValue() );
				times  = RbToStd( rho_times->getValue() );
			}

			// set the parameters
			tp_ptr->setMassSamplingEvents(times, params);
		}

	}

	// flag as clean
	rho_dirty = false;

}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updateXi(bool force)
{
	if ( force or xi_dirty )
	{
		// create empty vectors
		std::vector< std::vector<double> > params;
		std::vector<double>                times;

		if ( xi != NULL && xi_times != NULL )
		{
			// convert to std
			params = RbToStd( xi->getValue() );
			times  = RbToStd( xi_times->getValue() );

			// set the parameters
			tp_ptr->setMassDestrSamplingEvents(times, params);
		}

	}

	// flag as clean
	xi_dirty = false;

}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updateEta(bool force)
{
	if ( force or eta_dirty )
	{
		if ( eta_simple != NULL )
		{
			// create empty vectors
			std::vector< std::vector< std::vector<double> > > params;
			std::vector<double>                               times;

			// create the matrix
			std::vector< std::vector<double> > matrix;
			std::vector<double> row(num_states, eta_simple->getValue());
			double diagonal = -1.0 * eta_simple->getValue() * ((double)num_states - 1.0);
			for(size_t i = 0; i < num_states; ++i)
			{
				matrix.push_back( row );
				matrix[i][i] = diagonal;
			}

			// add the matrix to the vector
			params.push_back(matrix);

			// set the parameters
			tp_ptr->setEta(times, params);
		}
		else if ( eta_const != NULL )
		{
			// create intermediate parameters
			RbVector< RateGenerator > intermediate_parameters;
			intermediate_parameters.push_back( eta_const->getValue() );

			// create empty vectors
			std::vector< std::vector< std::vector<double> > > params = RbToStd( intermediate_parameters );
			std::vector<double>                               times;

			// set the parameters
			tp_ptr->setEta(times, params);
		}
		else if ( eta_var != NULL )
		{

			// check for ascending times
			checkTimesAreAscending(eta_times->getValue());

			// create empty vectors
			std::vector< std::vector< std::vector<double> > > params = RbToStd( eta_var->getValue() );
			std::vector<double>                               times  = RbToStd( eta_times->getValue() );

			// set the parameters
			tp_ptr->setEta(times, params);
		}

	}

	// flag as clean
	eta_dirty = false;

}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updateOmega(bool force)
{
	if ( force or omega_dirty )
	{
		if ( omega_const != NULL )
		{
			// create intermediate parameters
			RbVector< CladogeneticProbabilityMatrix > intermediate_parameters;
			intermediate_parameters.push_back( omega_const->getValue() );
			// create empty vectors
			std::vector< std::map< std::vector<unsigned>, double > > params = RbToStd( intermediate_parameters );
			std::vector<double>                                      times;
			// set the parameters
			tp_ptr->setOmega(num_states, times, params);
		}
		else if ( omega_var != NULL )
		{

			// check for ascending times
			checkTimesAreAscending(omega_times->getValue());

			// create empty vectors
			std::vector< std::map< std::vector<unsigned>, double > > params = RbToStd( omega_var->getValue() );
			std::vector<double>                                      times  = RbToStd( omega_times->getValue() );

			// set the parameters
			tp_ptr->setOmega(num_states, times, params);
		}

	}

	// flag as clean
	omega_dirty = false;

}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updateZeta(bool force)
{
	if ( force or zeta_dirty )
	{
		if ( zeta != NULL )
		{
			// convert to std
			std::vector< std::vector< std::vector<double> > > params;

			// convert to std
			params = RbToStd( zeta->getValue() );

			// set the parameters
			tp_ptr->setMassExtinctionStateChangeProb(params);
		}

	}

	// flag as clean
	zeta_dirty = false;

}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
	if ( oldP == age )
	{
		age = static_cast<const TypedDagNode<double>* >( newP );
	}
	if ( oldP == root_frequency )
	{
		root_frequency = static_cast<const TypedDagNode<Simplex>* >( newP );
		root_frequency_dirty = true;
	}
	if ( oldP == lambda_const )
	{
		lambda_const = static_cast<const TypedDagNode<RbVector<double>  >* >( newP );
		lambda_dirty = true;
	}
	if ( oldP == lambda_var )
	{
		lambda_var = static_cast<const TypedDagNode<RbVector<RbVector<double> > >* >( newP );
		lambda_dirty = true;
	}
	if ( oldP == lambda_times )
	{
		lambda_times = static_cast<const TypedDagNode<RbVector<double>  >* >( newP );
		lambda_dirty = true;
	}
	if ( oldP == mu_const )
	{
		mu_const = static_cast<const TypedDagNode<RbVector<double>  >* >( newP );
		mu_dirty = true;
	}
	if ( oldP == mu_var )
	{
		mu_var = static_cast<const TypedDagNode<RbVector<RbVector<double> > >* >( newP );
		mu_dirty = true;
	}
	if ( oldP == mu_times )
	{
		mu_times = static_cast<const TypedDagNode<RbVector<double>  >* >( newP );
		mu_dirty = true;
	}
	if ( oldP == phi_const )
	{
		phi_const = static_cast<const TypedDagNode<RbVector<double>  >* >( newP );
		phi_dirty = true;
	}
	if ( oldP == phi_var )
	{
		phi_var = static_cast<const TypedDagNode<RbVector<RbVector<double> > >* >( newP );
		phi_dirty = true;
	}
	if ( oldP == phi_times )
	{
		phi_times = static_cast<const TypedDagNode<RbVector<double>  >* >( newP );
		phi_dirty = true;
	}
	if ( oldP == delta_const )
	{
		delta_const = static_cast<const TypedDagNode<RbVector<double>  >* >( newP );
		delta_dirty = true;
	}
	if ( oldP == delta_var )
	{
		delta_var = static_cast<const TypedDagNode<RbVector<RbVector<double> > >* >( newP );
		delta_dirty = true;
	}
	if ( oldP == delta_times )
	{
		delta_times = static_cast<const TypedDagNode<RbVector<double>  >* >( newP );
		delta_dirty = true;
	}
	if ( oldP == upsilon )
	{
		upsilon = static_cast<const TypedDagNode<RbVector<RbVector<double> > >* >( newP );
		upsilon_dirty = true;
	}
	if ( oldP == upsilon_times )
	{
		upsilon_times = static_cast<const TypedDagNode<RbVector<double>  >* >( newP );
		upsilon_dirty = true;
	}
	if ( oldP == gamma )
	{
		gamma = static_cast<const TypedDagNode<RbVector<RbVector<double> > >* >( newP );
		gamma_dirty = true;
	}
	if ( oldP == gamma_times )
	{
		gamma_times = static_cast<const TypedDagNode<RbVector<double>  >* >( newP );
		gamma_dirty = true;
	}
	if ( oldP == rho_simple )
	{
		rho_simple = static_cast<const TypedDagNode<double>* >( newP );
		rho_dirty = true;
	}
	if ( oldP == rho )
	{
		rho = static_cast<const TypedDagNode<RbVector<RbVector<double> > >* >( newP );
		rho_dirty = true;
	}
	if ( oldP == rho_times )
	{
		rho_times = static_cast<const TypedDagNode<RbVector<double>  >* >( newP );
		rho_dirty = true;
	}
	if ( oldP == xi )
	{
		xi = static_cast<const TypedDagNode<RbVector<RbVector<double> > >* >( newP );
		xi_dirty = true;
	}
	if ( oldP == xi_times )
	{
		xi_times = static_cast<const TypedDagNode<RbVector<double>  >* >( newP );
		xi_dirty = true;
	}
	if ( oldP == eta_const )
	{
		eta_const = static_cast<const TypedDagNode<RateGenerator >* >( newP );
		eta_dirty = true;
	}
	if ( oldP == eta_simple )
	{
		eta_simple = static_cast<const TypedDagNode<double>* >( newP );
		eta_dirty = true;
	}
	if ( oldP == eta_var )
	{
		eta_var = static_cast<const TypedDagNode<RbVector<RateGenerator > >* >( newP );
		eta_dirty = true;
	}
	if ( oldP == eta_times )
	{
		eta_times = static_cast<const TypedDagNode<RbVector<double>  >* >( newP );
		eta_dirty = true;
	}
	if ( oldP == omega_const )
	{
		omega_const = static_cast<const TypedDagNode<CladogeneticProbabilityMatrix >* >( newP );
		omega_dirty = true;
	}
	if ( oldP == omega_var )
	{
		omega_var = static_cast<const TypedDagNode<RbVector<CladogeneticProbabilityMatrix > >* >( newP );
		omega_dirty = true;
	}
	if ( oldP == omega_times )
	{
		omega_times = static_cast<const TypedDagNode<RbVector<double>  >* >( newP );
		omega_dirty = true;
	}
	if ( oldP == zeta )
	{
		zeta = static_cast<const TypedDagNode<RbVector<MatrixReal > >* >( newP );
		zeta_dirty = true;
	}
}


void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<double> &rv) const
{

}

RevLanguage::RevPtr<RevLanguage::RevVariable> GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::executeProcedure(const std::string &name, const std::vector<DagNode *> args, bool &found)
{
    if (name == "clampCharData")
    {
        found = true;

    	// get the incoming data
        const AbstractHomologousDiscreteCharacterData& v = static_cast<const TypedDagNode<AbstractHomologousDiscreteCharacterData > *>( args[0] )->getValue();

        // check the number of states
        if(v.getNumberOfStates() != num_states)
        {
        	throw RbException("Number of discrete states in the clamped data does not match the number of discrete states for the distribution.");
        }

        // check the taxon labels in the new data are consistent
        std::vector<std::string> param_taxa;
        size_t num_taxa = taxa.size();

        std::vector<Taxon> input_tax = v.getTaxa();
        std::vector<std::string> input_taxa(input_tax.size());
        for(size_t i = 0; i < input_tax.size(); ++i)
        {
        	input_taxa[i] = input_tax[i].getName();
        }

        // check that the number of taxa match
        if ( input_taxa.size() != num_taxa )
        {
        	throw RbException("Number of taxa in the character data does not match number of taxa in the tree.");
        }

        // sort the labels
        std::sort(param_taxa.begin(), param_taxa.end());
        std::sort(input_taxa.begin(), input_taxa.end());
        for(size_t i = 0; i < num_taxa; ++i)
        {
        	if ( input_taxa[i] != input_taxa[i] )
        	{
        		throw RbException("Taxon labels in the new character data do not match with taxon labels in the tree.");
        	}
        }

        // copy the input data
        static_cast<TreeDiscreteCharacterData*>(this->value)->setCharacterData( v.clone() );

        // update the kernel
        updateData(true);

        // mark the likelihood dirty
        probability_dirty = true;

        return NULL;
    }

    if (name == "getCharData")
    {
        found = true;
        RevLanguage::AbstractHomologousDiscreteCharacterData *tip_states = new RevLanguage::AbstractHomologousDiscreteCharacterData( getCharacterData() );
        return new RevLanguage::RevVariable( tip_states );
    }

    if (name == "dumpModel")
    {
    	found = true;

    	// get the filename
    	std::string file_name = args[0]->getValueAsString();

    	// write the file
    	tp_ptr->writeStateToFile(file_name);

    	return NULL;
    }

	return TypedDistribution<Tree>::executeProcedure( name, args, found );
}
