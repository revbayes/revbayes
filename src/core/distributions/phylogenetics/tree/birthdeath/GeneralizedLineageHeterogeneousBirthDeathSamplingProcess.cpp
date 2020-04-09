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
	size_t                                                          n_proc_
) : TypedDistribution<Tree>( new TreeDiscreteCharacterData() ),
	current_ln_prob(0.0),
	old_ln_prob(0.0),
	probability_dirty(true),
	n_proc(n_proc_),
	taxa(taxa_),
	age(age_),
	condition_type(condition_type_),
	num_states(num_states_),
	root_frequency(root_frequency_),
	lambda_const(NULL),
	lambda_var(NULL),
	lambda_times(NULL),
	mu_const(NULL),
	mu_var(NULL),
	mu_times(NULL),
	phi_const(NULL),
	phi_var(NULL),
	phi_times(NULL),
	delta_const(NULL),
	delta_var(NULL),
	delta_times(NULL),
	upsilon(NULL),
	upsilon_times(NULL),
	gamma(NULL),
	gamma_times(NULL),
	rho_simple(NULL),
	rho(NULL),
	rho_times(NULL),
	xi(NULL),
	xi_times(NULL),
	eta_simple(NULL),
	eta_const(NULL),
	eta_var(NULL),
	eta_times(NULL),
	omega_const(NULL),
	omega_var(NULL),
	omega_times(NULL),
	zeta(NULL),
	use_origin(use_origin_),
	tree_dirty(true),
	root_freq_dirty(true),
	lambda_dirty(true),
	mu_dirty(true),
	phi_dirty(true),
	delta_dirty(true),
	upsilon_dirty(true),
	gamma_dirty(true),
	rho_dirty(true),
	xi_dirty(true),
	eta_dirty(true),
	omega_dirty(true),
	zeta_dirty(true)
{

	//assert(Plugin::loader().isTensorPhyloLoaded());
	try {
		tp_ptr = Plugin::loader().createTensorPhyloLik();
	} catch (...) {
		throw RbException("TensorPhylo is not loaded (use loadPlugin(...)).");
	}
//	std::cout << "tensorphylo version: " << tp_ptr->getVersion() << std::endl;
	// turn on/off debug
	tp_ptr->setDebugMode(TensorPhylo::Interface::DBG_FILE, "debug.txt");
	tp_ptr->setDebugMode(TensorPhylo::Interface::DBG_PRINT);
	tp_ptr->setConditionalProbCompatibilityMode(false); // FIXME Here you go Mike!

	// add the parameters
	addParameter(age);
	addParameter(root_frequency);

	// add this object to the tree change event handler
    value->getTreeChangeEventHandler().addListener( this );

    // simulate the tree
    simulateTree();

    // update the kernel
    updateRootFrequency(true);


//    tp_ptr->setApplyTreeLikCorrection(false);

    // set the condition type
    tp_ptr->setConditionalProbabilityType(TensorPhylo::Interface::TIME);
    if ( condition_type != "time" )
    {
        if ( use_origin )
        {
        	if ( condition_type == "survival" | condition_type == "sampledExtant" )
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
        	if ( condition_type == "survival" | condition_type == "sampledExtant" )
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
		parent->setNodeType( false, false, true );
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

	// store the old likelihood
    old_ln_prob = current_ln_prob;

    // calculate a likelihood!
//	tp_ptr->writeStateToFile("params.dat");
		//tp_ptr->loadStateFromFile("state.dat");
    current_ln_prob = tp_ptr->computeLogLikelihood();

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


void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::fireTreeChangeEvent(const TopologyNode &n, const unsigned& m)
{
	tree_dirty = true;
	updateTree();
}

const RevBayesCore::AbstractHomologousDiscreteCharacterData& GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::getCharacterData() const
{
    return static_cast<TreeDiscreteCharacterData*>(this->value)->getCharacterData();
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::redrawValue(void)
{
	// just simulate with a uniform tree and all missing data
    simulateTree();
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

	// dispatch an update
	updateLambda(true);

	// mark clean
	lambda_dirty = false;
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::setLambda(const TypedDagNode< RbVector< RbVector<double> > >* param, const TypedDagNode< RbVector<double> >* times)
{
	if ( lambda_const != NULL || lambda_var != NULL )
	{
		throw RbException("Tried to set lambda twice.");
	}

	// set the value
	lambda_var   = param;
	lambda_times = times;

	// include the parameter
	addParameter(lambda_var);
	addParameter(lambda_times);

	// dispatch an update
	updateLambda(true);

	// mark clean
	lambda_dirty = false;
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

	// dispatch an update
	updateMu(true);

	// mark clean
	mu_dirty = false;
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::setMu(const TypedDagNode< RbVector< RbVector<double> > >* param, const TypedDagNode< RbVector<double> >* times)
{
	if ( mu_const != NULL || mu_var != NULL )
	{
		throw RbException("Tried to set mu twice.");
	}

	// set the value
	mu_var   = param;
	mu_times = times;

	// include the parameter
	addParameter(mu_var);
	addParameter(mu_times);

	// dispatch an update
	updateMu(true);

	// mark clean
	mu_dirty = false;
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

	// dispatch an update
	updatePhi(true);

	// mark clean
	phi_dirty = false;
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::setPhi(const TypedDagNode< RbVector< RbVector<double> > >* param, const TypedDagNode< RbVector<double> >* times)
{
	if ( phi_const != NULL || phi_var != NULL )
	{
		throw RbException("Tried to set phi twice.");
	}

	// set the value
	phi_var   = param;
	phi_times = times;

	// include the parameter
	addParameter(phi_var);
	addParameter(phi_times);

	// dispatch an update
	updatePhi(true);

	// mark clean
	phi_dirty = false;
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

	// dispatch an update
	updateDelta(true);

	// mark clean
	delta_dirty = false;
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::setDelta(const TypedDagNode< RbVector< RbVector<double> > >* param, const TypedDagNode< RbVector<double> >* times)
{
	if ( delta_const != NULL || delta_var != NULL )
	{
		throw RbException("Tried to set delta twice.");
	}

	// set the value
	delta_var   = param;
	delta_times = times;

	// include the parameter
	addParameter(delta_var);
	addParameter(delta_times);

	// dispatch an update
	updateDelta(true);

	// mark clean
	delta_dirty = false;
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

	// dispatch an update
	updateUpsilon(true);

	// mark clean
	upsilon_dirty = false;
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

	// dispatch an update
	updateGamma(true);

	// mark clean
	gamma_dirty = false;
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

	// dispatch an update
	updateRho(true);

	// mark clean
	rho_dirty = false;
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

	// dispatch an update
	updateRho(true);

	// mark clean
	rho_dirty = false;
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

	// dispatch an update
	updateXi(true);

	// mark clean
	xi_dirty = false;
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

	// dispatch an update
	updateEta(true);

	// mark clean
	eta_dirty = false;
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

	// dispatch an update
	updateEta(true);

	// mark clean
	eta_dirty = false;
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::setEta(const TypedDagNode< RbVector< RateGenerator > >* param, const TypedDagNode< RbVector<double> >* times)
{
	if ( eta_simple != NULL ||  eta_const != NULL || eta_var != NULL )
	{
		throw RbException("Tried to set eta twice.");
	}

	// set the value
	eta_var   = param;
	eta_times = times;

	// include the parameter
	addParameter(eta_var);
	addParameter(eta_times);

	// dispatch an update
	updateEta(true);

	// mark clean
	eta_dirty = false;
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

	// dispatch an update
	updateOmega(true);

	// mark clean
	omega_dirty = false;
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::setOmega(const TypedDagNode< RbVector< CladogeneticProbabilityMatrix > >* param, const TypedDagNode< RbVector<double> >* times)
{
	if ( omega_const != NULL || omega_var != NULL )
	{
		throw RbException("Tried to set omega twice.");
	}

	// set the value
	omega_var   = param;
	omega_times = times;

	// include the parameter
	addParameter(omega_var);
	addParameter(omega_times);

	// dispatch an update
	updateOmega(true);

	// mark clean
	omega_dirty = false;
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

	// dispatch an update
	updateZeta(true);

	// mark clean
	zeta_dirty = false;
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::setValue(Tree *v, bool f)
{
	// check that tree is binary
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

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::keepSpecialization(DagNode* affecter)
{

    if ( affecter == age )
    {
        dag_node->keepAffected();
    }

    // clear all flags
    probability_dirty = false;
    tree_dirty        = false;
    root_freq_dirty   = false;
    lambda_dirty      = false;
    mu_dirty          = false;
    phi_dirty         = false;
    delta_dirty       = false;
    upsilon_dirty     = false;
    gamma_dirty       = false;
    rho_dirty         = false;
    xi_dirty          = false;
    eta_dirty         = false;
    omega_dirty       = false;
    zeta_dirty        = false;
    old_ln_prob       = current_ln_prob; // make sure we don't accidently restore the outdated likelihood

}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::prepareParameters(bool force)
{

//	std::cout << "Preparing parameters" << std::endl;

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

//	std::cout << "Done preparing parameters" << std::endl;

}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::restoreSpecialization(DagNode *restorer)
{

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

        // update the tree
        updateTree();
        tree_dirty = false;
    }

    if ( restorer != this->dag_node )
    {
    	if ( restorer == lambda_const || restorer == lambda_var || restorer == lambda_times  )
    	{
    		updateLambda();
    		lambda_dirty = false;
    	}
    	else if ( restorer == mu_const || restorer == mu_var || restorer == mu_times  )
    	{
    		updateMu();
    		mu_dirty = false;
    	}
    	else if ( restorer == phi_const || restorer == phi_var || restorer == phi_times  )
    	{
    		updatePhi();
    		phi_dirty = false;
    	}
    	else if ( restorer == delta_const || restorer == delta_var || restorer == delta_times )
    	{
    		updateDelta();
    		delta_dirty = false;
    	}
    	else if ( restorer == upsilon || restorer == upsilon_times  )
    	{
    		updateUpsilon();
    		upsilon_dirty = false;
    	}
    	else if ( restorer == gamma || restorer == gamma_times )
    	{
    		updateGamma();
    		gamma_dirty = false;
    	}
    	else if ( restorer == rho_simple || restorer == rho || restorer == rho_times  )
    	{
    		updateRho();
    		rho_dirty = false;
    	}
    	else if ( restorer == xi || restorer == xi_times  )
    	{
    		updateXi();
    		xi_dirty = false;
    	}
    	else if ( restorer == eta_simple || restorer == eta_const || restorer == eta_var || restorer == eta_times  )
    	{
    		updateEta();
    		eta_dirty = false;
    	}
    	else if ( restorer == omega_const || restorer == omega_var || restorer == omega_times  )
    	{
    		updateOmega();
    		omega_dirty = false;
    	}
    	else if ( restorer == zeta  )
    	{
    		updateZeta();
    		omega_dirty = false;
    	}
    }

    // clean the likelihood
    probability_dirty = false;
    current_ln_prob   = old_ln_prob;

}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::simulateTree(void)
{
	// Warning: simulating tree under uniform model.
	RBOUT("Warning: simulating tree under uniform model.");

	// simulating a tree

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
        node->setNodeType( true, false, false );
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

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::touchSpecialization(DagNode *affecter, bool touchAll)
{

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
        updateTree();

        // make sure we update the likelihood
        probability_dirty = true;

    }

    if ( affecter != this->dag_node )
    {
    	if ( affecter == lambda_const || affecter == lambda_var || affecter == lambda_times  )
    	{
    		lambda_dirty = true;
    		updateLambda();
    		probability_dirty = true;
    	}
    	else if ( affecter == mu_const || affecter == mu_var || affecter == mu_times  )
    	{
    		mu_dirty = true;
    		updateMu();
    		probability_dirty = true;
    	}
    	else if ( affecter == phi_const || affecter == phi_var || affecter == phi_times  )
    	{
    		phi_dirty = true;
    		updatePhi();
    		probability_dirty = true;
    	}
    	else if ( affecter == delta_const || affecter == delta_var || affecter == delta_times )
    	{
    		delta_dirty = true;
    		updateDelta();
    		probability_dirty = true;
    	}
    	else if ( affecter == upsilon || affecter == upsilon_times  )
    	{
    		upsilon_dirty = true;
    		updateUpsilon();
    		probability_dirty = true;
    	}
    	else if ( affecter == gamma || affecter == gamma_times )
    	{
    		gamma_dirty = true;
    		updateGamma();
    		probability_dirty = true;
    	}
    	else if ( affecter == rho_simple || affecter == rho || affecter == rho_times  )
    	{
    		rho_dirty = true;
    		updateRho();
    		probability_dirty = true;
    	}
    	else if ( affecter == xi || affecter == xi_times  )
    	{
    		xi_dirty = true;
    		updateXi();
    		probability_dirty = true;
    	}
    	else if ( affecter == eta_simple || affecter == eta_const || affecter == eta_var || affecter == eta_times  )
    	{
    		eta_dirty = true;
    		updateEta();
    		probability_dirty = true;
    	}
    	else if ( affecter == omega_const || affecter == omega_var || affecter == omega_times  )
    	{
    		omega_dirty = true;
    		updateOmega();
    		probability_dirty = true;
    	}
    	else if ( affecter == zeta  )
    	{
    		zeta_dirty = true;
    		updateZeta();
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
			for(size_t iK=0; iK<key.size(); ++iK) key[iK] -= 1;
			std_object[i][key] = val;
		}
	}
	return std_object;
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updateTree(bool force)
{
	if ( force | tree_dirty )
	{
		// get the newick string
		std::string var = this->getValue().getNewickRepresentation();

		if ( use_origin )
		{
			// strip out trailing zeros
			std::string pattern = ":";
			while ( true )
			{
				// if we found a colon stop
				if ( &var.back() == pattern )
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
			std::string pattern = ":";
			while ( true )
			{
				// if we found a colon stop
				if ( &var.back() == pattern )
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

		// std::cout << var << std::endl;

		// set the tree
		tp_ptr->setTree(var);
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
	if ( force | root_freq_dirty )
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
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updateLambda(bool force)
{
	if ( force | lambda_dirty )
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
			// create empty vectors
			std::vector< std::vector<double> > params = RbToStd( lambda_var->getValue() );
			std::vector<double>                times  = RbToStd( lambda_times->getValue() );

			// handle some errors
			if ( params.size() != times.size() + 1 )
			{
				throw RbException( "Number of lambda vectors does not match the number of intervals." );
			}

			// set the parameters
			tp_ptr->setLambda(times, params);
		}
	}
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updateMu(bool force)
{
	if ( force | mu_dirty )
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
			// create empty vectors
			std::vector< std::vector<double> > params = RbToStd( mu_var->getValue() );
			std::vector<double>                times  = RbToStd( mu_times->getValue() );

			// handle some errors
			if ( params.size() != times.size() + 1 )
			{
				throw RbException( "Number of mu vectors does not match the number of intervals." );
			}

			// set the parameters
			tp_ptr->setMu(times, params);
		}
	}
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updatePhi(bool force)
{
	if ( force | phi_dirty )
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
			// create empty vectors
			std::vector< std::vector<double> > params = RbToStd( phi_var->getValue() );
			std::vector<double>                times  = RbToStd( phi_times->getValue() );

			// handle some errors
			if ( params.size() != times.size() + 1 )
			{
				throw RbException( "Number of phi vectors does not match the number of intervals." );
			}

			// set the parameters
			tp_ptr->setPhi(times, params);
		}
	}
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updateDelta(bool force)
{
	if ( force | delta_dirty )
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
			// create empty vectors
			std::vector< std::vector<double> > params = RbToStd( delta_var->getValue() );
			std::vector<double>                times  = RbToStd( delta_times->getValue() );

			// handle some errors
			if ( params.size() != times.size() + 1 )
			{
				throw RbException( "Number of delta vectors does not match the number of intervals." );
			}

			// set the parameters
			tp_ptr->setDelta(times, params);
		}
	}

}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updateUpsilon(bool force)
{
	if ( force | upsilon_dirty )
	{
		// create empty vectors
		std::vector< std::vector<double> > params;
		std::vector<double>                times;

		if ( upsilon != NULL && upsilon_times != NULL )
		{
			// convert to std
			params = RbToStd( upsilon->getValue() );
			times  = RbToStd( upsilon_times->getValue() );

			// handle some errors
			if ( params.size() != times.size() )
			{
				throw RbException( "Number of upsilon events does not match the number of times." );
			}

			// set the parameters
			tp_ptr->setMassSpeciationEvents(times, params);
		}
	}
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updateGamma(bool force)
{
	if ( force | gamma_dirty )
	{
		// create empty vectors
		std::vector< std::vector<double> > params;
		std::vector<double>                times;

		if ( gamma != NULL && gamma_times != NULL )
		{
			// convert to std
			params = RbToStd( gamma->getValue() );
			times  = RbToStd( gamma_times->getValue() );

			// handle some errors
			if ( params.size() != times.size() )
			{
				throw RbException( "Number of gamma events does not match the number of times." );
			}

			// set the parameters
			tp_ptr->setMassExtinctionEvents(times, params);
		}
	}
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updateRho(bool force)
{
	if ( force | rho_dirty )
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

				// handle some errors
				if ( params.size() != times.size() )
				{
					throw RbException( "Number of rho events does not match the number of times." );
				}
			}

			// set the parameters
			tp_ptr->setMassSamplingEvents(times, params);
		}
	}
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updateXi(bool force)
{
	if ( force | xi_dirty )
	{
		// create empty vectors
		std::vector< std::vector<double> > params;
		std::vector<double>                times;

		if ( xi != NULL && xi_times != NULL )
		{
			// convert to std
			params = RbToStd( xi->getValue() );
			times  = RbToStd( xi_times->getValue() );

			// handle some errors
			if ( params.size() != times.size() )
			{
				throw RbException( "Number of xi events does not match the number of times." );
			}

			// set the parameters
			tp_ptr->setMassDestrSamplingEvents(times, params);
		}

	}
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updateEta(bool force)
{
	if ( force | eta_dirty )
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
			// create empty vectors
			std::vector< std::vector< std::vector<double> > > params = RbToStd( eta_var->getValue() );
			std::vector<double>                               times  = RbToStd( eta_times->getValue() );

			// handle some errors
			if ( params.size() != times.size() + 1 )
			{
				throw RbException( "Number of eta matrices does not match the number of intervals." );
			}

			// set the parameters
			tp_ptr->setEta(times, params);
		}
	}
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updateOmega(bool force)
{
	if ( force | omega_dirty )
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
			// create empty vectors
			std::vector< std::map< std::vector<unsigned>, double > > params = RbToStd( omega_var->getValue() );
			std::vector<double>                                      times  = RbToStd( omega_times->getValue() );

			// handle some errors
			if ( params.size() != times.size() + 1 )
			{
				throw RbException( "Number of omega matrices does not match the number of intervals." );
			}

			// set the parameters
			tp_ptr->setOmega(num_states, times, params);
		}
	}
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updateZeta(bool force)
{
	if ( force | zeta_dirty )
	{
		if ( zeta != NULL )
		{
			// convert to std
			std::vector< std::vector< std::vector<double> > > params;

			// convert to std
			params = RbToStd( zeta->getValue() );

			// handle some errors
			if ( params.size() != gamma_times->getValue().size() )
			{
				throw RbException( "Number of zeta matrices does not match the number of times." );
			}

			// set the parameters
			tp_ptr->setMassExtinctionStateChangeProb(params);
		}

	}
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
	}
	if ( oldP == lambda_const )
	{
		lambda_const = static_cast<const TypedDagNode<RbVector<double>  >* >( newP );
	}
	if ( oldP == lambda_var )
	{
		lambda_var = static_cast<const TypedDagNode<RbVector<RbVector<double> > >* >( newP );
	}
	if ( oldP == lambda_times )
	{
		lambda_times = static_cast<const TypedDagNode<RbVector<double>  >* >( newP );
	}
	if ( oldP == mu_const )
	{
		mu_const = static_cast<const TypedDagNode<RbVector<double>  >* >( newP );
	}
	if ( oldP == mu_var )
	{
		mu_var = static_cast<const TypedDagNode<RbVector<RbVector<double> > >* >( newP );
	}
	if ( oldP == mu_times )
	{
		mu_times = static_cast<const TypedDagNode<RbVector<double>  >* >( newP );
	}
	if ( oldP == phi_const )
	{
		phi_const = static_cast<const TypedDagNode<RbVector<double>  >* >( newP );
	}
	if ( oldP == phi_var )
	{
		phi_var = static_cast<const TypedDagNode<RbVector<RbVector<double> > >* >( newP );
	}
	if ( oldP == phi_times )
	{
		phi_times = static_cast<const TypedDagNode<RbVector<double>  >* >( newP );
	}
	if ( oldP == delta_const )
	{
		delta_const = static_cast<const TypedDagNode<RbVector<double>  >* >( newP );
	}
	if ( oldP == delta_var )
	{
		delta_var = static_cast<const TypedDagNode<RbVector<RbVector<double> > >* >( newP );
	}
	if ( oldP == delta_times )
	{
		delta_times = static_cast<const TypedDagNode<RbVector<double>  >* >( newP );
	}
	if ( oldP == upsilon )
	{
		upsilon = static_cast<const TypedDagNode<RbVector<RbVector<double> > >* >( newP );
	}
	if ( oldP == upsilon_times )
	{
		upsilon_times = static_cast<const TypedDagNode<RbVector<double>  >* >( newP );
	}
	if ( oldP == gamma )
	{
		gamma = static_cast<const TypedDagNode<RbVector<RbVector<double> > >* >( newP );
	}
	if ( oldP == gamma_times )
	{
		gamma_times = static_cast<const TypedDagNode<RbVector<double>  >* >( newP );
	}
	if ( oldP == rho_simple )
	{
		rho_simple = static_cast<const TypedDagNode<double>* >( newP );
	}
	if ( oldP == rho )
	{
		rho = static_cast<const TypedDagNode<RbVector<RbVector<double> > >* >( newP );
	}
	if ( oldP == rho_times )
	{
		rho_times = static_cast<const TypedDagNode<RbVector<double>  >* >( newP );
	}
	if ( oldP == xi )
	{
		xi = static_cast<const TypedDagNode<RbVector<RbVector<double> > >* >( newP );
	}
	if ( oldP == xi_times )
	{
		xi_times = static_cast<const TypedDagNode<RbVector<double>  >* >( newP );
	}
	if ( oldP == eta_const )
	{
		eta_const = static_cast<const TypedDagNode<RateGenerator >* >( newP );
	}
	if ( oldP == eta_simple )
	{
		eta_simple = static_cast<const TypedDagNode<double>* >( newP );
	}
	if ( oldP == eta_var )
	{
		eta_var = static_cast<const TypedDagNode<RbVector<RateGenerator > >* >( newP );
	}
	if ( oldP == eta_times )
	{
		eta_times = static_cast<const TypedDagNode<RbVector<double>  >* >( newP );
	}
	if ( oldP == omega_const )
	{
		omega_const = static_cast<const TypedDagNode<CladogeneticProbabilityMatrix >* >( newP );
	}
	if ( oldP == omega_var )
	{
		omega_var = static_cast<const TypedDagNode<RbVector<CladogeneticProbabilityMatrix > >* >( newP );
	}
	if ( oldP == omega_times )
	{
		omega_times = static_cast<const TypedDagNode<RbVector<double>  >* >( newP );
	}
	if ( oldP == zeta )
	{
		zeta = static_cast<const TypedDagNode<RbVector<MatrixReal > >* >( newP );
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
