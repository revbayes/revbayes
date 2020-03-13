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

using namespace RevBayesCore;

/**
 * Constructor.
*/
GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::GeneralizedLineageHeterogeneousBirthDeathSamplingProcess(
	const std::vector<Taxon>&                                       taxa_,
	const TypedDagNode<double>*                                     age_,
	const std::string&                                              condition_type_,
	const TypedDagNode<Simplex >*                                   root_frequency_,
	const TypedDagNode<RbVector< RbVector< double > > >*            lambda_,
	const TypedDagNode<RbVector< double > >*                        lambda_times_,
	const TypedDagNode<RbVector< RbVector< double > > >*            mu_,
	const TypedDagNode<RbVector< double > >*                        mu_times_,
	const TypedDagNode<RbVector< RbVector< double > > >*            phi_,
	const TypedDagNode<RbVector< double > >*                        phi_times_,
	const TypedDagNode<RbVector< RbVector< double > > >*            delta_,
	const TypedDagNode<RbVector< double > >*                        delta_times_,
	const TypedDagNode<RbVector< RbVector< double > > >*            upsilon_,
	const TypedDagNode<RbVector< double > >*                        upsilon_times_,
	const TypedDagNode<RbVector< RbVector< double > > >*            gamma_,
	const TypedDagNode<RbVector< double > >*                        gamma_times_,
	const TypedDagNode<RbVector< RbVector< double > > >*            rho_,
	const TypedDagNode<RbVector< double > >*                        rho_times_,
	const TypedDagNode<RbVector< RbVector< double > > >*            xi_,
	const TypedDagNode<RbVector< double > >*                        xi_times_,
	const TypedDagNode<RbVector< RateGenerator > >*                 eta_,
	const TypedDagNode<RbVector< double > >*                        eta_times_,
	const TypedDagNode<RbVector< CladogeneticProbabilityMatrix > >* omega_,
	const TypedDagNode<RbVector< double > >*                        omega_times_,
	const TypedDagNode<RbVector< MatrixReal > >*                    zeta_,
	bool                                                            use_origin_
) : TypedDistribution<Tree>( new TreeDiscreteCharacterData() ),
	current_ln_prob(0.0),
	old_ln_prob(0.0),
	probability_dirty(true),
	taxa(taxa_),
	age(age_),
	condition_type(condition_type_),
	root_frequency(root_frequency_),
	lambda(lambda_),
	lambda_times(lambda_times_),
	mu(mu_),
	mu_times(mu_times_),
	phi(phi_),
	phi_times(phi_times_),
	delta(delta_),
	delta_times(delta_times_),
	upsilon(upsilon_),
	upsilon_times(upsilon_times_),
	gamma(gamma_),
	gamma_times(gamma_times_),
	rho(rho_),
	rho_times(rho_times_),
	xi(xi_),
	xi_times(xi_times_),
	eta(eta_),
	eta_times(eta_times_),
	omega(omega_),
	omega_times(omega_times_),
	zeta(zeta_),
	use_origin(use_origin_),
	tree_dirty(false),
	root_freq_dirty(false),
	lambda_dirty(false),
	mu_dirty(false),
	phi_dirty(false),
	delta_dirty(false),
	upsilon_dirty(false),
	gamma_dirty(false),
	rho_dirty(false),
	xi_dirty(false),
	eta_dirty(false),
	omega_dirty(false),
	zeta_dirty(false)
{

	//assert(Plugin::loader().isTensorPhyloLoaded());
	try {
		tp_ptr = Plugin::loader().createTensorPhyloLik();
	} catch (...) {
		throw RbException("TensorPhylo is not loaded (use loadPlugin(...)).");
	}
//	std::cout << "tensorphylo version: " << tp_ptr->getVersion() << std::endl;

	// add the parameters
	addParameter(age);
	addParameter(root_frequency);
	addParameter(lambda);
	addParameter(lambda_times);
	addParameter(mu);
	addParameter(mu_times);
	addParameter(phi);
	addParameter(phi_times);
	addParameter(delta);
	addParameter(delta_times);
	addParameter(upsilon);
	addParameter(upsilon_times);
	addParameter(gamma);
	addParameter(gamma_times);
	addParameter(rho);
	addParameter(rho_times);
	addParameter(xi);
	addParameter(xi_times);
	addParameter(eta);
	addParameter(eta_times);
	addParameter(omega);
	addParameter(omega_times);
	addParameter(zeta);

	// add this object to the tree change event handler
    value->getTreeChangeEventHandler().addListener( this );

    // simulate the tree
    simulateTree();

    // update the parameters
    prepareParameters(true);

    // turn on/off debug
//    tp_ptr->setDebugMode(TensorPhylo::Interface::DBG_PRINT);

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

		// make sure ages are consistent

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
    current_ln_prob = tp_ptr->computeLogLikelihood();

    // flag the likelihood as up-to-date
    probability_dirty = false;

    // NOTE: this likelihood differs from the one computed in other revbayes
    // distributions because we include the probability density of the root node.
    // that probability density drops out if we condition on anything but time, in which
    // case this likelihood will be equivalent to a likelihood from other revbayes
    // birth-death models

    return current_ln_prob;
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::fireTreeChangeEvent(const TopologyNode &n, const unsigned& m)
{
	// nothing smart to do yet
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
	// FIXME the memory bug seems to be coming from this line or the next one (double delete) and most probably happen in setTree (i.e., delete root;)
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

    // get the number of states in the model
    size_t num_states = lambda->getValue()[0].size();

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
    	if ( restorer == lambda || restorer == lambda_times  )
    	{
    		updateLambda();
    		lambda_dirty = false;
    	}
    	else if ( restorer == mu || restorer == mu_times  )
    	{
    		updateMu();
    		mu_dirty = false;
    	}
    	else if ( restorer == phi || restorer == phi_times  )
    	{
    		updatePhi();
    		phi_dirty = false;
    	}
    	else if ( restorer == delta || restorer == delta_times )
    	{
    		updateDelta();
    		delta_dirty = false;
    	}
    	else if ( restorer == upsilon || restorer == upsilon_times  )
    	{
    		updateUpsilon();
    		upsilon_dirty = false;
    	}
    	else if ( restorer == gamma || restorer == gamma_times || restorer == zeta )
    	{
    		updateGamma();
    		updateZeta();
    		gamma_dirty = false;
    		zeta_dirty = false;
    	}
    	else if ( restorer == rho || restorer == rho_times  )
    	{
    		updateRho();
    		rho_dirty = false;
    	}
    	else if ( restorer == xi || restorer == xi_times  )
    	{
    		updateXi();
    		xi_dirty = false;
    	}
    	else if ( restorer == eta || restorer == eta_times  )
    	{
    		updateEta();
    		eta_dirty = false;
    	}
    	else if ( restorer == omega || restorer == omega_times  )
    	{
    		updateOmega();
    		omega_dirty = false;
    	}
    }

    // clean the likelihood
    probability_dirty = false;
    current_ln_prob   = old_ln_prob;

}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::simulateTree(void)
{
	// throw a warning that this is not a correct simulation
	// throw RbException(RbException::DEFAULT, "Warning: simulating tree under uniform model.");

	// simulating a tree
    delete value;

    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;

    // Create the time tree object (topology + times)
    Tree *psi = new Tree();

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
    value = psi;

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
    	if ( affecter == lambda || affecter == lambda_times  )
    	{
    		lambda_dirty = true;
    		updateLambda();
    		probability_dirty = true;
    	}
    	else if ( affecter == mu || affecter == mu_times  )
    	{
    		mu_dirty = true;
    		updateMu();
    		probability_dirty = true;
    	}
    	else if ( affecter == phi || affecter == phi_times  )
    	{
    		phi_dirty = true;
    		updatePhi();
    		probability_dirty = true;
    	}
    	else if ( affecter == delta || affecter == delta_times )
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
    	else if ( affecter == gamma || affecter == gamma_times || affecter == zeta )
    	{
    		gamma_dirty = true;
    		zeta_dirty = true;
    		updateGamma();
    		updateZeta();
    		probability_dirty = true;
    	}
    	else if ( affecter == rho || affecter == rho_times  )
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
    	else if ( affecter == eta || affecter == eta_times  )
    	{
    		eta_dirty = true;
    		updateEta();
    		probability_dirty = true;
    	}
    	else if ( affecter == omega || affecter == omega_times  )
    	{
    		omega_dirty = true;
    		updateOmega();
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
		size_t num_states = this_rate_generator.getNumberOfStates();
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
			double origin_age  = age->getValue();
			double root_age    = this->getValue().getRoot().getAge();
			double tail_length = origin_age - root_age;
			// TODO: append the tail to the rooted tree
		}

		// make sure there's a closing semicolon
		var += ";";

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
			RbBitSet this_bit_set = this_taxon_data[0].getState();
			size_t num_states = this_bit_set.size();
			std::vector<double> this_taxon_states(num_states);
			for(size_t j = 0; j < num_states; ++j)
			{
				this_taxon_states[j] = this_bit_set[j] == true ? 1.0 : 0.0;
			}
			taxon_map[taxa[i].getName()] = this_taxon_states;

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
		std::vector<double> var = RbToStd( root_frequency->getValue() );

		// set the parameters
		tp_ptr->setRootPrior(var);
	}
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updateLambda(bool force)
{
	if ( force | lambda_dirty )
	{
		// convert to std
		std::vector< std::vector<double> > params = RbToStd( lambda->getValue() );
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

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updateMu(bool force)
{
	if ( force | mu_dirty )
	{
		// convert to std
		std::vector< std::vector<double> > params = RbToStd( mu->getValue() );
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

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updatePhi(bool force)
{
	if ( force | phi_dirty )
	{
		// convert to std
		std::vector< std::vector<double> > params = RbToStd( phi->getValue() );
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

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updateDelta(bool force)
{
	if ( force | delta_dirty )
	{
		// convert to std
		std::vector< std::vector<double> > params = RbToStd( delta->getValue() );
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

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updateUpsilon(bool force)
{
	if ( force | upsilon_dirty )
	{
		// convert to std
		std::vector< std::vector<double> > params = RbToStd( upsilon->getValue() );
		std::vector<double>                times  = RbToStd( upsilon_times->getValue() );

		// handle some errors
		if ( params.size() != times.size() )
		{
			throw RbException( "Number of upsilon vectors does not match the number of intervals." );
		}

		// set the parameters
		tp_ptr->setMassSpeciationEvents(times, params);
	}
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updateGamma(bool force)
{
	if ( force | gamma_dirty )
	{
		// convert to std
		std::vector< std::vector<double> > params = RbToStd( gamma->getValue() );
		std::vector<double>                times  = RbToStd( gamma_times->getValue() );

		// handle some errors
		if ( params.size() != times.size() )
		{
			throw RbException( "Number of gamma vectors does not match the number of intervals." );
		}

		// set the parameters
		tp_ptr->setMassExtinctionEvents(times, params);
	}
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updateRho(bool force)
{
	if ( force | rho_dirty )
	{
		// convert to std
		std::vector< std::vector<double> > params = RbToStd( rho->getValue() );
		std::vector<double>                times  = RbToStd( rho_times->getValue() );

		// handle some errors
		if ( params.size() != times.size() )
		{
			throw RbException( "Number of rho vectors does not match the number of intervals." );
		}

		// set the parameters
		tp_ptr->setMassSamplingEvents(times, params);
	}
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updateXi(bool force)
{
	if ( force | xi_dirty )
	{
		// convert to std
		std::vector< std::vector<double> > params = RbToStd( xi->getValue() );
		std::vector<double>                times  = RbToStd( xi_times->getValue() );

		// handle some errors
		if ( params.size() != times.size() )
		{
			throw RbException( "Number of xi vectors does not match the number of intervals." );
		}

		// set the parameters
		tp_ptr->setMassDestrSamplingEvents(times, params);
	}
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updateEta(bool force)
{
	if ( force | eta_dirty )
	{
		// convert to std
		std::vector< std::vector< std::vector<double> > > params = RbToStd( eta->getValue() );
		std::vector<double>                               times  = RbToStd( eta_times->getValue() );

		// handle some errors
		if ( params.size() != times.size() + 1 )
		{
			throw RbException( "Number of eta vectors does not match the number of intervals." );
		}

		// set the parameters
		tp_ptr->setEta(times, params);
	}
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updateOmega(bool force)
{
	if ( force | omega_dirty )
	{
		// convert to std
		std::vector< std::map< std::vector<unsigned>, double > > params = RbToStd( omega->getValue() );
		std::vector<double>                                      times  = RbToStd( omega_times->getValue() );

		const AbstractHomologousDiscreteCharacterData& char_data = getCharacterData();
		size_t num_states = char_data.getNumberOfStates();

		// handle some errors
		if ( params.size() != times.size() + 1 )
		{
			throw RbException( "Number of omega vectors does not match the number of intervals." );
		}

		// set the parameters
		tp_ptr->setOmega(num_states, times, params);
	}

}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updateZeta(bool force)
{
	if ( force | zeta_dirty )
	{
		// convert to std
		std::vector< std::vector< std::vector<double> > > params = RbToStd( zeta->getValue() );

		// handle some errors
		if ( params.size() != gamma_times->getValue().size() )
		{
			throw RbException( "Number of zeta vectors does not match the number of intervals." );
		}

		// set the parameters
		tp_ptr->setMassExtinctionStateChangeProb(params);
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
	if ( oldP == lambda )
	{
		lambda = static_cast<const TypedDagNode<RbVector<RbVector<double> > >* >( newP );
	}
	if ( oldP == lambda_times )
	{
		lambda_times = static_cast<const TypedDagNode<RbVector<double>  >* >( newP );
	}
	if ( oldP == mu )
	{
		mu = static_cast<const TypedDagNode<RbVector<RbVector<double> > >* >( newP );
	}
	if ( oldP == mu_times )
	{
		mu_times = static_cast<const TypedDagNode<RbVector<double>  >* >( newP );
	}
	if ( oldP == phi )
	{
		phi = static_cast<const TypedDagNode<RbVector<RbVector<double> > >* >( newP );
	}
	if ( oldP == phi_times )
	{
		phi_times = static_cast<const TypedDagNode<RbVector<double>  >* >( newP );
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
	if ( oldP == eta )
	{
		eta = static_cast<const TypedDagNode<RbVector<RateGenerator > >* >( newP );
	}
	if ( oldP == eta_times )
	{
		eta_times = static_cast<const TypedDagNode<RbVector<double>  >* >( newP );
	}
	if ( oldP == omega )
	{
		omega = static_cast<const TypedDagNode<RbVector<CladogeneticProbabilityMatrix > >* >( newP );
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

	return TypedDistribution<Tree>::executeProcedure( name, args, found );
}
