#include <string>

#include "AbstractHomologousDiscreteCharacterData.h"
#include "GeneralizedLineageHeterogeneousBirthDeathSamplingProcess.h"
#include "RlAbstractHomologousDiscreteCharacterData.h"
#include "RlString.h"
#include "TopologyNode.h"
#include "Taxon.h"
#include "Tree.h"

using namespace RevBayesCore;

/**
 * Constructor.
*/
GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::GeneralizedLineageHeterogeneousBirthDeathSamplingProcess(
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


double GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::computeLnProbability(void)
{


	return 0.0;
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::fireTreeChangeEvent(const TopologyNode &n, const unsigned& m)
{
	// TODO: update the tensorphylo object with the new tree
}

const RevBayesCore::AbstractHomologousDiscreteCharacterData& GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::getCharacterData() const
{
    return static_cast<TreeDiscreteCharacterData*>(this->value)->getCharacterData();
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::redrawValue(void)
{

}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::setValue(Tree *v, bool f)
{

    if (v->isBinary() == false)
    {
        throw RbException("The process is only implemented for binary trees.");
    }
    value->getTreeChangeEventHandler().removeListener( this );

    // set the tree
    static_cast<TreeDiscreteCharacterData *>(this->value)->setTree( *v );

    // clear memory
    delete v;

    // add this object to the tree change event handler
    value->getTreeChangeEventHandler().addListener( this );

    // TODO: update the tensorphylo object with the new tree
    // we should probability initialize the tip data to ?


}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::getAffected(RbOrderedSet<DagNode *>& affected, DagNode* affecter)
{
    if ( affecter == age )
    {
        dag_node->initiateGetAffectedNodes( affected );
    }
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::keepSpecialization(DagNode* affecter)
{

    if ( affecter == age )
    {
        dag_node->keepAffected();
    }

}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::prepareParameters(bool force)
{
	// make sure all the parametes are up-to-date
	updateTree(force);
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



void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::restoreSpecialization(DagNode *restorer)
{

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
    }

    if ( affecter != this->dag_node )
    {

    	// TODO: update tensorphlyo object

    }

}

std::vector<double> GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::RbToStd(const Simplex &obj)
{
	std::vector<double> std_object;
	for(size_t i = 0; i < std_object.size(); ++i) {
		std_object.push_back(obj[i]);
	}
	return std_object;
}

std::vector<double> GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::RbToStd(const RbVector<double> &obj)
{
	std::vector<double> std_object;
	for(size_t i = 0; i < std_object.size(); ++i) {
		std_object.push_back(obj[i]);
	}
	return std_object;
}

std::vector< std::vector<double> > GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::RbToStd(const RbVector< RbVector<double>> &obj)
{
	std::vector< std::vector<double> > std_object;
	for(size_t i = 0; i < std_object.size(); ++i) {
		std::vector<double> sub;
		for(size_t j = 0; j < std_object[i].size(); ++j) {
			sub.push_back(obj[i][j]);
		}
		std_object.push_back(sub);
	}
	return std_object;
}

std::vector< std::vector< std::vector<double> > > GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::RbToStd(const RbVector< RateGenerator > &obj)
{
	std::vector< std::vector< std::vector<double> > > std_object;
	for(size_t i = 0; i < std_object.size(); ++i) {
		std::vector< std::vector<double> > sub;
		const RateGenerator &this_rate_generator = obj[i];
		for(size_t j = 0; j < std_object[i].size(); ++j) {
			std::vector<double> sub_sub;
			for(size_t k = 0; k < std_object[i][j].size(); ++k)
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
	for(size_t i = 0; i < std_object.size(); ++i) {
		std::vector< std::vector<double> > sub;
		const MatrixReal &this_matrix = obj[i];
		for(size_t j = 0; j < std_object[i].size(); ++j) {
			std::vector<double> sub_sub;
			for(size_t k = 0; k < std_object[i][j].size(); ++k)
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
	std::vector< std::map< std::vector<unsigned >, double >> std_object;
	for(size_t i = 0; i < std_object.size(); ++i)
	{
		std_object.push_back( obj[i].getEventMap() );
	}
	return std_object;
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updateTree(bool force)
{
	if ( force | tree_dirty )
	{
		std::string var = this->getValue().getNewickRepresentation();
		// TODO: send var to tensorphylo
	}
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updateRootFrequency(bool force)
{
	if ( force | root_freq_dirty )
	{
		std::vector<double> var = RbToStd( root_frequency->getValue() );
		// TODO: send var to tensorphylo
	}
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updateLambda(bool force)
{
	if ( force | lambda_dirty )
	{
		std::vector< std::vector<double> > params = RbToStd( lambda->getValue() );
		std::vector<double>                times  = RbToStd( lambda_times->getValue() );
		// TODO: send var to tensorphylo
	}
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updateMu(bool force)
{
	if ( force | mu_dirty )
	{
		std::vector< std::vector<double> > params = RbToStd( mu->getValue() );
		std::vector<double>                times  = RbToStd( mu_times->getValue() );
		// TODO: send var to tensorphylo
	}
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updatePhi(bool force)
{
	if ( force | phi_dirty )
	{
		std::vector< std::vector<double> > params = RbToStd( phi->getValue() );
		std::vector<double>                times  = RbToStd( phi_times->getValue() );
		// TODO: send var to tensorphylo
	}
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updateDelta(bool force)
{
	if ( force | delta_dirty )
	{
		std::vector< std::vector<double> > params = RbToStd( delta->getValue() );
		std::vector<double>                times  = RbToStd( delta_times->getValue() );
		// TODO: send var to tensorphylo
	}
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updateUpsilon(bool force)
{
	if ( force | upsilon_dirty )
	{
		std::vector< std::vector<double> > params = RbToStd( upsilon->getValue() );
		std::vector<double>                times  = RbToStd( upsilon_times->getValue() );
		// TODO: send var to tensorphylo
	}
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updateGamma(bool force)
{
	if ( force | gamma_dirty )
	{
		std::vector< std::vector<double> > params = RbToStd( gamma->getValue() );
		std::vector<double>                times  = RbToStd( gamma_times->getValue() );
		// TODO: send var to tensorphylo
	}
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updateRho(bool force)
{
	if ( force | rho_dirty )
	{
		std::vector< std::vector<double> > params = RbToStd( rho->getValue() );
		std::vector<double>                times  = RbToStd( rho_times->getValue() );
		// TODO: send var to tensorphylo
	}
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updateXi(bool force)
{
	if ( force | xi_dirty )
	{
		std::vector< std::vector<double> > params = RbToStd( xi->getValue() );
		std::vector<double>                times  = RbToStd( xi_times->getValue() );
		// TODO: send var to tensorphylo
	}
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updateEta(bool force)
{
	if ( force | mu_dirty )
	{
		std::vector< std::vector< std::vector<double> > > params = RbToStd( eta->getValue() );
		std::vector<double>                               times  = RbToStd( eta_times->getValue() );
		// TODO: send var to tensorphylo
	}
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updateOmega(bool force)
{
	if ( force | mu_dirty )
	{
		std::vector< std::map< std::vector<unsigned>, double > > params = RbToStd( omega->getValue() );
		std::vector<double>                                      times  = RbToStd( omega_times->getValue() );
		// TODO: send var to tensorphylo
	}
}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::updateZeta(bool force)
{
	if ( force | mu_dirty )
	{
		std::vector< std::vector< std::vector<double> > > params = RbToStd( zeta->getValue() );
		// TODO: send var to tensorphylo
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

        const AbstractHomologousDiscreteCharacterData& v = static_cast<const TypedDagNode<AbstractHomologousDiscreteCharacterData > *>( args[0] )->getValue();

        // check if the tip names match
        bool match = true;
        std::vector<std::string> tips = value->getTipNames();
        for (size_t i = 0; i < tips.size(); i++)
        {
            found = false;
            for (size_t j = 0; j < v.getNumberOfTaxa(); j++)
            {
                if (tips[i] == v[j].getTaxonName())
                {
                    found = true;
                    break;
                }
            }
            if (found == false)
            {
                match = false;
                break;
            }
        }
        if (match == false)
        {
            throw RbException("To clamp a character data object all taxa present in the tree must be present in the character data.");
        }

        static_cast<TreeDiscreteCharacterData*>(this->value)->setCharacterData( v.clone() );

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



















