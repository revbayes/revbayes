#include "GeneralizedLineageHeterogeneousBirthDeathSamplingProcess.h"


using namespace RevBayesCore;

/**
 * Constructor.
*/
GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::GeneralizedLineageHeterogeneousBirthDeathSamplingProcess() : TypedDistribution<Tree>( new TreeDiscreteCharacterData() )
{

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

}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::redrawValue(void)
{

}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::setValue(Tree *v, bool f)
{

}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::getAffected(RbOrderedSet<DagNode *>& affected, DagNode* affecter)
{

}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::keepSpecialization(DagNode* affecter)
{

}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::restoreSpecialization(DagNode *restorer)
{

}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::touchSpecialization(DagNode *toucher, bool touchAll)
{

}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{

}

void GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<double> &rv) const
{

}

RevLanguage::RevPtr<RevLanguage::RevVariable> GeneralizedLineageHeterogeneousBirthDeathSamplingProcess::executeProcedure(const std::string &name, const std::vector<DagNode *> args, bool &found)
{

}



















