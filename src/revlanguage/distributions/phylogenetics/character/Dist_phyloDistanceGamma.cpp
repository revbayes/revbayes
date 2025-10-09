#include "Dist_phyloDistanceGamma.h"

#include <cstddef>
#include <ostream>

#include "PhyloDistanceGamma.h"
#include "RlString.h"
#include "RlTree.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ModelObject.h"
#include "ModelVector.h"
#include "RbVector.h"
#include "TypeSpec.h"

namespace RevBayesCore { class Tree; }

using namespace RevLanguage;


Dist_phyloDistanceGamma::Dist_phyloDistanceGamma() : TypedDistribution< DistanceMatrix >()
{
    
}


Dist_phyloDistanceGamma::~Dist_phyloDistanceGamma()
{
    
}



Dist_phyloDistanceGamma* Dist_phyloDistanceGamma::clone( void ) const
{
    
    return new Dist_phyloDistanceGamma(*this);
}


RevBayesCore::TypedDistribution< RevBayesCore::DistanceMatrix >* Dist_phyloDistanceGamma::createDistribution( void ) const
{
    
    // get the parameters tau, nam, varianceNode and distanceNode that will be used to create the actual distribution.
    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* tau = static_cast<const Tree &>( tree->getRevObject() ).getDagNode();
    
    const std::vector<std::string>      &nam  = static_cast<const ModelVector<RlString> &>( names->getRevObject() ).getDagNode()->getValue();
    
    RevBayesCore::TypedDagNode< RevBayesCore::DistanceMatrix >* varianceNode = NULL;
    varianceNode = static_cast<const DistanceMatrix &>( varianceMatrix->getRevObject() ).getDagNode();
    
    RevBayesCore::TypedDagNode< RevBayesCore::DistanceMatrix >* distanceNode = NULL;
    distanceNode = static_cast<const DistanceMatrix &>( distanceMatrix->getRevObject() ).getDagNode();
    
    RevBayesCore::PhyloDistanceGamma* d = new RevBayesCore::PhyloDistanceGamma( tau );
    
    d->setNames( nam );
    d->setVarianceMatrix( varianceNode );
    d->setDistanceMatrix( distanceNode );
    
    d->redrawValue();
    
    return d;
}



/* Get Rev type of object */
const std::string& Dist_phyloDistanceGamma::getClassType(void)
{
    
    static std::string rev_type = "Dist_phyloDistanceGamma";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_phyloDistanceGamma::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution< DistanceMatrix >::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_phyloDistanceGamma::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "PhyloDistanceGamma";
    
    return d_name;
}


/** Return member rules (no members) */
const MemberRules& Dist_phyloDistanceGamma::getParameterRules(void) const
{
    
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        dist_member_rules.push_back( new ArgumentRule( "tree"             , Tree::getClassTypeSpec()                  , "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "distanceMatrix"   , DistanceMatrix::getClassTypeSpec()      , "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "varianceMatrix"   , DistanceMatrix::getClassTypeSpec()      , "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "names"            , ModelVector<RlString>::getClassTypeSpec() , "", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return dist_member_rules;
}


const TypeSpec& Dist_phyloDistanceGamma::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}


/** Print value for user */
void Dist_phyloDistanceGamma::printValue(std::ostream& o) const
{
    
    o << "Distance-Matrix-Generation-Along-Tree Process(tree=";
    if ( tree != NULL ) {
        o << tree->getName();
    } else {
        o << "?";
    }
    o << ", distanceMatrix=";
    if ( distanceMatrix != NULL ) {
        o << distanceMatrix->getName();
    } else {
        o << "?";
    }
    o << ", varianceMatrix=";
    if ( varianceMatrix != NULL ) {
        o << varianceMatrix->getName();
    } else {
        o << "?";
    }
    o << ", names=";
    if ( names != NULL ) {
        o << names->getName();
    } else {
        o << "?";
    }
    o << ")";
    
}


/** Set a member variable */
void Dist_phyloDistanceGamma::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "tree" )
    {
        tree = var;
    }
    else if ( name == "distanceMatrix" )
    {
        distanceMatrix = var;
    }
    else if ( name == "varianceMatrix" )
    {
        varianceMatrix = var;
    }
    else if ( name == "names" )
    {
        names = var;
    }
    else
    {
        TypedDistribution< DistanceMatrix >::setConstParameter(name, var);
    }
    
}

