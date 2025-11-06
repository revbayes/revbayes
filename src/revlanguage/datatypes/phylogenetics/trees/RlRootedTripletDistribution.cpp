#include <cstddef>
#include <sstream>
#include <string>
#include <vector>

#include "ModelVector.h"
#include "Natural.h"
#include "RlRootedTripletDistribution.h"
#include "RlString.h"
#include "RlTaxon.h"
#include "RlTimeTree.h"
#include "RootedTripletDistribution.h"
#include "Taxon.h"
#include "TypeSpec.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ConstantNode.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "MemberProcedure.h"
#include "MethodTable.h"
#include "ModelObject.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlConstantNode.h"
#include "RlUtils.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

using namespace RevLanguage;

/** Default constructor */
RootedTripletDistribution::RootedTripletDistribution(void) : ModelObject<RevBayesCore::RootedTripletDistribution>()
{
    ArgumentRules* nTripletsArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "nTriplets", Natural::getClassTypeSpec(),         nTripletsArgRules   ) );
    
    ArgumentRules* nTreesArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "nTrees", Natural::getClassTypeSpec(),          nTreesArgRules   ) );
    
    ArgumentRules* SpeciesArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "species", ModelVector<RlString>::getClassTypeSpec(), SpeciesArgRules   ) );
    
    ArgumentRules* TaxaArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "taxa", ModelVector<Taxon>::getClassTypeSpec(), TaxaArgRules ) );
    
    ArgumentRules* setTaxaArgRules = new ArgumentRules();
    setTaxaArgRules->push_back( new ArgumentRule( "taxa", ModelVector<Taxon>::getClassTypeSpec(), "", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    methods.addFunction( new MemberProcedure( "setTaxa", RlUtils::Void, setTaxaArgRules ) );
    
    ArgumentRules* setSpeciesArgRules = new ArgumentRules();
    setSpeciesArgRules->push_back( new ArgumentRule( "species", ModelVector<RlString>::getClassTypeSpec(), "", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    methods.addFunction( new MemberProcedure( "setSpecies", RlUtils::Void, setSpeciesArgRules ) );
    
    ArgumentRules* extractTripletsArgRules = new ArgumentRules();
    extractTripletsArgRules->push_back( new ArgumentRule( "trees", ModelVector<TimeTree>::getClassTypeSpec(), "", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    methods.addFunction( new MemberProcedure( "extractTriplets", RlUtils::Void,            extractTripletsArgRules  ) );
    
}

/** Construct from bool */
RootedTripletDistribution::RootedTripletDistribution(RevBayesCore::RootedTripletDistribution *t) : ModelObject<RevBayesCore::RootedTripletDistribution>( t )
{
    ArgumentRules* nTripletsArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "nTriplets", Natural::getClassTypeSpec(),         nTripletsArgRules   ) );
    
    ArgumentRules* nTreesArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "nTrees", Natural::getClassTypeSpec(),          nTreesArgRules   ) );
    
    ArgumentRules* SpeciesArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "species", ModelVector<RlString>::getClassTypeSpec(), SpeciesArgRules   ) );
    
    ArgumentRules* TaxaArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "taxa", ModelVector<Taxon>::getClassTypeSpec(), TaxaArgRules ) );
    
    ArgumentRules* setTaxaArgRules = new ArgumentRules();
    setTaxaArgRules->push_back( new ArgumentRule( "taxa", ModelVector<Taxon>::getClassTypeSpec(), "", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    methods.addFunction( new MemberProcedure( "setTaxa", RlUtils::Void, setTaxaArgRules ) );
    
    ArgumentRules* setSpeciesArgRules = new ArgumentRules();
    setSpeciesArgRules->push_back( new ArgumentRule( "species", ModelVector<RlString>::getClassTypeSpec(), "", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    methods.addFunction( new MemberProcedure( "setSpecies", RlUtils::Void, setSpeciesArgRules ) );
    
    ArgumentRules* extractTripletsArgRules = new ArgumentRules();
    extractTripletsArgRules->push_back( new ArgumentRule( "trees", ModelVector<TimeTree>::getClassTypeSpec(), "", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    methods.addFunction( new MemberProcedure( "extractTriplets", RlUtils::Void,            extractTripletsArgRules  ) );
    
}

/** Construct from bool */
RootedTripletDistribution::RootedTripletDistribution(const RevBayesCore::RootedTripletDistribution &t) : ModelObject<RevBayesCore::RootedTripletDistribution>( new RevBayesCore::RootedTripletDistribution( t ) )
{
    ArgumentRules* nTripletsArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "nTriplets", Natural::getClassTypeSpec(),         nTripletsArgRules   ) );
    
    ArgumentRules* nTreesArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "nTrees", Natural::getClassTypeSpec(),          nTreesArgRules   ) );
    
    ArgumentRules* SpeciesArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "species", ModelVector<RlString>::getClassTypeSpec(), SpeciesArgRules   ) );
    
    ArgumentRules* TaxaArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "taxa", ModelVector<Taxon>::getClassTypeSpec(), TaxaArgRules ) );
    
    ArgumentRules* setTaxaArgRules = new ArgumentRules();
    setTaxaArgRules->push_back( new ArgumentRule( "taxa", ModelVector<Taxon>::getClassTypeSpec(), "", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    methods.addFunction( new MemberProcedure( "setTaxa", RlUtils::Void, setTaxaArgRules ) );
    
    ArgumentRules* setSpeciesArgRules = new ArgumentRules();
    setSpeciesArgRules->push_back( new ArgumentRule( "species", ModelVector<RlString>::getClassTypeSpec(), "", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    methods.addFunction( new MemberProcedure( "setSpecies", RlUtils::Void, setSpeciesArgRules ) );
    
    ArgumentRules* extractTripletsArgRules = new ArgumentRules();
    extractTripletsArgRules->push_back( new ArgumentRule( "trees", ModelVector<TimeTree>::getClassTypeSpec(), "", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    methods.addFunction( new MemberProcedure( "extractTriplets", RlUtils::Void,            extractTripletsArgRules  ) );
    
}

/** Construct from bool */
RootedTripletDistribution::RootedTripletDistribution(RevBayesCore::TypedDagNode<RevBayesCore::RootedTripletDistribution> *n) : ModelObject<RevBayesCore::RootedTripletDistribution>( n )
{
    ArgumentRules* nTripletsArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "nTriplets", Natural::getClassTypeSpec(),         nTripletsArgRules   ) );
    
    ArgumentRules* nTreesArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "nTrees", Natural::getClassTypeSpec(),          nTreesArgRules   ) );
    
    ArgumentRules* SpeciesArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "species", ModelVector<RlString>::getClassTypeSpec(), SpeciesArgRules   ) );
    
    ArgumentRules* TaxaArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "taxa", ModelVector<Taxon>::getClassTypeSpec(), TaxaArgRules ) );
    
    ArgumentRules* setTaxaArgRules = new ArgumentRules();
    setTaxaArgRules->push_back( new ArgumentRule( "taxa", ModelVector<Taxon>::getClassTypeSpec(), "", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    methods.addFunction( new MemberProcedure( "setTaxa", RlUtils::Void, setTaxaArgRules ) );
    
    ArgumentRules* setSpeciesArgRules = new ArgumentRules();
    setSpeciesArgRules->push_back( new ArgumentRule( "species", ModelVector<RlString>::getClassTypeSpec(), "", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    methods.addFunction( new MemberProcedure( "setSpecies", RlUtils::Void, setSpeciesArgRules ) );
    
    ArgumentRules* extractTripletsArgRules = new ArgumentRules();
    extractTripletsArgRules->push_back( new ArgumentRule( "trees", ModelVector<TimeTree>::getClassTypeSpec(), "", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    methods.addFunction( new MemberProcedure( "extractTriplets", RlUtils::Void,            extractTripletsArgRules  ) );
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
RootedTripletDistribution* RootedTripletDistribution::clone(void) const
{
    
	return new RootedTripletDistribution(*this);
}


/* Map calls to member methods */
RevLanguage::RevPtr<RevLanguage::RevVariable> RootedTripletDistribution::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
    
    if (name == "nTriplets")
    {
        found = true;
        
        size_t n = this->dag_node->getValue().getNumberOfTriplets();
        return new RevVariable( new Natural( n ) );
    }
    else if (name == "nTrees")
    {
        found = true;
        
        size_t n = this->dag_node->getValue().getNumberOfTrees();
        return new RevVariable( new Natural( n ) );
    }
    else if (name == "species")
    {
        found = true;
        
        const std::vector<std::string>& n = this->dag_node->getValue().getSpecies();
        return new RevVariable( new ModelVector<RlString>( n ) );
    }
    else if (name == "taxa")
    {
        found = true;
        
        const std::vector<RevBayesCore::Taxon>& n = this->dag_node->getValue().getTaxa();
        return new RevVariable( new ModelVector<Taxon>( n ) );
    }
    else if (name == "setTaxa")
    {
        found = true;
        
        std::vector<RevBayesCore::Taxon> t = static_cast<const ModelVector<Taxon>  &>( args[0].getVariable()->getRevObject() ).getValue();
        RevBayesCore::RootedTripletDistribution &dist = dag_node->getValue();
        dist.setTaxa(t);
        
        return NULL;
    }
    else if (name == "setSpecies")
    {
        found = true;
        
        std::vector<std::string> t = static_cast<const ModelVector<RlString>  &>( args[0].getVariable()->getRevObject() ).getValue();
        RevBayesCore::RootedTripletDistribution &dist = dag_node->getValue();
        dist.setSpecies(t);
        
        return NULL;
    }
/*    else if (name == "extractTriplets")
    {
        std::vector<RevBayesCore::Tree> t = static_cast<const ModelVector<TimeTree>  &>( args[0].getVariable()->getRevObject() ).getValue();
        RevBayesCore::RootedTripletDistribution &dist = dagNode->getValue();
        dist.extractTriplets(t);
        
        return NULL;
    }*/

    return ModelObject<RevBayesCore::RootedTripletDistribution>::executeMethod( name, args, found );
}


/** Get Rev type of object */
const std::string& RootedTripletDistribution::getClassType(void) {
    
    static std::string rev_type = "RootedTripletDistribution";
    
	return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& RootedTripletDistribution::getClassTypeSpec(void) {
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( RevObject::getClassTypeSpec() ) );
    
	return rev_type_spec;
}


/** Get type spec */
const TypeSpec& RootedTripletDistribution::getTypeSpec( void ) const {
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


