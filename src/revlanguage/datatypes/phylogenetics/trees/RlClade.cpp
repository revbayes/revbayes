#include "ConstantNode.h"
#include "Ellipsis.h"
#include "ModelVector.h"
#include "RlClade.h"
#include "RbUtil.h"
#include "RealPos.h"
#include "RlString.h"
#include "RlTaxon.h"
#include "TypeSpec.h"

#include <sstream>

using namespace RevLanguage;

/** Default constructor */
Clade::Clade(void) : ModelObject<RevBayesCore::Clade>()
{
    
    initMethods();

}

/** Construct from core Clade */
Clade::Clade(RevBayesCore::Clade *c) : ModelObject<RevBayesCore::Clade>( c )
{
    
    initMethods();

}

/** Construct from core Clade */
Clade::Clade(const RevBayesCore::Clade &t) : ModelObject<RevBayesCore::Clade>( new RevBayesCore::Clade( t ) )
{
    
    initMethods();

}

/** Construct from DAG node */
Clade::Clade(RevBayesCore::TypedDagNode<RevBayesCore::Clade> *n) : ModelObject<RevBayesCore::Clade>( n )
{
    
    initMethods();

}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Clade* Clade::clone(void) const
{

  return new Clade(*this);
}


void Clade::constructInternalObject( void )
{
  // we free the memory first
    if ( dag_node != NULL )
    {
        if ( dag_node->decrementReferenceCount() == 0 )
        {
            delete dag_node;
        }
    }

    // get clade constraint flags
    bool match = static_cast<const RlBoolean &>( optional_match->getRevObject() ).getValue();
    bool neg = static_cast<const RlBoolean &>( is_negative_constraint->getRevObject() ).getValue();

    // now allocate a new Clade
    std::set<RevBayesCore::Taxon> n_set;
    for (std::vector<RevPtr<const RevVariable> >::iterator it = names.begin(); it != names.end(); ++it)
    {
        RevBayesCore::Taxon t = RevBayesCore::Taxon(static_cast<const RlString &>( (*it)->getRevObject() ).getValue());
        n_set.insert( t );
    }

    // now add names to new clade
    for (std::vector<RevPtr<const RevVariable> >::iterator it = names_vector.begin(); it != names_vector.end(); ++it)
    {
        const ModelVector<RlString> &tmp = static_cast<const ModelVector<RlString> &>( (*it)->getRevObject() );

        for (size_t i=0; i<tmp.size(); ++i)
        {
            RevBayesCore::Taxon t = RevBayesCore::Taxon( tmp.getElement(i)->getValue() );
            n_set.insert( t );
        }
    }

    // now allocate a new Clade
    std::vector<RevBayesCore::Clade> optional_constraints;
    for (std::vector<RevPtr<const RevVariable> >::iterator it = clades.begin(); it != clades.end(); ++it)
    {

        const RevBayesCore::Clade &c = static_cast<const Clade &>( (*it)->getRevObject() ).getValue();

        if (match)
        {
            optional_constraints.push_back(c);
        }

        const std::vector<RevBayesCore::Taxon> &taxa = c.getTaxa();
        for (size_t i=0; i<taxa.size(); ++i)
        {
            const RevBayesCore::Taxon &t = taxa[i];
            n_set.insert( t );
        }
    }
    
    
    // convert to vector
    std::vector<RevBayesCore::Taxon> n;
    for (std::set<RevBayesCore::Taxon>::iterator it = n_set.begin(); it != n_set.end(); it++)
    {
        n.push_back( *it );
    }


    RevBayesCore::Clade *c = new RevBayesCore::Clade(n);

    // set the age if provided
    if ( age->getRevObject() != RevNullObject::getInstance() )
    {
        double a = static_cast<const RealPos &>( age->getRevObject() ).getValue();
        c->setAge( a );
    }

    // set the number of missing if provided
    if ( missing->getRevObject() != RevNullObject::getInstance() )
    {
        long n = static_cast<const Natural &>( missing->getRevObject() ).getValue();
        c->setNumberMissingTaxa( (int)n );
    }

    // set optional clade constraints if provided
    if (match && optional_constraints.size() > 0)
    {
        c->setOptionalConstraints( optional_constraints );
    }


    // set optional match clade constraint
    c->setOptionalMatch( match );

    // set negative clade constraint
    c->setNegativeConstraint( neg );


    dag_node = new RevBayesCore::ConstantNode<RevBayesCore::Clade>("", c);
    dag_node->incrementReferenceCount();

}

/* Map calls to member methods */
RevLanguage::RevPtr<RevLanguage::RevVariable> Clade::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
    
    if (name == "size")
    {
        found = true;
        
        return new RevVariable( new Natural( this->dag_node->getValue().size() ) );
    }
    else if (name == "getTaxonName" )
    {
        found = true;
        int index = (int)static_cast<const Natural&>( args[0].getVariable()->getRevObject() ).getValue() - 1;
        
        std::string t = this->dag_node->getValue().getTaxonName( index );
        return new RevVariable( new RlString( t ) );
    }
    else if (name == "getTaxon" )
    {
        found = true;
        int index = (int)static_cast<const Natural&>( args[0].getVariable()->getRevObject() ).getValue() - 1;
        
        RevBayesCore::Taxon t = this->dag_node->getValue().getTaxon( index );
        return new RevVariable( new Taxon( t ) );
    }
    else if (name == "getNumberOfTaxaMissing" )
    {
        found = true;
        
        int n = this->dag_node->getValue().getNumberMissingTaxa();
        return new RevVariable( new Natural( n ) );
    }
    else if (name == "setNumberOfTaxaMissing" )
    {
        found = true;
        int n = (int)static_cast<const Natural&>( args[0].getVariable()->getRevObject() ).getValue();
        
        this->dag_node->getValue().setNumberMissingTaxa( n );
        return NULL;
    }
    
    
    return ModelObject<RevBayesCore::Clade>::executeMethod( name, args, found );
}



/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Clade::getConstructorFunctionName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "clade";

    return c_name;
}


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Clade::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Will Pett" );
    authors.push_back( "Sebastian Hoehna" );
    authors.push_back( "Michael Landis" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Clade::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::string description = "A clade is a subtree within a phylogeny.";
    
    return description;
}


/**
 * Get the more detailed description of the function
 */
std::string Clade::getHelpDetails(void) const
{
    // create a variable for the description of the function
    std::string details = "";
    details += "Clades are defined in terms of a taxon set and a shared tree topology. In phylogenetic ";
    details += "analyses, clades are generally used (a) to constrain tree topologies to match provided ";
    details += "taxon relationships, (b) to identify the most recent common ancestor of a taxon set ";
    details += "within a phylogeny, or (c) to apply node age calibrations on particular nodes in the phylogeny.";
    
    return details;
}


/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
std::string Clade::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";
    
    example += "# read in a tree\n";
    example += "phy = readTrees(\"primates.tre\")[1]\n";
    example += "# get taxa from the tree\n";
    example += "taxa = phy.taxa()\n";
    example += "# create a clade for (1,2) using taxon objects\n";
    example += "clade_12 = clade( taxa[1], taxa[2] )\n";
    example += "# create a clade for (1,2,3)\n";
    example += "clade_123 = clade( taxa[3], clade_12 )\n";
    example += "# create a clade for (4,5) using taxon names\n";
    example += "clade_45 = clade( \"Taxon_4\", \"Taxon_5\" )\n";
    example += "# create a negative clade constraint\n";
    example += "clade_not_23 = clade( taxa[2], taxa[3], negative=true )\n";
    example += "# create an optional clade constraint\n";
    example += "clade_67 = clade( taxa[6], taxa[7] )\n";
    example += "clade_68 = clade( taxa[6], taxa[8] )\n";
    example += "clade_67_or_68 = clade( clade_67, clade_68, optional_match=true )\n";
    
    return example;
}


/**
 * Get some references/citations for this function
 *
 */
std::vector<RevBayesCore::RbHelpReference> Clade::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    
    
    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Clade::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    see_also.push_back( "dnConstrainedTopology" );
    see_also.push_back( "tmrca" );
    see_also.push_back( "mrcaIndex" );
    
    return see_also;
}


/**
 * Get the title of this help entry
 */
std::string Clade::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "Clade";
    
    return title;
}


/** Return member rules (no members) */
const MemberRules& Clade::getParameterRules(void) const
{

    static MemberRules memberRules;
    static bool rules_set = false;

    if ( !rules_set )
    {

        memberRules.push_back( new Ellipsis( "Taxon names as string values.", RlString::getClassTypeSpec() ) );
        memberRules.push_back( new Ellipsis("Taxon names as a vector of string values.", ModelVector<RlString>::getClassTypeSpec() ) );
        memberRules.push_back( new Ellipsis( "Taxa as clade objects.", Clade::getClassTypeSpec() ) );
        memberRules.push_back( new Ellipsis( "Taxon names as taxon values", Taxon::getClassTypeSpec() ) );
        memberRules.push_back( new Ellipsis( "Taxon names as a vector of taxons", ModelVector<Taxon>::getClassTypeSpec() ) );
        memberRules.push_back( new ArgumentRule("age", RealPos::getClassTypeSpec(), "The age of the clade (optional).", ArgumentRule::BY_VALUE, ArgumentRule::ANY, NULL ) );
        memberRules.push_back( new ArgumentRule("missing", Natural::getClassTypeSpec(), "Number of missing taxa in the clade (optional).", ArgumentRule::BY_VALUE, ArgumentRule::ANY, NULL ) );
        memberRules.push_back( new ArgumentRule("negative", RlBoolean::getClassTypeSpec(), "Is this a negative clade constraint?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ) );
        memberRules.push_back( new ArgumentRule("optional_match", RlBoolean::getClassTypeSpec(), "Clade constraint satisfied when any Clade argument matched", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ) );

        rules_set = true;
    }

    return memberRules;
}


/** Get Rev type of object */
const std::string& Clade::getClassType(void)
{

  static std::string rev_type = "Clade";

  return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& Clade::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( ModelObject<RevBayesCore::Clade>::getClassTypeSpec() ) );

    return rev_type_spec;
}


/** Get type spec */
const TypeSpec& Clade::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}


/**
 * Initialize the member methods.
 */
void Clade::initMethods( void )
{
    
    
    ArgumentRules* ntipsArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "size", RlUtils::Void, ntipsArgRules ) );
    
    ArgumentRules* namesArgRules = new ArgumentRules();
    namesArgRules->push_back( new ArgumentRule( "node", Natural::getClassTypeSpec(), "The index of the node.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
    methods.addFunction( new MemberProcedure( "getTaxonName", RlString::getClassTypeSpec(), namesArgRules ) );
    
    ArgumentRules* taxaArgRules = new ArgumentRules();
    taxaArgRules->push_back( new ArgumentRule( "node", Natural::getClassTypeSpec(), "The index of the node.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
    methods.addFunction( new MemberProcedure( "getTaxon", Taxon::getClassTypeSpec(), taxaArgRules ) );

    ArgumentRules* get_num_missing_arg_rules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "getNumberOfTaxaMissing", Natural::getClassTypeSpec(), get_num_missing_arg_rules ) );

    ArgumentRules* set_num_missing_arg_rules = new ArgumentRules();
    set_num_missing_arg_rules->push_back( new ArgumentRule( "n", Natural::getClassTypeSpec(), "The number of missing taxa.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
    methods.addFunction( new MemberProcedure( "setNumberOfTaxaMissing", RlUtils::Void, set_num_missing_arg_rules ) );

    
}


/** Set a member variable */
void Clade::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{

    if ( name == "taxonName" || (name == "" && var->getRevObject().getTypeSpec() == RlString::getClassTypeSpec() ) )
    {
        names.push_back( var );
    }
    else if ( name == "" && var->getRevObject().getTypeSpec() == Clade::getClassTypeSpec() )
    {
        clades.push_back( var );
    }
    else if ( name == "" && var->getRevObject().getTypeSpec() == ModelVector<RlString>::getClassTypeSpec() )
    {
        names_vector.push_back( var );
    }
    //Alternatively, one could provide taxons, not names
    else if ( name == "taxonName" || (name == "" && var->getRevObject().getTypeSpec() == Taxon::getClassTypeSpec() ) )
    {
        std::string name = static_cast<const Taxon &>( var->getRevObject() ).getDagNode()->getValue().getName();
        const RevPtr<const RevVariable> rlName = new const RevVariable(new RlString( name )) ;
        names.push_back( rlName );
    }
    else if ( name == "" && var->getRevObject().getTypeSpec() == ModelVector<Taxon>::getClassTypeSpec() )
    {
        ModelVector<RlString> rlvector = ModelVector<RlString>();
        for (size_t i = 0; i<static_cast<const ModelVector<Taxon> &>( var->getRevObject()).size(); ++i)
        {
            std::string name = static_cast<const ModelVector<Taxon> &>( var->getRevObject())[i].getSpeciesName();
            //const RevPtr<const RevVariable> rlName = new const RevVariable(new RlString( name )) ;
            const RlString rlName  = RlString (name);
            rlvector.push_back(rlName);
        }
        const RevPtr<const RevVariable>  rvptr = new const RevVariable(new ModelVector<RlString> (rlvector) );
        names_vector.push_back( rvptr );
    }
    else if ( name == "age")
    {
        age = var;
    }
    else if ( name == "missing")
    {
        missing = var;
    }
    else if ( name == "negative" )
    {
        is_negative_constraint = var;
    }
    else if ( name == "optional_match" )
    {
        optional_match = var;
    }
  //    else if ( name == "optional_constraints" )
  //    {
  //        optional_constraints = var;
  //    }
    else
    {
        RevObject::setConstParameter(name, var);
    }
}
