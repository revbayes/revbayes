// #include "ConstantNode.h"
// #include "Ellipsis.h"
// #include "ModelVector.h"
// #include "RlClade.h"
// #include "RlSubsplit.h"
// #include "RbUtil.h"
// #include "RealPos.h"
// #include "RlString.h"
// #include "RlTaxon.h"
// #include "TypeSpec.h"
//
// #include <sstream>
//
// using namespace RevLanguage;
//
// /** Default constructor */
// Subsplit::Subsplit(void) : ModelObject<RevBayesCore::Subsplit>()
// {
//
//     initMethods();
//
// }
//
// /** Construct from core Subsplit */
// Subsplit::Subsplit(RevBayesCore::Subsplit *c) : ModelObject<RevBayesCore::Subsplit>( c )
// {
//
//     initMethods();
//
// }
//
// /** Construct from core Subsplit */
// Subsplit::Subsplit(const RevBayesCore::Subsplit &t) : ModelObject<RevBayesCore::Subsplit>( new RevBayesCore::Subsplit( t ) )
// {
//
//     initMethods();
//
// }
//
// /** Construct from DAG node */
// Subsplit::Subsplit(RevBayesCore::TypedDagNode<RevBayesCore::Subsplit> *n) : ModelObject<RevBayesCore::Subsplit>( n )
// {
//
//     initMethods();
//
// }
//
//
// /**
//  * The clone function is a convenience function to create proper copies of inherited objected.
//  * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
//  *
//  * \return A new copy of the process.
//  */
// Subsplit* Subsplit::clone(void) const
// {
//
//   return new Subsplit(*this);
// }
//
//
// void Subsplit::constructInternalObject( void )
// {
//   // we free the memory first
//     if ( dag_node != NULL )
//     {
//         if ( dag_node->decrementReferenceCount() == 0 )
//         {
//             delete dag_node;
//         }
//     }
//
//     // Make taxa useful
//     // std::vector<RevBayesCore::Taxon> t = static_cast<const ModelVector<Taxon> &>( taxa->getRevObject());
//     const ModelVector<Taxon> &rb_taxa = static_cast<const ModelVector<Taxon> &>( taxa->getRevObject());
//
//     std::vector<RevBayesCore::Taxon> t;
//
//     for (size_t i = 0; i<static_cast<const ModelVector<Taxon> &>( taxa->getRevObject()).size(); ++i)
//     {
//         t.push_back(rb_taxa[i]);
//     }
//
//     // ModelVector<RlString> rlvector = ModelVector<RlString>();
//     // for (size_t i = 0; i<static_cast<const ModelVector<Taxon> &>( var->getRevObject()).size(); ++i)
//     // {
//     //     std::string name = static_cast<const ModelVector<Taxon> &>( var->getRevObject())[i].getSpeciesName();
//     //     const RlString rlName  = RlString (name);
//     //     rlvector.push_back(rlName);
//     // }
//     // const RevPtr<const RevVariable>  rvptr = new const RevVariable(new ModelVector<RlString> (rlvector) );
//     // taxa.push_back( rvptr );
//
//     // make sure subsplit has valid number of constituent clades
//     if ( clades.size() == 1 )
//     {
//       // Fake subsplit
//       const RevBayesCore::Clade &c = static_cast<const Clade &>( clades[0]->getRevObject() ).getValue();
//       RevBayesCore::Subsplit *s = new RevBayesCore::Subsplit(c, t);
//       dag_node = new RevBayesCore::ConstantNode<RevBayesCore::Subsplit>("",s);
//     }
//     else if ( clades.size() == 2 )
//     {
//       // Real subsplit, core functions will handle ordering
//       const RevBayesCore::Clade &c1 = static_cast<const Clade &>( clades[0]->getRevObject() ).getValue();
//       const RevBayesCore::Clade &c2 = static_cast<const Clade &>( clades[1]->getRevObject() ).getValue();
//       RevBayesCore::Subsplit *s = new RevBayesCore::Subsplit(c1, c2, t);
//       dag_node = new RevBayesCore::ConstantNode<RevBayesCore::Subsplit>("",s);
//     }
//     else
//     {
//       throw(RbException("Can only initialize subsplit object with 1 or 2 constituent clades."));
//     }
//
//     dag_node->incrementReferenceCount();
//
// }
//
// /* Map calls to member methods */
// RevLanguage::RevPtr<RevLanguage::RevVariable> Subsplit::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
// {
//
//     if (name == "size")
//     {
//         found = true;
//
//         return new RevVariable( new Natural( this->dag_node->getValue().size() ) );
//     }
//     else if (name == "getCladeY" )
//     {
//         found = true;
//
//         RevBayesCore::Clade c = this->dag_node->getValue().getY(  );
//         return new RevVariable( new Clade( c ) );
//     }
//     else if (name == "getCladeZ" )
//     {
//         found = true;
//
//         RevBayesCore::Clade c = this->dag_node->getValue().getZ(  );
//         return new RevVariable( new Clade( c ) );
//     }
//
//
//     return ModelObject<RevBayesCore::Subsplit>::executeMethod( name, args, found );
// }
//
//
//
// /**
//  * Get the Rev name for the constructor function.
//  *
//  * \return Rev name of constructor function.
//  */
// std::string Subsplit::getConstructorFunctionName( void ) const
// {
//     // create a constructor function name variable that is the same for all instance of this class
//     std::string s_name = "subsplit";
//
//     return s_name;
// }
//
//
// /**
//  * Get the author(s) of this function so they can receive credit (and blame) for it.
//  */
// std::vector<std::string> Subsplit::getHelpAuthor(void) const
// {
//     // create a vector of authors for this function
//     std::vector<std::string> authors;
//     authors.push_back( "Andrew Magee" );
//
//     return authors;
// }
//
//
// /**
//  * Get the (brief) description for this function
//  */
// std::string Subsplit::getHelpDescription(void) const
// {
//     // create a variable for the description of the function
//     std::string description = "A subsplits are pairs of clades, equivalently, a subdivision of a clade into its constituent clades.";
//
//     return description;
// }
//
//
// /**
//  * Get the more detailed description of the function
//  */
// std::string Subsplit::getHelpDetails(void) const
// {
//     // create a variable for the description of the function
//     std::string details = "";
//     details += "Subsplits are defined in terms of clades, and represent the splitting of a clade into ";
//     details += "two sister clades. Internally, these clades are sorted alphabetically. Subsplits are ";
//     details += "primarily for use with the SBN distributions on trees.";
//
//     return details;
// }
//
//
// /**
//  * Get an executable and instructive example.
//  * These example should help the users to show how this function works but
//  * are also used to test if this function still works.
//  */
// std::string Subsplit::getHelpExample(void) const
// {
//     // create an example as a single string variable.
//     std::string example = "";
//     example += "# create taxa\n";
//     example += "for (i in 1:4) {\n";
//     example += "taxa[i] = taxon(\"t\" + i)\n";
//     example += "}\n";
//     example += "# create two clades\n";
//     example += "clade_1 = clade( taxa[1], taxa[2] )\n";
//     example += "clade_2 = clade( taxa[3], taxa[4] )\n";
//     example += "# create a subsplit from the clades\n";
//     example += "subsplit_1 = subsplit( clade_1, clade_2 )\n";
//     example += "# note that we can reverse the ordering without changing the subsplit\n";
//     example += "# t1 comes first alphabetically, so the clade with t1 will be first\n";
//     example += "# we can also make fake subsplits\n";
//     example += "# these include a single clade that includes a single taxon\n";
//     example += "clade_3 = clade( 't4' )\n";
//     example += "fake_subsplit = subsplit( clade_3, taxa )\n";
//
//     return example;
// }
//
//
// /**
//  * Get some references/citations for this function
//  *
//  */
// std::vector<RevBayesCore::RbHelpReference> Subsplit::getHelpReferences(void) const
// {
//     // create an entry for each reference
//     std::vector<RevBayesCore::RbHelpReference> references;
//
//
//     return references;
// }
//
//
// /**
//  * Get the names of similar and suggested other functions
//  */
// std::vector<std::string> Subsplit::getHelpSeeAlso(void) const
// {
//     // create an entry for each suggested function
//     std::vector<std::string> see_also;
//     see_also.push_back( "clade" );
//     see_also.push_back( "dnUnconstrainedSBN" );
//
//     return see_also;
// }
//
//
// /**
//  * Get the title of this help entry
//  */
// std::string Subsplit::getHelpTitle(void) const
// {
//     // create a title variable
//     std::string title = "Subsplit";
//
//     return title;
// }
//
//
// /** Return member rules (no members) */
// const MemberRules& Subsplit::getParameterRules(void) const
// {
//
//     static MemberRules memberRules;
//     static bool rules_set = false;
//
//     if ( !rules_set )
//     {
//
//         memberRules.push_back( new Ellipsis( "Clades as clade objects.", Clade::getClassTypeSpec() ) );
//         memberRules.push_back( new ArgumentRule("taxa", ModelVector<Taxon>::getClassTypeSpec(), "Taxon names as a vector of taxons.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
//
//         rules_set = true;
//     }
//
//     return memberRules;
// }
//
//
// /** Get Rev type of object */
// const std::string& Subsplit::getClassType(void)
// {
//
//   static std::string rev_type = "Subsplit";
//
//   return rev_type;
// }
//
// /** Get class type spec describing type of object */
// const TypeSpec& Subsplit::getClassTypeSpec(void)
// {
//
//     static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( ModelObject<RevBayesCore::Subsplit>::getClassTypeSpec() ) );
//
//     return rev_type_spec;
// }
//
//
// /** Get type spec */
// const TypeSpec& Subsplit::getTypeSpec( void ) const
// {
//
//     static TypeSpec type_spec = getClassTypeSpec();
//
//     return type_spec;
// }
//
//
// /**
//  * Initialize the member methods.
//  */
// void Subsplit::initMethods( void )
// {
//
//
//     ArgumentRules* sizeRules = new ArgumentRules();
//     methods.addFunction( new MemberProcedure( "size", RlUtils::Void, sizeRules ) );
//
//     ArgumentRules* cladeYRules = new ArgumentRules();
//     methods.addFunction( new MemberProcedure( "getCladeY", RlUtils::Void, cladeYRules ) );
//
//     ArgumentRules* cladeZRules = new ArgumentRules();
//     methods.addFunction( new MemberProcedure( "getCladeZ", RlUtils::Void, cladeZRules ) );
//
//
// }
//
//
// /** Set a member variable */
// void Subsplit::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
// {
//     // Subsplits are initialized from clade objects
//     if ( name == "" && var->getRevObject().getTypeSpec() == Clade::getClassTypeSpec() )
//     {
//         clades.push_back( var );
//     }
//     else if ( name == "taxa" )
//     {
//
//       taxa = var;
//     }
//     else
//     {
//         RevObject::setConstParameter(name, var);
//     }
// }
