#include "Func_computeLt.h"
#include "ComputeLtFunction.h"

#include "ModelVector.h"
#include "RlTree.h"
#include "TypeSpec.h"

namespace RevBayesCore { template <class valueType> class RbVector; }
namespace RevBayesCore { class Tree; }

using namespace RevLanguage;

/** default constructor */
Func_computeLt::Func_computeLt (void) : TypedFunction<RealPos>( )
{

}

/** clone function */
Func_computeLt* Func_computeLt::clone( void ) const
{

    return new Func_computeLt( *this );
}

/** you probably have to change the first line of this function */
RevBayesCore::TypedFunction<double>* Func_computeLt::createFunction( void ) const
{

  RevBayesCore::DagNode* node = this->args[0].getVariable()->getRevObject().getDagNode();

  // vector of branching times or tree
  //RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* a = static_cast<const ModelVector<Real> &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
  RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* a = dynamic_cast<RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >*>(node);
  RevBayesCore::TypedDagNode< RevBayesCore::Tree>* alt = dynamic_cast< RevBayesCore::TypedDagNode< RevBayesCore::Tree >*>(node);

  RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* b = static_cast<const ModelVector<Real> &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
  RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* c = static_cast<const ModelVector<Real> &>( this->args[2].getVariable()->getRevObject() ).getDagNode();
  RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* d = static_cast<const ModelVector<Real> &>( this->args[3].getVariable()->getRevObject() ).getDagNode();
  RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* e = static_cast<const ModelVector<Real> &>( this->args[4].getVariable()->getRevObject() ).getDagNode();
  RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* f = static_cast<const ModelVector<Real> &>( this->args[5].getVariable()->getRevObject() ).getDagNode();

  //RevBayesCore::TypedDagNode< double >* mu = static_cast<const RealPos &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
  RevBayesCore::TypedDagNode< double >* t = static_cast<const RealPos &>( this->args[6].getVariable()->getRevObject() ).getDagNode();
  RevBayesCore::TypedDagNode< double >* l = static_cast<const RealPos &>( this->args[7].getVariable()->getRevObject() ).getDagNode();
  RevBayesCore::TypedDagNode< double >* m = static_cast<const RealPos &>( this->args[8].getVariable()->getRevObject() ).getDagNode();
  RevBayesCore::TypedDagNode< double >* p = static_cast<const RealPos &>( this->args[9].getVariable()->getRevObject() ).getDagNode();
  RevBayesCore::TypedDagNode< double >* o = static_cast<const RealPos &>( this->args[10].getVariable()->getRevObject() ).getDagNode();
  RevBayesCore::TypedDagNode< double >* rho = static_cast<const RealPos &>( this->args[11].getVariable()->getRevObject() ).getDagNode();
  RevBayesCore::TypedDagNode< double >* r = static_cast<const RealPos &>( this->args[12].getVariable()->getRevObject() ).getDagNode();

  RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* g = static_cast<const ModelVector<Real> &>( this->args[13].getVariable()->getRevObject() ).getDagNode();

  //RevBayesCore::ComputeLtFunction* fxn = new RevBayesCore::ComputeLtFunction( a, b, c, d, e, f, t, l, m, p, o, rho, r, g );
  RevBayesCore::ComputeLtFunction* fxn = a != NULL ? new RevBayesCore::ComputeLtFunction( a, b, c, d, e, f, t, l, m, p, o, rho, r, g ) : new RevBayesCore::ComputeLtFunction( alt, b, c, d, e, f, t, l, m, p, o, rho, r, g );

  return fxn;

}


/* Get argument rules */
/* you probably have to edit the argument rules here */
/* but after this I don't think you need to change anything else in this file */
const ArgumentRules& Func_computeLt::getArgumentRules( void ) const
{

    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {

        std::vector<TypeSpec> paramTypes;
        paramTypes.push_back( ModelVector<Real>::getClassTypeSpec() );
        paramTypes.push_back( Tree::getClassTypeSpec() );

        argumentRules.push_back( new ArgumentRule( "listA", paramTypes, "A vector of values/tree.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        argumentRules.push_back( new ArgumentRule( "listB", ModelVector<Real>::getClassTypeSpec(), "A vector of values.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "listC", ModelVector<Real>::getClassTypeSpec(), "A vector of values.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "listD", ModelVector<Real>::getClassTypeSpec(), "A vector of values.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "listE", ModelVector<Real>::getClassTypeSpec(), "A vector of values.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "listF", ModelVector<Real>::getClassTypeSpec(), "A vector of values.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        argumentRules.push_back( new ArgumentRule( "tor", RealPos::getClassTypeSpec(), "Origin time.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "l", RealPos::getClassTypeSpec(), "lambda.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "m", RealPos::getClassTypeSpec(), "mu.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "p", RealPos::getClassTypeSpec(), "psi.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "o", RealPos::getClassTypeSpec(), "omega.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "rho", RealPos::getClassTypeSpec(), "rho.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "r", RealPos::getClassTypeSpec(), "removel Pr.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        argumentRules.push_back( new ArgumentRule( "listF", ModelVector<Real>::getClassTypeSpec(), "A vector of values.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }

    return argumentRules;
}

const std::string& Func_computeLt::getClassType(void)
{

    static std::string rev_type = "Func_computeLt";

    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_computeLt::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_computeLt::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "computeLt";

    return f_name;
}

const TypeSpec& Func_computeLt::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}
