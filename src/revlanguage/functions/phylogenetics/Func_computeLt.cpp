#include "Func_computeLt.h"
#include "ComputeLtFunction.h"

// you can probably get rid of some of these
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ModelVector.h"
#include "OptionRule.h"
#include "RevObject.h"
#include "RlTree.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"

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

  //RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* a = static_cast<const ModelVector<Real> &>( this->args[0].getVariable()->getRevObject() ).getDagNode();

  RevBayesCore::DagNode* node = this->args[0].getVariable()->getRevObject().getDagNode();
  // vector of branching times or tree
  RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* a = dynamic_cast<RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >*>(node);
  RevBayesCore::TypedDagNode< RevBayesCore::Tree>* alt = dynamic_cast< RevBayesCore::TypedDagNode< RevBayesCore::Tree >*>(node);

  RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* b = static_cast<const ModelVector<Real> &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
  RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* c = static_cast<const ModelVector<Real> &>( this->args[2].getVariable()->getRevObject() ).getDagNode();
  RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* d = static_cast<const ModelVector<Real> &>( this->args[3].getVariable()->getRevObject() ).getDagNode();
  RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* e = static_cast<const ModelVector<Real> &>( this->args[4].getVariable()->getRevObject() ).getDagNode();
  RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* f = static_cast<const ModelVector<Real> &>( this->args[5].getVariable()->getRevObject() ).getDagNode();

  RevBayesCore::TypedDagNode< double >* t = static_cast<const RealPos &>( this->args[6].getVariable()->getRevObject() ).getDagNode();
  RevBayesCore::TypedDagNode< double >* l = static_cast<const RealPos &>( this->args[7].getVariable()->getRevObject() ).getDagNode();
  RevBayesCore::TypedDagNode< double >* m = static_cast<const RealPos &>( this->args[8].getVariable()->getRevObject() ).getDagNode();
  RevBayesCore::TypedDagNode< double >* p = static_cast<const RealPos &>( this->args[9].getVariable()->getRevObject() ).getDagNode();
  RevBayesCore::TypedDagNode< double >* o = static_cast<const RealPos &>( this->args[10].getVariable()->getRevObject() ).getDagNode();
  RevBayesCore::TypedDagNode< double >* rh = static_cast<const RealPos &>( this->args[11].getVariable()->getRevObject() ).getDagNode();
  RevBayesCore::TypedDagNode< double >* rm = static_cast<const RealPos &>( this->args[12].getVariable()->getRevObject() ).getDagNode();

  //RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* g = static_cast<const ModelVector<Real> &>( this->args[13].getVariable()->getRevObject() ).getDagNode();
  // times of interest
  RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* g = NULL;
  if ( listG->getRevObject() != RevNullObject::getInstance() )
  {
      g = static_cast<const ModelVector<Real> &>( listG->getRevObject() ).getDagNode();
  }

  //RevBayesCore::ComputeLtFunction* fxn = new RevBayesCore::ComputeLtFunction( a, b, c, d, e, f, t, l, m, p, o, rho, r, g );
  RevBayesCore::ComputeLtFunction* fxn = a != NULL ? new RevBayesCore::ComputeLtFunction( a, b, c, d, e, f, t, l, m, p, o, rh, rm, g ) : new RevBayesCore::ComputeLtFunction( alt, b, c, d, e, f, t, l, m, p, o, rh, rm, g );

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

        argumentRules.push_back( new ArgumentRule( "origin", RealPos::getClassTypeSpec(), "Origin time.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "lambda", RealPos::getClassTypeSpec(), "lambda.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "mu", RealPos::getClassTypeSpec(), "mu.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "psi", RealPos::getClassTypeSpec(), "psi.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "omega", RealPos::getClassTypeSpec(), "omega.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "rho", RealPos::getClassTypeSpec(), "rho.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "removal", RealPos::getClassTypeSpec(), "removel Pr.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        argumentRules.push_back( new ArgumentRule( "listG", ModelVector<Real>::getClassTypeSpec(), "A vector of values.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

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

/**
 * Set a member variable.
 *
 * Sets a member variable with the given name and store the pointer to the variable.
 * The value of the variable might still change but this function needs to be called again if the pointer to
 * the variable changes. The current values will be used to create the distribution object.
 *
 * \param[in]    name     Name of the member variable.
 * \param[in]    var      Pointer to the variable.
 */
void Func_computeLt::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    if ( name == "listG" )
    {
        listG = var;
    }
    else {
      Function::setConstParameter(name, var);
    }
}
