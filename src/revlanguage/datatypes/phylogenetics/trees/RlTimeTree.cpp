#include <cstddef>
#include <sstream>
#include <string>
#include <vector>

#include "ModelVector.h"
#include "Natural.h"
#include "RlBoolean.h"
#include "RlBranchLengthTree.h"
#include "RlTimeTree.h"
#include "RlMemberFunction.h"
#include "RlTaxon.h"
#include "RealPos.h"
#include "TopologyNode.h"
#include "TypeSpec.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ConstantNode.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "MemberFunction.h"
#include "MemberProcedure.h"
#include "MethodTable.h"
#include "ModelObject.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "Real.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlConstantNode.h"
#include "RlDeterministicNode.h"
#include "RlTree.h"
#include "RlTypedFunction.h"
#include "RlUtils.h"
#include "Taxon.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

using namespace RevLanguage;

/** Default constructor */
TimeTree::TimeTree(void) : Tree()
{

    // initialize the member methods
    initMethods();
    
}

/** Construct from core object */
TimeTree::TimeTree(RevBayesCore::Tree *t) : Tree( t )
{
    
    // initialize the member methods
    initMethods();

}

/** Construct from bool */
TimeTree::TimeTree(const RevBayesCore::Tree &t) : Tree( t )
{
    
    
    // initialize the member methods
    initMethods();

}

/** Construct from bool */
TimeTree::TimeTree(RevBayesCore::TypedDagNode<RevBayesCore::Tree> *n) : Tree( n )
{
    
    
    // initialize the member methods
    initMethods();

}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
TimeTree* TimeTree::clone(void) const
{
    
    return new TimeTree(*this);
}


/* Map calls to member methods */
RevLanguage::RevPtr<RevLanguage::RevVariable> TimeTree::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
    
    if (name == "collapseNegativeBranches")
    {
        found = true;
        
        double length = static_cast<const RealPos&>( args[0].getVariable()->getRevObject() ).getValue();
        this->dag_node->getValue().collapseNegativeBranchLengths(length);
        return NULL;
    }
    else if (name == "dropFossils")
    {
        found = true;

        std::vector<RevBayesCore::Taxon> t = this->dag_node->getValue().getFossilTaxa();
        for (size_t i = 0; i < t.size(); i++)
        {
            std::string taxon_name = t[i].getName();
            this->dag_node->getValue().dropTipNodeWithName( taxon_name );
        }
        return NULL;
    }
    else if (name == "getFossils")
    {
        found = true;
        
        std::vector<RevBayesCore::Taxon> t = this->dag_node->getValue().getFossilTaxa();
        return new RevVariable( new ModelVector<Taxon>( t ) );
    }
    else if (name == "isRoot")
    {
        found = true;
        
        int index = (int)static_cast<const Natural&>( args[0].getVariable()->getRevObject() ).getValue() - 1;
        
        bool tf = this->dag_node->getValue().getNode((size_t)index).isRoot();
        return new RevVariable( new RlBoolean( tf ) );
    }
    else if (name == "numSampledAncestors")
    {
        found = true;

        size_t n = this->dag_node->getValue().getNumberOfTips();

        size_t num = 0;
        for (size_t i=0; i<n; i++){
            RevBayesCore::TopologyNode &node = this->dag_node->getValue().getNode(i);
            num += node.isSampledAncestorTip();
        }
        return new RevVariable( new Natural( num ) );
    }
    else if (name == "setAge")
    {
        found = true;

        const RevObject& current = args[0].getVariable()->getRevObject();
        if ( current.isType( Natural::getClassTypeSpec() ) )
        {
            size_t index = static_cast<const Natural&>( args[0].getVariable()->getRevObject() ).getValue() - 1;
            const RevObject& new_value = args[1].getVariable()->getRevObject();
            if ( new_value.isType( RealPos::getClassTypeSpec() ) )
            {
                double value = static_cast<const RealPos&>( new_value ).getValue();
                RevBayesCore::Tree &tree = dag_node->getValue();
                tree.getNode(index).setAge(value);
            }
        }
        return NULL;
    }
    else if (name == "unroot")
    {
        found = true;
        
        RevBayesCore::Tree *unrooted = this->dag_node->getValue().clone();
        
        // now unroot the tree
        unrooted->unroot();
        
        return new RevVariable( new BranchLengthTree( unrooted ) );
    }
    
    return Tree::executeMethod( name, args, found );
}


/** Get Rev type of object */
const std::string& TimeTree::getClassType(void)
{
    
    static std::string rev_type = "TimeTree";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& TimeTree::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Tree::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/** Get type spec */
const TypeSpec& TimeTree::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


void TimeTree::initMethods( void )
{
    
    ArgumentRules* is_root_arg_rules = new ArgumentRules();
    is_root_arg_rules->push_back( new ArgumentRule( "node", Natural::getClassTypeSpec(), "The index of the node.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    methods.addFunction( new MemberProcedure( "isRoot", RlBoolean::getClassTypeSpec(), is_root_arg_rules ) );

    ArgumentRules* drop_fossils_arg_rules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "dropFossils", RlUtils::Void, drop_fossils_arg_rules ) );

    ArgumentRules* get_fossils_arg_rules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "getFossils", ModelVector<Taxon>::getClassTypeSpec(), get_fossils_arg_rules ) );

    ArgumentRules* collapse_negative_branches_arg_rules = new ArgumentRules();
    collapse_negative_branches_arg_rules->push_back( new ArgumentRule( "length", RealPos::getClassTypeSpec(), "The new length of all negative branches.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(0.0) ));
    methods.addFunction( new MemberProcedure( "collapseNegativeBranches", RlUtils::Void, collapse_negative_branches_arg_rules ) );

    ArgumentRules* n_sampled_ancestors_arg_rules = new ArgumentRules();
    methods.addFunction( new MemberFunction<TimeTree, Natural>( "numSampledAncestors", this, n_sampled_ancestors_arg_rules ) );

    // member functions
    ArgumentRules* height_arg_rules = new ArgumentRules();
    methods.addFunction( new MemberFunction<TimeTree, RealPos>( "rootAge", this, height_arg_rules   ) );
    
    ArgumentRules* node_age_arg_rules = new ArgumentRules();
    node_age_arg_rules->push_back( new ArgumentRule( "node", Natural::getClassTypeSpec(), "The index of the node.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
    methods.addFunction( new MemberFunction<TimeTree, RealPos>( "nodeAge", this, node_age_arg_rules   ) );

    ArgumentRules* colless_arg_rules = new ArgumentRules();
    methods.addFunction( new MemberFunction<TimeTree, Natural>( "colless", this, colless_arg_rules ) );
    
    ArgumentRules* gamma_arg_rules = new ArgumentRules();
    methods.addFunction( new MemberFunction<TimeTree, Real>( "gammaStatistic", this, gamma_arg_rules ) );

    ArgumentRules* set_age_arg_rules         = new ArgumentRules();
    set_age_arg_rules->push_back( new ArgumentRule("index"    , Natural::getClassTypeSpec(), "The index of the node.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    set_age_arg_rules->push_back( new ArgumentRule("value"    , RealPos::getClassTypeSpec(), "The new age value.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    methods.addFunction( new MemberProcedure( "setAge", RlUtils::Void, set_age_arg_rules ) );

    ArgumentRules* unroot_arg_rules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "unroot", BranchLengthTree::getClassTypeSpec(), unroot_arg_rules ) );

}
