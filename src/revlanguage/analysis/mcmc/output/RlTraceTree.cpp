#include <math.h>
#include <stddef.h>
#include <iostream>
#include <string>
#include <vector>

#include "ArgumentRules.h"
#include "MemberProcedure.h"
#include "MethodTable.h"
#include "ModelVector.h"
#include "Natural.h"
#include "Probability.h"
#include "RlBoolean.h"
#include "RlBranchLengthTree.h"
#include "RlClade.h"
#include "RlSBNParameters.h"
#include "RlString.h"
#include "RlTimeTree.h"
#include "RlTraceTree.h"
#include "RlTree.h"
#include "RlUtils.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "Clade.h"
#include "ConstantNode.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "Integer.h"
#include "ModelObject.h"
#include "RbBoolean.h"
#include "RbException.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlConstantNode.h"
#include "TraceTree.h"
#include "Tree.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"
#include "WorkspaceToCoreWrapperObject.h"


TraceTree::TraceTree(void) : WorkspaceToCoreWrapperObject<RevBayesCore::TraceTree>()
{

    // initialize the methods
    initMethods();

}


TraceTree::TraceTree(const RevBayesCore::TraceTree& x) : WorkspaceToCoreWrapperObject<RevBayesCore::TraceTree>( new RevBayesCore::TraceTree(x) )
{

    // initialize the methods
    initMethods();

}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */

TraceTree* TraceTree::clone(void) const
{

    return new TraceTree(*this);
}



void TraceTree::constructInternalObject( void )
{
    throw RbException("We do not support a constructor function for TraceTree.");
}


/* Map calls to member methods */

RevPtr<RevVariable> TraceTree::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{


    if ( name == "setBurnin" )
    {
        found = true;

        int burnin = 0;

        RevObject& b = args[0].getVariable()->getRevObject();
        if ( b.isType( Integer::getClassTypeSpec() ) )
        {
            burnin = (int)static_cast<const Integer &>(b).getValue();
        }
        else
        {
            double burninFrac = static_cast<const Probability &>(b).getValue();
            burnin = int( floor( value->size()*burninFrac ) );
        }

        this->value->setBurnin( burnin );

        return NULL;
    }
    else if ( name == "setOutgroup" )
    {
        found = true;

        const RevBayesCore::Clade &c    = static_cast<const Clade &>( args[0].getVariable()->getRevObject() ).getValue();

        this->value->setOutgroup(c);

        return NULL;

    }
    else if ( name == "summarize" )
    {
        found = true;

        double treeCI       = static_cast<const Probability &>( args[0].getVariable()->getRevObject() ).getValue();
        double minCladeProb = static_cast<const Probability &>( args[1].getVariable()->getRevObject() ).getValue();
        bool verbose = static_cast<const RlBoolean &>( args[2].getVariable()->getRevObject() ).getValue();

        this->value->printTreeSummary(std::cout, treeCI, verbose);
        this->value->printCladeSummary(std::cout, minCladeProb, verbose);

        return NULL;
    }
    else if ( name == "cladeProbability" )
    {
        found = true;

        const RevBayesCore::Clade &c    = static_cast<const Clade &>( args[0].getVariable()->getRevObject() ).getValue();
        bool verbose = static_cast<const RlBoolean &>( args[1].getVariable()->getRevObject() ).getValue();

        double p = this->value->cladeProbability( c, verbose );

        return new RevVariable( new Probability( p ) );

    }
    else if ( name == "computeEntropy" )
    {
        found = true;

        double tree_CI  = static_cast<const Probability &>( args[0].getVariable()->getRevObject() ).getValue();
        int num_taxa    = (int)static_cast<const Integer &>( args[1].getVariable()->getRevObject() ).getValue();
        bool verbose    = static_cast<const RlBoolean &>( args[2].getVariable()->getRevObject() ).getValue();

        double entropy = this->value->computeEntropy(tree_CI, num_taxa, verbose);

        return new RevVariable( new RealPos(entropy) );
    }
    else if ( name == "computePairwiseRFDistances" )
    {
        found = true;

        double tree_CI         = static_cast<const Probability &>( args[0].getVariable()->getRevObject() ).getValue();
        bool verbose           = static_cast<const RlBoolean &>( args[1].getVariable()->getRevObject() ).getValue();

        std::vector<double> distances = this->value->computePairwiseRFDistance(tree_CI, verbose);

        ModelVector<RealPos> *rl_dist = new ModelVector<RealPos>;
        for (size_t i=0; i<distances.size(); ++i)
        {
            rl_dist->push_back( distances[i] );
        }

        return new RevVariable( rl_dist );
    }
    else if ( name == "computeTreeLengths" )
    {
        found = true;

        std::vector<double> tree_lengths = this->value->computeTreeLengths();

        ModelVector<RealPos> *rl_tree_lengths = new ModelVector<RealPos>;
        for (size_t i=0; i<tree_lengths.size(); ++i)
        {
            rl_tree_lengths->push_back( tree_lengths[i] );
        }

        return new RevVariable( rl_tree_lengths );
    }
    else if ( name == "size" || name == "getNumberSamples" )
    {
        found = true;

        bool post = static_cast<const RlBoolean &>( args[0].getVariable()->getRevObject() ).getValue();

        int n = (int)this->value->size(post);

        return new RevVariable( new Natural( n ) );
    }
    else if ( name == "getBurnin" )
    {
        found = true;

        int n = this->value->getBurnin();

        return new RevVariable( new Natural( n ) );
    }
    else if ( name == "getTree" )
    {
        found = true;

        // get the index which is the only argument for this method
        long i    = static_cast<const Natural &>( args[0].getVariable()->getRevObject() ).getValue() - 1;

        bool post = static_cast<const RlBoolean &>( args[1].getVariable()->getRevObject() ).getValue();
        i += post * this->value->getBurnin();

        const RevBayesCore::Tree &current_tree = this->value->objectAt( i );

        Tree *rl_tree = NULL;
        if ( this->value->isClock() == true )
        {
            rl_tree = new TimeTree( current_tree );
        }
        else
        {
            rl_tree = new BranchLengthTree( current_tree );
        }
        return new RevVariable( rl_tree );
    }
    else if ( name == "getTopologyFrequency" )
    {
        found = true;

        // get the tree which is the only argument for this method
        const RevBayesCore::Tree &current_tree = static_cast<const Tree &>( args[0].getVariable()->getRevObject() ).getValue();
        bool verbose = static_cast<const RlBoolean &>( args[1].getVariable()->getRevObject() ).getValue();
        int f = this->value->getTopologyFrequency( current_tree, verbose );

        return new RevVariable( new Natural( f ) );
    }
    else if ( name == "getUniqueClades" )
    {
        found = true;

        double clade_CI       = static_cast<const Probability &>( args[0].getVariable()->getRevObject() ).getValue();
        bool non_trivial_only  = static_cast<const RlBoolean &>( args[1].getVariable()->getRevObject() ).getValue();
        bool verbose = static_cast<const RlBoolean &>( args[2].getVariable()->getRevObject() ).getValue();

        std::vector<RevBayesCore::Clade> clades = this->value->getUniqueClades(clade_CI, non_trivial_only, verbose);

        ModelVector<Clade> *rl_clades = new ModelVector<Clade>;
        for (size_t i=0; i<clades.size(); ++i)
        {
            rl_clades->push_back( clades[i] );
        }
        return new RevVariable( rl_clades );

    }
    else if ( name == "getUniqueTrees" )
    {
        found = true;

        double tree_CI       = static_cast<const Probability &>( args[0].getVariable()->getRevObject() ).getValue();
        bool verbose = static_cast<const RlBoolean &>( args[1].getVariable()->getRevObject() ).getValue();

        std::vector<RevBayesCore::Tree> trees = this->value->getUniqueTrees(tree_CI, verbose);

        if ( this->value->isClock() == true )
        {
            ModelVector<TimeTree> *rl_trees = new ModelVector<TimeTree>;
            for (size_t i=0; i<trees.size(); ++i)
            {
                rl_trees->push_back( trees[i] );
            }
            return new RevVariable( rl_trees );
        }
        else
        {
            ModelVector<BranchLengthTree> *rl_trees = new ModelVector<BranchLengthTree>;
            for (size_t i=0; i<trees.size(); ++i)
            {
                rl_trees->push_back( trees[i] );
            }
            return new RevVariable( rl_trees );
        }

    }
    else if ( name == "isTreeCovered" )
    {
        found = true;

        // get the tree which is the only argument for this method
        const RevBayesCore::Tree &current_tree = static_cast<const Tree &>( args[0].getVariable()->getRevObject() ).getValue();
        double ci_size = static_cast<const Probability &>( args[1].getVariable()->getRevObject() ).getValue();
        bool verbose = static_cast<const RlBoolean &>( args[2].getVariable()->getRevObject() ).getValue();
        bool cov = this->value->isCoveredInInterval(current_tree, ci_size, verbose);

        return new RevVariable( new RlBoolean( cov ) );
    }
    else if ( name == "getTrees" )
    {
        found = true;

        std::vector<RevBayesCore::Tree> trees = this->value->getValues();
        size_t start_index = this->value->getBurnin();

        if ( this->value->isClock() == true )
        {
            ModelVector<TimeTree> *rl_trees = new ModelVector<TimeTree>;
            for (size_t i=start_index; i<trees.size(); ++i)
            {
                rl_trees->push_back( trees[i] );
            }
            return new RevVariable( rl_trees );
        }
        else
        {
            ModelVector<BranchLengthTree> *rl_trees = new ModelVector<BranchLengthTree>;
            for (size_t i=start_index; i<trees.size(); ++i)
            {
                rl_trees->push_back( trees[i] );
            }
            return new RevVariable( rl_trees );
        }

    }
    else if ( name == "learnTimeCalibratedSBN" )
    {
        found = true;

        const std::string &branch_length_approximation  = static_cast<const RlString &>( args[0].getVariable()->getRevObject() ).getValue();

        // RevBayesCore::SBNParameters sbn = this->value->learnTimeCalibratedSBN(branch_length_approximation);
        RevBayesCore::SBNParameters sbn;
        return new RevVariable( new SBNParameters( sbn ) );

    }
    else if ( name == "learnUnconstrainedSBN" )
    {
        found = true;

        const std::string &branch_length_approximation  = static_cast<const RlString &>( args[0].getVariable()->getRevObject() ).getValue();

        RevBayesCore::SBNParameters sbn = this->value->learnUnconstrainedSBN(branch_length_approximation);
        return new RevVariable( new SBNParameters( sbn ) );

    }

    return RevObject::executeMethod( name, args, found );
}


/** Get Rev type of object */

const std::string& TraceTree::getClassType(void)
{

    static std::string rev_type = "TraceTree";

    return rev_type;
}

/** Get class type spec describing type of object */

const TypeSpec& TraceTree::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( WorkspaceToCoreWrapperObject<RevBayesCore::TraceTree>::getClassTypeSpec() ) );

    return rev_type_spec;
}



/** Return member rules (no members) */

const MemberRules& TraceTree::getParameterRules(void) const
{

    static MemberRules modelMemberRules;
    static bool rules_set = false;

    if ( !rules_set )
    {

        rules_set = true;
    }

    return modelMemberRules;
}


/** Get type spec */

const TypeSpec& TraceTree::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}



void TraceTree::initMethods( void )
{

    ArgumentRules* burninArgRules = new ArgumentRules();
    std::vector<TypeSpec> burninTypes;
    burninTypes.push_back( Probability::getClassTypeSpec() );
    burninTypes.push_back( Integer::getClassTypeSpec() );
    burninArgRules->push_back( new ArgumentRule("burnin",      burninTypes, "The fraction/number of samples to disregard as burnin.", ArgumentRule::BY_VALUE, ArgumentRule::ANY) );
    this->methods.addFunction( new MemberProcedure( "setBurnin", RlUtils::Void, burninArgRules) );

    ArgumentRules* getBurninArgRules = new ArgumentRules();
    this->methods.addFunction( new MemberProcedure( "getBurnin", Natural::getClassTypeSpec(), getBurninArgRules) );

    ArgumentRules* outgroupArgRules = new ArgumentRules();
    outgroupArgRules->push_back( new ArgumentRule("clade", Clade::getClassTypeSpec(), "The (monophyletic) outgroup.", ArgumentRule::BY_VALUE, ArgumentRule::ANY) );
    this->methods.addFunction( new MemberProcedure( "setOutgroup", RlUtils::Void, outgroupArgRules) );

    ArgumentRules* summarizeArgRules = new ArgumentRules();
    summarizeArgRules->push_back( new ArgumentRule("credibleTreeSetSize", Probability::getClassTypeSpec(), "The size of the credible set to print.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Probability(0.95)) );
    summarizeArgRules->push_back( new ArgumentRule("minCladeProbability", Probability::getClassTypeSpec(), "The minimum clade probability used when printing.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Probability(0.05)) );
    summarizeArgRules->push_back( new ArgumentRule("verbose", RlBoolean::getClassTypeSpec(), "Printing verbose output.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true)) );
    this->methods.addFunction( new MemberProcedure( "summarize", RlUtils::Void, summarizeArgRules) );

    ArgumentRules* cladeProbArgRules = new ArgumentRules();
    cladeProbArgRules->push_back( new ArgumentRule("clade", Clade::getClassTypeSpec(), "The (monophyletic) clade.", ArgumentRule::BY_VALUE, ArgumentRule::ANY) );
    cladeProbArgRules->push_back( new ArgumentRule("verbose", RlBoolean::getClassTypeSpec(), "Printing verbose output.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true)) );
    this->methods.addFunction( new MemberProcedure( "cladeProbability", Probability::getClassTypeSpec(), cladeProbArgRules) );

    ArgumentRules* learnTimeCalibratedSBNRules = new ArgumentRules();
    learnTimeCalibratedSBNRules->push_back( new ArgumentRule("nodeAgeApproximationMethod", RlString::getClassTypeSpec(), "gammaRootNodePropBeta", ArgumentRule::BY_VALUE, ArgumentRule::ANY) );
    this->methods.addFunction( new MemberProcedure( "learnTimeCalibratedSBN", SBNParameters::getClassTypeSpec(), learnTimeCalibratedSBNRules) );

    ArgumentRules* learnUnconstrainedSBNRules = new ArgumentRules();
    learnUnconstrainedSBNRules->push_back( new ArgumentRule("branchLengthApproximationMethod", RlString::getClassTypeSpec(), "gammaMOM|lognormalML|lognormalMOM.", ArgumentRule::BY_VALUE, ArgumentRule::ANY) );
    this->methods.addFunction( new MemberProcedure( "learnUnconstrainedSBN", SBNParameters::getClassTypeSpec(), learnUnconstrainedSBNRules) );

    ArgumentRules* getNumberSamplesArgRules = new ArgumentRules();
    getNumberSamplesArgRules->push_back( new ArgumentRule("post", RlBoolean::getClassTypeSpec(), "Get the post-burnin number of samples?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false)) );
    this->methods.addFunction( new MemberProcedure( "getNumberSamples", Natural::getClassTypeSpec(), getNumberSamplesArgRules) );

    ArgumentRules* getSizeArgRules = new ArgumentRules();
    getSizeArgRules->push_back( new ArgumentRule("post", RlBoolean::getClassTypeSpec(), "Get the post-burnin number of samples?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false)) );
    this->methods.addFunction( new MemberProcedure( "size", Natural::getClassTypeSpec(), getSizeArgRules) );

    ArgumentRules* getTreeArgRules = new ArgumentRules();
    getTreeArgRules->push_back( new ArgumentRule("index", Natural::getClassTypeSpec(), "The index of the tree.", ArgumentRule::BY_VALUE, ArgumentRule::ANY) );
    getTreeArgRules->push_back( new ArgumentRule("post", RlBoolean::getClassTypeSpec(), "Use post-burnin indices?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false)) );
    this->methods.addFunction( new MemberProcedure( "getTree", Tree::getClassTypeSpec(), getTreeArgRules) );

    ArgumentRules* getUniqueTreesArgRules = new ArgumentRules();
    getUniqueTreesArgRules->push_back( new ArgumentRule("credibleTreeSetSize", Probability::getClassTypeSpec(), "The size of the credible set.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Probability(0.95)) );
    getUniqueTreesArgRules->push_back( new ArgumentRule("verbose", RlBoolean::getClassTypeSpec(), "Printing verbose output.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true)) );
    this->methods.addFunction( new MemberProcedure( "getUniqueTrees", ModelVector<Tree>::getClassTypeSpec(), getUniqueTreesArgRules) );

    ArgumentRules* get_trees_arg_rules = new ArgumentRules();
    this->methods.addFunction( new MemberProcedure( "getTrees", ModelVector<Tree>::getClassTypeSpec(), get_trees_arg_rules) );

    ArgumentRules* get_unique_clades_arg_rules = new ArgumentRules();
    get_unique_clades_arg_rules->push_back( new ArgumentRule("credibleTreeSetSize", Probability::getClassTypeSpec(), "The size of the credible set.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Probability(0.95)) );
    get_unique_clades_arg_rules->push_back( new ArgumentRule("nonTrivial", RlBoolean::getClassTypeSpec(), "Retrieve only the non-trivial clades.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true)) );
    get_unique_clades_arg_rules->push_back( new ArgumentRule("verbose", RlBoolean::getClassTypeSpec(), "Inlcude only non-trivial clades.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true)) );
    this->methods.addFunction( new MemberProcedure( "getUniqueClades", ModelVector<Clade>::getClassTypeSpec(), get_unique_clades_arg_rules) );

    ArgumentRules* getTopologyFrequencyArgRules = new ArgumentRules();
    getTopologyFrequencyArgRules->push_back( new ArgumentRule("tree", Tree::getClassTypeSpec(), "The tree.", ArgumentRule::BY_VALUE, ArgumentRule::ANY) );
    getTopologyFrequencyArgRules->push_back( new ArgumentRule("verbose", RlBoolean::getClassTypeSpec(), "Printing verbose output.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true)) );
    this->methods.addFunction( new MemberProcedure( "getTopologyFrequency", Natural::getClassTypeSpec(), getTopologyFrequencyArgRules) );


    ArgumentRules* is_covered_arg_rules = new ArgumentRules();
    is_covered_arg_rules->push_back( new ArgumentRule("tree", Tree::getClassTypeSpec(), "The tree.", ArgumentRule::BY_VALUE, ArgumentRule::ANY) );
    is_covered_arg_rules->push_back( new ArgumentRule("ci_size", Probability::getClassTypeSpec(), "The size of the credible interval.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Probability(0.95)) );
    is_covered_arg_rules->push_back( new ArgumentRule("verbose", RlBoolean::getClassTypeSpec(), "Printing verbose output.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true)) );
    this->methods.addFunction( new MemberProcedure( "isTreeCovered", RlBoolean::getClassTypeSpec(), is_covered_arg_rules) );


    ArgumentRules* computePairwiseRFDistanceArgRules = new ArgumentRules();
    computePairwiseRFDistanceArgRules->push_back( new ArgumentRule("credibleTreeSetSize", Probability::getClassTypeSpec(), "The size of the credible set.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Probability(0.95)) );
    computePairwiseRFDistanceArgRules->push_back( new ArgumentRule("verbose", RlBoolean::getClassTypeSpec(), "Printing verbose output.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true)) );
    this->methods.addFunction( new MemberProcedure( "computePairwiseRFDistances", ModelVector<RealPos>::getClassTypeSpec(), computePairwiseRFDistanceArgRules) );

    ArgumentRules* computeTreeLengthsArgRules = new ArgumentRules();
    this->methods.addFunction( new MemberProcedure( "computeTreeLengths", ModelVector<RealPos>::getClassTypeSpec(), computeTreeLengthsArgRules) );

    ArgumentRules* computeEntropyArgRules = new ArgumentRules();
    computeEntropyArgRules->push_back( new ArgumentRule("credibleTreeSetSize", Probability::getClassTypeSpec(), "The size of the credible set.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Probability(0.95)) );
    computeEntropyArgRules->push_back( new ArgumentRule("num_taxa", Natural::getClassTypeSpec(), "The number of taxa in the dataset.", ArgumentRule::BY_VALUE, ArgumentRule::ANY) );
    computeEntropyArgRules->push_back( new ArgumentRule("verbose", RlBoolean::getClassTypeSpec(), "Printing verbose output.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true)) );
    this->methods.addFunction( new MemberProcedure( "computeEntropy", RealPos::getClassTypeSpec(), computeEntropyArgRules) );

}


/** Get type spec */

void TraceTree::printValue(std::ostream &o, bool user) const
{

    o << "TreeTrace (" << getValue().getFileName() << ")";
}


/** Set a member variable */

void TraceTree::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{

    if ( name == "xxx")
    {

    }
    else
    {
        RevObject::setConstParameter(name, var);
    }
}
