#include "Func_GammaRateModel.h"

#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "RlMixtureModel.h"
#include "RlRateMatrix.h"
#include "RlSimplex.h"
#include "GenericFunction.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "RbException.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "RlBoolean.h"
#include "Integer.h"
#include "Simplex.h"
#include "TypeSpec.h"

#include "DistributionChisq.h"
#include "RbMathFunctions.h"

// wait, so, how do we encode the operations like scale(model), rate(model), etc?
// well, ratemodel should perhaps be its OWN interface!

// how about things like the m3 model?
// also, could we create this using a FUNCTION?

namespace RevBayesCore
{
    struct RateMixtureModel: public MixtureModel
    {
        const RateMatrix* matrix;
        std::vector<double> fraction;
        std::vector<double> rates;

        RateMixtureModel* clone() const {return new RateMixtureModel(*this);}
        void calculateTransitionProbabilities(const Tree& tau,
                                              int node_index,
                                              int m,
                                              TransitionProbabilityMatrix& P) const
        {
            const TopologyNode* node = tau.getNodes()[node_index];

            if ( node->isRoot() == true )
            {
                throw RbException("RateMixtureModel called updateTransitionProbabilities for the root node\n");
            }

            double end_age = node->getAge();

            // if the tree is not a time tree, then the age will be not a number
            if ( RbMath::isFinite(end_age) == false )
            {
                // we assume by default that the end is at time 0
                end_age = 0.0;
            }
            double start_age = end_age + node->getBranchLength();

            matrix->calculateTransitionProbabilities( start_age, end_age,  rates[m], P);
        }

        std::vector<double> getRootFrequencies(int) const
        {
            return matrix->getStationaryFrequencies();
        }

        std::vector<double> componentProbs() const { return fraction;}

        void scale(double factor)
        {
            for(auto& r: rates)
                r *= factor;
        }

        double rate() const
        {
            double r = 0;
            for(int i=0;i<getNumberOfComponents();i++)
                r += fraction[i] * rates[i];
            return r;
        }
        
        RateMixtureModel(const RateMatrix& m, const std::vector<double>& f, const std::vector<double>& r)
            :MixtureModel(f.size(), m.getNumberOfStates()),
             matrix(m.clone()),
             fraction(f),
             rates(r)
        {
            if (f.size() != r.size())
                throw RbException()<<"RateMixtureModel: got "<<f.size()<<" weights but "<<r.size()<<" rates.  These should be equal.";
        }
    };
}

RevBayesCore::MixtureModel* GammaRateModelFunc(const RevBayesCore::RateGenerator& gen, double a, int nCats, RevBayesCore::Boolean median)
{
    using namespace RevBayesCore;

    auto m = dynamic_cast<const RevBayesCore::RateMatrix*>(&gen);
    if (not m)
        throw RbException()<<"GammaRateModel: needs a RateMatrix, not just a RateGenerator!";
    auto& matrix = *m;
    
    if (nCats < 1)
        throw RbException()<<"GammaRateModel: number of components must be at least 1, but got "<<nCats;

    if (a <= 0)
        throw RbException()<<"GammaRateModel: alpha should be positive, but alpha = "<<a;

    double b = a;
    
    double factor = a / b * nCats;

    std::vector<double> gamma_rates(nCats, 0);
    std::vector<double> fraction(nCats, 1.0/nCats);

    if (median) {
        /* the median value for each category is used to represent all of the values
        in that category */
        double interval = 1.0 / (2.0 * nCats);
        for (int i=0; i<nCats; i++) 
            gamma_rates[i] = RbStatistics::ChiSquare::quantile((i * 2.0 + 1.0) * interval, 2.0 * a) / (2.0 * b);
        double t = 0.0;
        for (int i=0; i<nCats; i++) 
            t += gamma_rates[i];
        for (int i=0; i<nCats; i++)     
            gamma_rates[i] *= factor / t;
    }
    else
    {
        /* the mean value for each category is used to represent all of the values
        in that category */
        /* calculate the points in the gamma distribution */
        for (int i=0; i<nCats-1; i++) 
            gamma_rates[i] = RbStatistics::ChiSquare::quantile((i + 1.0) / nCats, 2.0 * a) / (2.0 * b);
        /* calculate the cumulative values */
        for (int i=0; i<nCats-1; i++) 
            gamma_rates[i] = RbMath::incompleteGamma(gamma_rates[i] * b, a + 1.0);
        gamma_rates[nCats-1] = 1.0;
        /* calculate the relative values and rescale */
        for (int i=nCats-1; i>0; i--){
            gamma_rates[i] -= gamma_rates[i-1];
            gamma_rates[i] *= factor;
        }
        gamma_rates[0] *= factor;
    }
    
    return new RevBayesCore::RateMixtureModel(matrix, fraction, gamma_rates);
}

using namespace RevLanguage;


/** default constructor */
Func_GammaRateModel::Func_GammaRateModel( void ) : TypedFunction<MixtureModel>( )
{
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_GammaRateModel* Func_GammaRateModel::clone( void ) const
{
    return new Func_GammaRateModel( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::MixtureModel >* Func_GammaRateModel::createFunction( void ) const
{
    RevBayesCore::TypedDagNode< double >* a = static_cast<const RealPos &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< RevBayesCore::RateGenerator >* matrix = static_cast<const RateGenerator &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< long >* nCats = static_cast<const Integer &>( this->args[2].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< RevBayesCore::Boolean >* median = static_cast<const RlBoolean &>( this->args[3].getVariable()->getRevObject() ).getDagNode();

    return RevBayesCore::generic_function_ptr2< RevBayesCore::MixtureModel >( GammaRateModelFunc, matrix, a, nCats, median);
}


/* Get argument rules */
const ArgumentRules& Func_GammaRateModel::getArgumentRules( void ) const
{
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "submodel"       , RateMatrix::getClassTypeSpec(), "Sub-rate matrix.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "alpha"          , RealPos::getClassTypeSpec(), "The alpha parameter of the gamma distribution.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "n"              , Integer::getClassTypeSpec(), "The number of bins to approximate the gamma distirbution.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Integer(4) ) );
        argumentRules.push_back( new ArgumentRule( "median"         , RlBoolean::getClassTypeSpec(), "Should we use the median (or the mean)", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RlBoolean(true) ) );

        rules_set = true;
    }

    return argumentRules;
}


const std::string& Func_GammaRateModel::getClassType(void)
{

    static std::string rev_type = "Func_GammaRateModel";

    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_GammaRateModel::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_GammaRateModel::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnGamma";

    return f_name;
}


const TypeSpec& Func_GammaRateModel::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}
