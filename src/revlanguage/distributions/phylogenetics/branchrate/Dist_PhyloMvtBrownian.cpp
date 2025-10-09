#include <cstddef>
#include <ostream>
#include <string>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_PhyloMvtBrownian.h"
#include "PhyloMultivariateBrownianProcess.h"
#include "RlMatrixRealSymmetric.h"
#include "RlTimeTree.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlDistribution.h"
#include "TypeSpec.h"

namespace RevBayesCore { class MatrixReal; }
namespace RevBayesCore { class Tree; }
namespace RevBayesCore { template <class valueType> class TypedDagNode; }

using namespace RevLanguage;


Dist_PhyloMvtBrownian* Dist_PhyloMvtBrownian::clone( void ) const
{
    return new Dist_PhyloMvtBrownian(*this);
}


RevBayesCore::PhyloMultivariateBrownianProcess* Dist_PhyloMvtBrownian::createDistribution( void ) const
{
    // get the parameters

    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* tau = static_cast<const TimeTree &>( tree->getRevObject() ).getDagNode();
    
    RevBayesCore::TypedDagNode<RevBayesCore::MatrixReal>* sig  = static_cast<const MatrixRealSymmetric&>( sigma->getRevObject() ).getDagNode();
//    RevBayesCore::TypedDagNode< RbVector<double> >* r  = static_cast<const ModelVector<Real>&>( rootval->getRevObject() ).getDagNode();
    
//    RevBayesCore::PhyloMultivariateBrownianProcess* process    = new RevBayesCore::PhyloMultivariateBrownianProcess( tau, om, r );
    RevBayesCore::PhyloMultivariateBrownianProcess* process    = new RevBayesCore::PhyloMultivariateBrownianProcess( tau, sig );
    
    return process;

}



/* Get Rev type of object */
const std::string& Dist_PhyloMvtBrownian::getClassType(void) {
    
    static std::string rev_type = "Dist_PhyloMvtBrownian";
    
	return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_PhyloMvtBrownian::getClassTypeSpec(void) {
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Distribution::getClassTypeSpec() ) );
    
	return rev_type_spec;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_PhyloMvtBrownian::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "PhyloBrownianMultiVariate";
    
    return d_name;
}


/** Return member rules (no members) */
const MemberRules& Dist_PhyloMvtBrownian::getParameterRules(void) const
{
    
    static MemberRules dist;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        dist.push_back( new ArgumentRule( "tree" , TimeTree::getClassTypeSpec()           , "The tree along which the process evolves.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist.push_back( new ArgumentRule( "sigma", MatrixRealSymmetric::getClassTypeSpec(), "The variance-covariance matrix.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
//        dist.push_back( new ArgumentRule( "rootval", true, Vector<Real>::getClassTypeSpec() ) );
        rules_set = true;
    }
    
    return dist;
}


const TypeSpec& Dist_PhyloMvtBrownian::getTypeSpec( void ) const {
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}




/** Print value for user */

 void Dist_PhyloMvtBrownian::printValue(std::ostream& o) const {
    
    o << " MultivariateBrownian(";
    
    o << "tau=";
    if ( tree != NULL ) {
        o << tree->getName();
    } else {
        o << "?";
    }

     o << ",";
     
     o << "sigma=";
     if ( sigma != NULL ) {
         o << sigma->getName();
     } else {
         o << "?";
     }

     o << ",";
     
     /*
     o << "rootval=";
     if ( rootval != NULL ) {
         o << rootval->getName();
     } else {
         o << "?";
     }
     */
     o << ")";
}



/** Set a member variable */
void Dist_PhyloMvtBrownian::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {

    if ( name == "tree" )
    {
        tree = var;
    }
    else if ( name == "sigma" )
    {
        sigma = var;
    }
    /*
    else if ( name == "rootval" )
    {
        rootval = var;
    }
    */
    else {
        Distribution::setConstParameter(name, var);
    }
}


