#include "RlContinuousStochasticNode.h"

#include <math.h>
#include <cstddef>
#include <ostream>
#include <string>

#include "RealPos.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MemberProcedure.h"
#include "Real.h"
#include "RevVariable.h"
#include "RlContinuousDistribution.h"
#include "RlDistribution.h"
#include "RlPositiveContinuousDistribution.h"
#include "RlTypedDistribution.h"
#include "RlUtils.h"
#include "StringUtilities.h"
#include "Probability.h" // IWYU pragma: keep

namespace RevBayesCore { class ContinuousDistribution; }

using namespace RevLanguage;

ContinuousStochasticNode::ContinuousStochasticNode( const std::string& n, RevBayesCore::ContinuousDistribution* dist, ContinuousDistribution* rlDist ) :
    RevBayesCore::ContinuousStochasticNode( n, dist ),
    rlDistribution( rlDist )
{
    
    ArgumentRules* clampArgRules = new ArgumentRules();
    clampArgRules->push_back( new ArgumentRule("x", Real::getClassTypeSpec(), "The observed value.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    this->methods.addFunction( new MemberProcedure( "clamp", RlUtils::Void, clampArgRules) );
    
    ArgumentRules* redrawArgRules = new ArgumentRules();
    this->methods.addFunction( new MemberProcedure( "redraw", RlUtils::Void, redrawArgRules) );
    
    ArgumentRules* probArgRules = new ArgumentRules();
    this->methods.addFunction( new MemberProcedure( "probability", RealPos::getClassTypeSpec(), probArgRules) );
    
    ArgumentRules* lnprobArgRules = new ArgumentRules();
    this->methods.addFunction( new MemberProcedure( "lnProbability", Real::getClassTypeSpec(), lnprobArgRules) );
    
    ArgumentRules* setValueArgRules = new ArgumentRules();
    setValueArgRules->push_back( new ArgumentRule("x", Real::getClassTypeSpec(), "The value", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    this->methods.addFunction( new MemberProcedure( "setValue", RlUtils::Void, setValueArgRules) );
    
    ArgumentRules* unclampArgRules = new ArgumentRules();
    this->methods.addFunction( new MemberProcedure( "unclamp", RlUtils::Void, unclampArgRules) );
    
}


ContinuousStochasticNode::ContinuousStochasticNode( const std::string& n, RevBayesCore::ContinuousDistribution* dist, PositiveContinuousDistribution* rlDist ) :
    RevBayesCore::ContinuousStochasticNode( n, dist ),
    rlDistribution( rlDist )
{
    
    ArgumentRules* clampArgRules = new ArgumentRules();
    clampArgRules->push_back( new ArgumentRule("x", Real::getClassTypeSpec(), "The observed value.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    this->methods.addFunction( new MemberProcedure( "clamp", RlUtils::Void, clampArgRules) );
    
    ArgumentRules* redrawArgRules = new ArgumentRules();
    this->methods.addFunction( new MemberProcedure( "redraw", RlUtils::Void, redrawArgRules) );
    
    ArgumentRules* probArgRules = new ArgumentRules();
    this->methods.addFunction( new MemberProcedure( "probability", RealPos::getClassTypeSpec(), probArgRules) );
    
    ArgumentRules* lnprobArgRules = new ArgumentRules();
    this->methods.addFunction( new MemberProcedure( "lnProbability", Real::getClassTypeSpec(), lnprobArgRules) );
    
    ArgumentRules* setValueArgRules = new ArgumentRules();
    setValueArgRules->push_back( new ArgumentRule("x", Real::getClassTypeSpec(), "The value", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    this->methods.addFunction( new MemberProcedure( "setValue", RlUtils::Void, setValueArgRules) );
    
    ArgumentRules* unclampArgRules = new ArgumentRules();
    this->methods.addFunction( new MemberProcedure( "unclamp", RlUtils::Void, unclampArgRules) );
    
}


ContinuousStochasticNode::ContinuousStochasticNode( const std::string& n, RevBayesCore::ContinuousDistribution* dist, TypedDistribution<Probability>* rlDist ) :
    RevBayesCore::ContinuousStochasticNode( n, dist ),
    rlDistribution( rlDist )
{
    
    ArgumentRules* clampArgRules = new ArgumentRules();
    clampArgRules->push_back( new ArgumentRule("x", Real::getClassTypeSpec(), "The observed value.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    this->methods.addFunction( new MemberProcedure( "clamp", RlUtils::Void, clampArgRules) );
    
    ArgumentRules* redrawArgRules = new ArgumentRules();
    this->methods.addFunction( new MemberProcedure( "redraw", RlUtils::Void, redrawArgRules) );
    
    ArgumentRules* probArgRules = new ArgumentRules();
    this->methods.addFunction( new MemberProcedure( "probability", RealPos::getClassTypeSpec(), probArgRules) );
    
    ArgumentRules* lnprobArgRules = new ArgumentRules();
    this->methods.addFunction( new MemberProcedure( "lnProbability", Real::getClassTypeSpec(), lnprobArgRules) );
    
    ArgumentRules* setValueArgRules = new ArgumentRules();
    setValueArgRules->push_back( new ArgumentRule("x", Real::getClassTypeSpec(), "The value", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    this->methods.addFunction( new MemberProcedure( "setValue", RlUtils::Void, setValueArgRules) );
    
    ArgumentRules* unclampArgRules = new ArgumentRules();
    this->methods.addFunction( new MemberProcedure( "unclamp", RlUtils::Void, unclampArgRules) );
    
}


ContinuousStochasticNode::ContinuousStochasticNode( const ContinuousStochasticNode &n ) :
    RevBayesCore::ContinuousStochasticNode( n ),
    rlDistribution( n.rlDistribution->clone() ),
    methods( n.methods )
{
    
}


ContinuousStochasticNode::~ContinuousStochasticNode( void )
{
    
    delete rlDistribution;
    
}


ContinuousStochasticNode& ContinuousStochasticNode::operator=( const ContinuousStochasticNode &n )
{
    
    // check for self-assignment
    if ( this != &n )
    {
        RevBayesCore::ContinuousStochasticNode::operator=(n);
        
        delete rlDistribution;
        
        rlDistribution  = n.rlDistribution->clone();
        methods         = n.methods;
    }
    
    return *this;
}


RevLanguage::ContinuousStochasticNode* RevLanguage::ContinuousStochasticNode::clone( void ) const
{
    
    return new ContinuousStochasticNode( *this );
}


/* Execute calls to member methods */
RevLanguage::RevPtr<RevLanguage::RevVariable> RevLanguage::ContinuousStochasticNode::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
    
    if (name == "clamp")
    {
        
        // we found the corresponding member method
        found = true;
        
        // get the observation
        const double &observation = static_cast<const Real &>( args[0].getVariable()->getRevObject() ).getValue();
        
        // clamp
        this->clamp( new double(observation) );
        
        return NULL;
    }
    else if (name == "lnProbability")
    {
        
        // we found the corresponding member method
        found = true;
        
        return RevPtr<RevVariable>( new RevVariable( new Real( this->getLnProbability() ), "" ) );
    }
    else if (name == "probability")
    {
        
        // we found the corresponding member method
        found = true;
        
        return RevPtr<RevVariable>( new RevVariable( new RealPos( exp( this->getLnProbability() ) ), "" ) );
    }
    else if (name == "redraw")
    {
        
        // we found the corresponding member method
        found = true;
        
        // manually calling redraw allows value to be set
        this->setIgnoreRedraw( false );
        
        // redraw the value
        this->redraw();
        
        return NULL;
    }
    else if (name == "setValue")
    {
        
        // we found the corresponding member method
        found = true;
        
        // get the observation
        const double &observation = static_cast<const Real &>( args[0].getVariable()->getRevObject() ).getValue();
        
        // set value
        this->setValue( new double(observation) );
        
        // mark this node to ignore redraws
        this->setIgnoreRedraw( true );
        
        return NULL;
    }
    else if (name == "unclamp")
    {
        
        // we found the corresponding member method
        found = true;
        
        // Unclamp
        this->unclamp();
        
        return NULL;
    }
    
    found = false;
    
    return NULL;
}


/**
 * Get common member methods.
 */
const RevLanguage::MethodTable& RevLanguage::ContinuousStochasticNode::getMethods( void ) const
{
    
    // return the internal value
    return methods;
}


RevLanguage::Distribution& RevLanguage::ContinuousStochasticNode::getRlDistribution( void )
{
    
    return *rlDistribution;
}


const RevLanguage::Distribution& RevLanguage::ContinuousStochasticNode::getRlDistribution( void ) const
{
    
    return *rlDistribution;
}


/** Print struct for user */
void RevLanguage::ContinuousStochasticNode::printStructureInfo( std::ostream& o, bool verbose ) const
{
    
    o << "_dagType      = Stochastic node (distribution)" << std::endl;
    o << "_distribution = " << rlDistribution->getRevDeclaration() << std::endl;
    o << "_clamped      = " << ( this->clamped ? "TRUE" : "FALSE" ) << std::endl;
    o << "_lnProb       = " << const_cast< ContinuousStochasticNode* >( this )->getLnProbability() << std::endl;
    
    if ( this->touched == true && verbose == true)
    {
        o << "_stored_ln_prob = " << this->stored_ln_prob << std::endl; // const_cast< ContinuousStochasticNode* >( this )->getLnProbability() << std::endl;
    }

    
    o << "_parents      = ";
    this->printParents( o, 16, 70, verbose );
    o << std::endl;
    
    o << "_children     = ";
    this->printChildren( o, 16, 70, verbose );
    o << std::endl;
    
    if ( verbose == true )
    {
        o << "_dagNode      = " << this->name << " <" << this << ">" << std::endl;
        o << "_refCount     = " << this->getReferenceCount() << std::endl;
    }
    
    o << std::endl;
}
