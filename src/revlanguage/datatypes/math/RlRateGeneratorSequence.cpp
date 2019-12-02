#include "RlRateGeneratorSequence.h"

#include "ArgumentRule.h"
#include "MemberProcedure.h"
#include "Natural.h"
#include "RlBoolean.h"

using namespace RevLanguage;

RateGeneratorSequence::RateGeneratorSequence(void) : ModelObject<RevBayesCore::RateGeneratorSequence>()
{
    
    // add method for call "x.nChars()" as a function
    ArgumentRules* nCharsArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "nChars", Natural::getClassTypeSpec(), nCharsArgRules) );
    
    // add method for call "x.nStates()" as a function
    ArgumentRules* nStatesArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "nStates", Natural::getClassTypeSpec(), nStatesArgRules) );

}


RateGeneratorSequence::RateGeneratorSequence( RevBayesCore::RateGeneratorSequence *v) : ModelObject<RevBayesCore::RateGeneratorSequence>( v )
{
    
    // add method for call "x.nChars()" as a function
    ArgumentRules* nCharsArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "nChars", Natural::getClassTypeSpec(), nCharsArgRules) );
    
    // add method for call "x.nStates()" as a function
    ArgumentRules* nStatesArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "nStates", Natural::getClassTypeSpec(), nStatesArgRules) );
    
}


RateGeneratorSequence::RateGeneratorSequence( RevBayesCore::TypedDagNode<RevBayesCore::RateGeneratorSequence> *m) : ModelObject<RevBayesCore::RateGeneratorSequence>( m )
{
    
    // add method for call "x.nChars()" as a function
    ArgumentRules* nCharsArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "nChars", Natural::getClassTypeSpec(), nCharsArgRules) );
    
    // add method for call "x.nStates()" as a function
    ArgumentRules* nStatesArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "nStates", Natural::getClassTypeSpec(), nStatesArgRules) );
    
}


RateGeneratorSequence* RateGeneratorSequence::clone() const
{
    return new RateGeneratorSequence( *this );
}


/* Map calls to member methods */
RevPtr<RevVariable> RateGeneratorSequence::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
   
    if (name == "nChars")
    {
        found = true;
        
        int n = (int)this->dag_node->getValue().getNumberOfCharacters();
        return new RevVariable( new Natural(n) );
    }
    else if (name == "nStates")
    {
        found = true;
        
        int n = (int)this->dag_node->getValue().getNumberOfStates();
        return new RevVariable( new Natural(n) );
    }

    return ModelObject<RevBayesCore::RateGeneratorSequence>::executeMethod( name, args, found );
}


/* Get Rev type of object */
const std::string& RateGeneratorSequence::getClassType(void) {
    
    static std::string revType = "RateGeneratorSequence";
    
	return revType;
}

/* Get class type spec describing type of object */
const TypeSpec& RateGeneratorSequence::getClassTypeSpec(void) {
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( RevObject::getClassTypeSpec() ) );
    
	return rev_type_spec;
}


/** Get the type spec of this class. We return a member variable because instances might have different element types. */
const TypeSpec& RateGeneratorSequence::getTypeSpec(void) const {
    
    static TypeSpec type_spec = getClassTypeSpec();
    return type_spec;
}
