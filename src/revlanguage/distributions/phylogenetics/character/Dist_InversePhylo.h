#ifndef Dist_InversePhylo_H
#define Dist_InversePhylo_H

#include "AbstractHomologousDiscreteCharacterData.h"
#include "InversePhyloDistribution.h"
#include "TypedDistribution.h"

namespace RevLanguage {

    class Dist_InversePhylo : public TypedDistribution< AbstractHomologousDiscreteCharacterData > {

    public:
        Dist_InversePhylo( void ) : TypedDistribution< AbstractHomologousDiscreteCharacterData >(), base_distribution( NULL ) {}
        
        virtual ~Dist_InversePhylo() {}

        // Basic utility functions
        // Clone the object
        Dist_InversePhylo* clone(void) const {
            return new Dist_InversePhylo(*this);
        }

        // Get Rev type
        static const std::string& getClassType(void) {
            static std::string rev_type = "Dist_InversePhylo";
            return rev_type;
        }
        
        // Get class type spec
        static const TypeSpec& getClassTypeSpec(void) {
            static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution< AbstractHomologousDiscreteCharacterData >::getClassTypeSpec() ) );
            return rev_type_spec;
        }

        // Get the alternative Rev names (aliases) for the constructor function
        std::vector<std::string> getDistributionFunctionAliases( void ) const {
            // create alternative constructor function names variable that is the same for all instance of this class
            std::vector<std::string> a_names;
            a_names.push_back( "inv" );
            
            return a_names;
        }

        // Get the Rev-name for this distribution
        std::string getDistributionFunctionName(void) const {
            return "inverse";
        }  

        // Get the type spec of the instance
        const TypeSpec& getTypeSpec(void) const {
            static TypeSpec type_spec = getClassTypeSpec();
            return type_spec;
        }  

        // Get member rules (const)
        const MemberRules& getParameterRules(void) const {
            static MemberRules dist_member_rules;
            static bool rules_set = false;
            
            if ( !rules_set ) {
                dist_member_rules.push_back( new ArgumentRule( "distribution", 
                    TypedDistribution<AbstractHomologousDiscreteCharacterData>::getClassTypeSpec(), 
                    "The distribution to invert.", 
                    ArgumentRule::BY_CONSTANT_REFERENCE, 
                    ArgumentRule::ANY ) );
                rules_set = true;
            }
            
            return dist_member_rules;
        }

        // Distribution functions you have to override      
        RevBayesCore::TypedDistribution< RevBayesCore::AbstractHomologousDiscreteCharacterData >* createDistribution(void) const {
            const Distribution& orig_dist = static_cast<const Distribution &>( base_distribution->getRevObject() );
            
            RevBayesCore::TypedDistribution<RevBayesCore::AbstractHomologousDiscreteCharacterData>* typedDistPtr = 
                static_cast<RevBayesCore::TypedDistribution<RevBayesCore::AbstractHomologousDiscreteCharacterData>* >( orig_dist.createDistribution() );

            RevBayesCore::InversePhyloDistribution<RevBayesCore::AbstractHomologousDiscreteCharacterData>* d = new RevBayesCore::InversePhyloDistribution<RevBayesCore::AbstractHomologousDiscreteCharacterData>(*typedDistPtr);

            delete typedDistPtr;
            
            // Here we return an InversePhyloDistribution, but 'disguised' as a TypedDistribution
            // in order to match the base class's function signature.
            return d;
        }

    protected:
        void setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
            if ( name == "distribution" ) {
                base_distribution = var;
            }
            else {
                TypedDistribution< AbstractHomologousDiscreteCharacterData >::setConstParameter(name, var);
            }
        }

    private:
        RevPtr<const RevVariable> base_distribution;
    };

}

#endif