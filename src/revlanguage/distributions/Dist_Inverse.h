#ifndef Dist_Inverse_H
#define Dist_Inverse_H

// Need to include headers for all possible values of valType
#include "Integer.h"
#include "Natural.h"
#include "Probability.h"
#include "Real.h"
#include "RealPos.h"
#include "Simplex.h"
#include "AbstractHomologousDiscreteCharacterData.h"
#include "IidDistribution.h"
#include "ModelVector.h"
#include "TypedDistribution.h"
#include "InverseDistribution.h"
#include "TypeSpec.h"

namespace RevLanguage {

    template<typename valType>
    class Dist_Inverse : public TypedDistribution< valType > {
        
    public:
        Dist_Inverse( void ) : TypedDistribution< valType >(), dist( NULL ) {}
        
        virtual ~Dist_Inverse() {}
        
        // Basic utility functions
        // Clone the object
        Dist_Inverse* clone(void) const {
            return new Dist_Inverse(*this);
        }
        
        // Get Rev type
        static const std::string& getClassType(void) {
            static std::string rev_type = "Dist_Inverse";
            return rev_type;
        }
        
        // Get class type spec
        static const TypeSpec& getClassTypeSpec(void) {
            static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution< valType >::getClassTypeSpec() ) );
            return rev_type_spec;
        }
        
        // Get the Rev-name for this distribution
        std::string getDistributionFunctionName(void) const {
            return "inverse";
        }
        
        // Get the type spec of the instance
        const TypeSpec& getTypeSpec(void) const {
            static TypeSpec ts = getClassTypeSpec();
            return ts;
        }
        
        // Get member rules (const)
        const MemberRules& getParameterRules(void) const {
            static MemberRules dist_member_rules;
            static bool rules_set = false;
            
            if ( !rules_set ) {
                dist_member_rules.push_back( new ArgumentRule( "distribution", 
                    TypedDistribution<valType>::getClassTypeSpec(), 
                    "The distribution to invert.", 
                    ArgumentRule::BY_CONSTANT_REFERENCE, 
                    ArgumentRule::ANY ) );
                rules_set = true;
            }
            
            return dist_member_rules;
        }
        
        // Distribution functions you have to override
        RevBayesCore::InverseDistribution<typename valType::valueType>* createDistribution(void) const {
            // get the parameters
            const Distribution& orig_dist = static_cast<const Distribution &>( dist->getRevObject() );
            
            // Cast the single distribution to the specific type expected by InverseDistribution
            RevBayesCore::TypedDistribution<typename valType::valueType>* typedDistPtr = 
                static_cast<RevBayesCore::TypedDistribution<typename valType::valueType>* >( orig_dist.createDistribution() );

            // Create an instance of InverseDistribution using the single distribution
            RevBayesCore::InverseDistribution<typename valType::valueType>* d = 
                new RevBayesCore::InverseDistribution<typename valType::valueType>(*typedDistPtr);

            delete typedDistPtr;
            
            return d;
        }
        
    protected:
        void setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
            if ( name == "distribution" ) {
                dist = var;
            }
            else {
                TypedDistribution< valType >::setConstParameter(name, var);
            }
        }
        
    private:
        RevPtr<const RevVariable> dist;
    };
    
}

#endif // Dist_Inverse_H