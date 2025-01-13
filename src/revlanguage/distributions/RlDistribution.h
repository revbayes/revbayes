/**
 * @file
 * This file contains the declaration of Distribution, which is
 * the interface and abstract base class for RevLanguage distribution objects.
 *
 * @brief Declaration of Distribution
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-08-03 17:51:49 +0200 (Fri, 03 Aug 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @interface Distribution
 * @package distributions
 * @since Version 1.0, 2012-08-06
 *
 * $Id: Function.h 1734 2012-08-03 15:51:49Z hoehna $
 */

#ifndef RlDistribution_H
#define RlDistribution_H

#include "Distribution.h"
#include "RbHelpDistribution.h"
#include "RevObject.h"

namespace RevLanguage {

    class Distribution : public RevObject {
    
    public:
        virtual                                            ~Distribution(void);                                                                 //!< Destructor
    
        // Basic utility functions you have to override
        virtual Distribution*                               clone(void) const = 0;                                                              //!< Clone object
        static const std::string&                           getClassType(void);                                                                 //!< Get Rev type
        static const TypeSpec&                              getClassTypeSpec(void);                                                             //!< Get class type spec
        virtual MethodTable                                 getDistributionMethods( void ) const;                                                                       //!< Get the member methods


        std::vector<std::string>                            getConstructorFunctionAliases(void) const;                                              //!< Get the alternative names used for the constructor function in Rev.
        std::string                                         getConstructorFunctionName(void) const;
        virtual std::vector<std::string>                    getDistributionFunctionAliases(void) const { return std::vector<std::string>(); }   //!< Get the alternative names used for the constructor function in Rev.
        virtual std::string                                 getDistributionFunctionName(void) const = 0;

        
        // Basic utility functions
        std::string                                         getRevDeclaration(void) const;
        virtual const TypeSpec&                             getVariableTypeSpec(void) const = 0;                                                //!< Get the variable type spec of this distribution
        //void                                                printValue(std::ostream& o) const;                                                  //!< Print value for user

        // GUI methods
        std::string                                         getGuiDistributionName(void) { return guiDistributionName; }
        std::string                                         getGuiDistributionToolTip(void) { return guiDistributionToolTip; }
        void                                                setGuiDistributionName(std::string s) { guiDistributionName = s; }
        void                                                setGuiDistributionToolTip(std::string s) { guiDistributionToolTip = s; }

        // Distribution functions you have to override
        virtual RevObject*                                  createRandomVariable(void) const = 0;                                               //!< Create a random variable from this distribution
        virtual RevBayesCore::Distribution*                 createDistribution(void) const = 0;                                                 //!< Create a random variable from this distribution
    
    
    protected:
        Distribution(void);
        
	void                                                markAsTransform();
        // virtual method that may be overwritten
        virtual void                                        addSpecificHelpFields(RevBayesCore::RbHelpEntry* e) const;
        virtual RevBayesCore::RbHelpDistribution*           constructTypeSpecificHelp(void) const;
                
    private:
        
	bool                                                is_transform = false;
	
        std::string                                         guiDistributionName;
        std::string                                         guiDistributionToolTip;

    };
    
}

#endif

