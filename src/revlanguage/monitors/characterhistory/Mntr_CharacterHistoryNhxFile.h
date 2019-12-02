/**
 * @file
 * This file contains the declaration of the RevLanguage wrapper of an extended newick file monitor.
 *
 * @brief Declaration of Mntr_CharacterHistoryNhxFile
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date: 2012-08-06 20:14:22 +0200 (Mon, 06 Aug 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2009-11-20, version 1.0
 * @extends RbObject
 *
 * $Id: Real.h 1746 2012-08-06 18:14:22Z hoehna $
 */

#ifndef Mntr_CharacterHistoryNhxFile_H
#define Mntr_CharacterHistoryNhxFile_H

#include "TreeCharacterHistoryNhxMonitor.h"
#include "RlMonitor.h"
#include "TypedDagNode.h"

#include <ostream>
#include <string>

namespace RevLanguage {
    
    class Mntr_CharacterHistoryNhxFile : public Monitor {
        
    public:
        
        Mntr_CharacterHistoryNhxFile(void);                                                                                                 //!< Default constructor
        
        // Basic utility functions
        virtual Mntr_CharacterHistoryNhxFile*       clone(void) const;                                                                      //!< Clone object
        void                                        constructInternalObject(void);                                                          //!< We construct the a new internal monitor.
        static const std::string&                   getClassType(void);                                                                     //!< Get class name
        static const TypeSpec&                      getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                 getMonitorName(void) const;                                                             //!< Get the name used for the constructor function in Rev.
        const MemberRules&                          getParameterRules(void) const;                                                          //!< Get member rules (const)
        virtual const TypeSpec&                     getTypeSpec(void) const;                                                                //!< Get language type of the object
        virtual void                                printValue(std::ostream& o) const;                                                      //!< Print value (for user)
        
    protected:
        
        void                                        setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        std::set<RevPtr<const RevVariable> >        vars;
        RevPtr<const RevVariable>                   filename;
        RevPtr<const RevVariable>                   tree;
        RevPtr<const RevVariable>                   ctmc;
        RevPtr<const RevVariable>                   atlas;
        RevPtr<const RevVariable>                   maxgen;
        RevPtr<const RevVariable>                   samplegen;
        RevPtr<const RevVariable>                   burnin;
        RevPtr<const RevVariable>                   separator;
        RevPtr<const RevVariable>                   prior;
        RevPtr<const RevVariable>                   posterior;
        RevPtr<const RevVariable>                   likelihood;
        RevPtr<const RevVariable>                   type;
        
    };
    
}

#endif
