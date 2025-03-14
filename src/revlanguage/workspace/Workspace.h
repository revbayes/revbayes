#ifndef Workspace_H
#define Workspace_H

#include "Environment.h"

#include <map>
#include <ostream>
#include <string>

namespace RevLanguage {

    class AbstractModelObject;
    class AbstractWorkspaceObject;
    class RevObject;
    class Container;
    class Distribution;
    class RandomNumberGenerator;
    class Function;

    typedef std::map<std::string, RevObject*> TypeTable;

    /**
     * @brief Workspace
     *
     * The Workspace class is used for two singleton instances, the global workspace and the user
     * workspace.
     *
     * The global workspace is the base environment for all other environments. It contains all
     * builtin functions and types, as well as any builtin (system) variables. The user workspace
     * is enclosed within the global workspace, and contains all functions, types and global
     * variables defined by the user. Local variables defined by the user are kept in local
     * environments.
     *
     * The workspace has a variable table and a function table, which it inherits from Environment. In addition, it
     * keeps a type table. It provides various types of functionality for
     * storing and retrieving functions, types and member variables and their initializers.
     *
     * The workspace ensures that symbol names are unique for functions and variables by cross-
     * checking the tables. Class names need to be unique but can overlap with function or variable
     * names. Names are not allowed to clash with reserved words like 'for' or 'while'. As far as
     * user-defined functions are concerned, this is guaranteed by the parser.
     *
     * The class declaration contains code to create the global workspace and the user workspace.
     * Copy, assignment and construction of workspaces are prevented by making constructors and
     * assignment operator private.
     *
     * A distribution is a special complex of functions related to a particular probability
     * distribution. When addDistribution is called, the relevant distribution functions are added
     * to the function table.
     *
     */

    class Workspace : public Environment {

    public:
        virtual ~Workspace(void);                                                                                       //!< Destrcutor
        

        // Environment (frame) functions you have to override
        Workspace*                          clone(void) const;                                                          //!< Clone frame
        void                                printValue(std::ostream& o) const;                                          //!< Print table for user

        // Workspace functions
        bool                                addDistribution(Distribution *dist);               //!< Add distribution
        bool                                addType(RevObject *exampleObj);                                             //!< Add type (auto-generated name = rlType)
        bool                                areTypesInitialized(void) const { return typesInitialized; }                //!< Is type table initialized?
        bool                                existsType(const std::string& name) const;                                  //!< Does the type exist in the type table?
        const TypeSpec&                     getClassTypeSpecOfType(const std::string& type) const;                      //!< Get reference to class vector of type
        const TypeTable&                    getTypeTable(void) const;                                                   //!< Get the type table
        void                                initializeGlobalWorkspace(void);                                            //!< Initialize global workspace for types
        RevObject*                          makeNewDefaultObject(const std::string& type) const;                        //!< Make a clone of the template type object
        void                                updateVectorVariables(void);
        
        static std::shared_ptr<Workspace>   globalWorkspacePtr(void) //!< Get global workspace
        {
            static std::shared_ptr<Workspace> globalSpacePtr(new Workspace("GlobalWorkspace"));
            return globalSpacePtr;
        }
        static Workspace&                   globalWorkspace(void) //!< Get global workspace
        {
            return *globalWorkspacePtr();
        }

        static std::shared_ptr<Workspace>   userWorkspacePtr(void)  //!< Get user workspace
        {
            static std::shared_ptr<Workspace> userSpacePtr(new Workspace(globalWorkspacePtr(),"UserWorkspace"));
            return userSpacePtr;
        }

        static Workspace&                   userWorkspace(void) //!< Get user workspace
        {
            return *userWorkspacePtr();
        }


    private:
                                            Workspace(const std::string &n);                                            //!< Workspace without parent
                                            Workspace(const std::shared_ptr<Environment>& parentSpace, const std::string &n);                  //!< Workspace with parent
                                            Workspace(const Workspace& w);                                              //!< Prevent copy

        Workspace&                          operator=(const Workspace& w);                                              //!< Prevent assignment

        void                                initializeBasicGlobalWorkspace(void);                                       //!< Initialize global workspace for basic procedures and IO
        void                                initializeBasicTypeGlobalWorkspace(void);                                   //!< Initialize global workspace for types
        void                                initializeDemographicFunctionGlobalWorkspace(void);                         //!< Initialize global workspace for demographic functions
        void                                initializeDistGlobalWorkspace(void);                                        //!< Initialize global workspace for distributions
        void                                initializeExtraHelp(void);                                                  //!< Initialize the extra help entries
        void                                initializeFuncGlobalWorkspace(void);                                        //!< Initialize global workspace for functions (fnXXX)
        void                                initializeMonitorGlobalWorkspace(void);                                     //!< Initialize global workspace for monitors
        void                                initializeMoveGlobalWorkspace(void);                                        //!< Initialize global workspace for moves
        void                                initializeTypeGlobalWorkspace(void);                                        //!< Initialize global workspace for types
        void                                initializeVectorTypeGlobalWorkspace(void);                                  //!< Initialize global workspace for types

        TypeTable                           typeTable;                                                                  //!< Type table
        bool                                typesInitialized;                                                           //!< Are types initialized?

    };

    
}


#endif

