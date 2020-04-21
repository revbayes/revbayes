#ifndef TreeDiscreteCharacterData_h
#define TreeDiscreteCharacterData_h

#include "Cloneable.h"
#include "AbstractHomologousDiscreteCharacterData.h"
#include "Serializable.h"
#include "Tree.h"


namespace RevBayesCore {
    
    /**
     * Class that contains both a tree and discrete character data object.
     *
     * The usage of this class is as the data type for distributions that are 
     * defined jointly on the tree and a character data object.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2016-02-18, version 1.0
     */
    class TreeDiscreteCharacterData : public Tree  {
        
    public:
        TreeDiscreteCharacterData();                                                                                                    //!< Default constructor
//        ~TreeDiscreteCharacterData();                                                                                                    //!< Default constructor
        
        // implemented methods of the Cloneable interface
        TreeDiscreteCharacterData*                              clone(void) const;
        
        // implemented methods of the Serializable interface
        void                                                    initFromFile(const std::string &dir, const std::string &fn);            //!< Read and resurrect this object from a file in its default format.
        void                                                    initFromString(const std::string &s);                                   //!< Serialize (resurrect) the object from a string value
        virtual void                                            writeToFile(const std::string &dir, const std::string &fn) const;

        
        AbstractHomologousDiscreteCharacterData&                getCharacterData(void);
        const AbstractHomologousDiscreteCharacterData&          getCharacterData(void) const;
        std::vector<double>                                     getTimeInStates(void);
        bool                                                    hasCharacterData(void) const;
        void                                                    setCharacterData(AbstractHomologousDiscreteCharacterData *d);
        void                                                    setTimeInStates(std::vector<double> t); 
        void                                                    setTree(const Tree &t);


    protected:
        // Utility functions
        
        // Member variables

        
    private:
        // Utility functions
        
        // Member variables
        AbstractHomologousDiscreteCharacterData                 *character_data;
        std::vector<double>                                     time_in_states;

    };
    
    
}



#endif /* TreeDiscreteCharacterData_h */
