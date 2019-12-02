/**
 * @file
 * This file contains the declaration of UserInterface, which is
 * the base class for different user interfaces.
 *
 * @brief Declaration of UserInterface
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @extends Frame
 * @package parser
 * @version 1.0
 * @since version 1.0 2009-09-02
 *
 * $Id$
 */

#ifndef UserInterface_H
#define UserInterface_H

#include "RlUserInterfaceOutputStream.h"

#include <iostream>
#include <sstream>

namespace RevLanguage {

#define RBOUT(m) RevLanguage::UserInterface::userInterface().output((m))

class UserInterface {

    public:
        bool                        ask(std::string msg);                                       //!< Ask user a question
        bool                        initialize(void) { return true; }                           //!< Initialize interface
        void                        output(std::string msg);                                    //!< Display message from string
        void                        output(std::string msg, const bool hasPadding);             //!< Display message from string with control of padding
        void                        output(std::ostringstream msg);                             //!< Display message from stringstream
        void                        setOutputStream(UserInterfaceOutputStream *o);              //!< Set the output stream for the current interface
        static UserInterface&       userInterface(void)                                         //!< Get the user interface
		                               {
		                               static UserInterface theInterface = UserInterface();
		                               return theInterface;
		                               }

    protected:
                                    UserInterface(void);                                        //!< Prevent construction
                                    UserInterface(const UserInterface& x);                      //!< Prevent copy construction
        virtual                    ~UserInterface(void) {}                                      //!< Destructor
        UserInterface&              operator=(const UserInterface& w) { return (*this); }       //!< Prevent assignment
    
        int                         process_id;
        UserInterfaceOutputStream*  output_stream;
};
    
}

#endif

