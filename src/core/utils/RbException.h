#ifndef RbException_H
#define RbException_H


#include <iostream>
#include <sstream>

class RbException {

    public:
        // Exception types
        enum                        ExceptionType { DEFAULT,
                                                    BUG,
                                                    MATH_ERROR,
                                                    SKIP_PROPOSAL,
                                                    MISSING_VARIABLE,
                                                    QUIT };         //!< Exception types

                                    RbException(void) = default;                            //!< Default constructor
                                    RbException(const std::string& msg);                    //!< Default with message 
                                    RbException(ExceptionType type, const std::string& msg="");//!< General constructor
                                    RbException(RbException&&) = default;
                                    RbException(const RbException&);

        // Regular functions
        ExceptionType               getExceptionType(void) const;                           //!< Get exception type
        void                        setMessage(const std::string& msg);
        std::string                 getMessage(void) const;
        void                        prepend(const std::string& s);
        void                        print(std::ostream &o) const;                           //!< Print the exception

        template <typename T>
        RbException&                operator<<(const T&);

    private:
	ExceptionType               exception_type = DEFAULT;                               //!< Exception type
	std::ostringstream          message;                                                //!< Error message
    
};

template <typename T>
RbException& RbException::operator<<(const T& t) {
    message<<t;
    return *this;
}


#endif

