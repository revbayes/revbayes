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
                                                    MISSING_VARIABLE,
                                                    QUIT };         //!< Exception types

                                    RbException(void);                                      //!< Default constructor
                                    RbException(const std::string& msg);                    //!< Default with message 
                                    RbException(ExceptionType type, const std::string& msg="");//!< General constructor

        // Regular functions
        ExceptionType               getExceptionType(void) const;                           //!< Get exception type
        void                        setMessage(const std::string& msg);
        std::string                 getMessage(void) const;
        void                        prepend(const std::string& s);
        void                        print(std::ostream &o) const;                           //!< Print the exception

        template <typename T>
        RbException&                operator<<(const T&);

    private:
	    ExceptionType               exception_type;                                         //!< Exception type
	    std::string                 message;                                                //!< Error message
    
};

template <typename T>
RbException& RbException::operator<<(const T& t) {
  std::ostringstream oss;
  oss<<message<<t;
  message = oss.str();
  return *this;
}


#endif

