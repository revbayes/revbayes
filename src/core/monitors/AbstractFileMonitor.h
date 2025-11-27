#ifndef AbstractFileMonitor_H
#define AbstractFileMonitor_H

#include <fstream>
#include <vector>
#include "RbFileManager.h"

#include "Monitor.h"

namespace RevBayesCore {
class DagNode;

    /** @brief Base abstract class for all file monitors
    *
    * File monitors save information to a file about one or several variable DAG node(s).
    */
    class AbstractFileMonitor : public Monitor {
        
    public:
        // Constructors and Destructors
        AbstractFileMonitor(DagNode *n, std::uint64_t g, const path &fname, bool ap=false, bool wv=true);  //!< Constructor with single DAG node
        AbstractFileMonitor(const std::vector<DagNode *> &n, std::uint64_t g, const path &fname, bool ap=false, bool wv=true);  //!< Constructor with vector of DAG node
        AbstractFileMonitor(const AbstractFileMonitor& f);
        
        virtual ~AbstractFileMonitor(void);
        
        // basic methods
        virtual AbstractFileMonitor*        clone(void) const = 0;
        
        // Monitor functions
        void                                addFileExtension(const std::string &s, bool dir);
        virtual void                        monitor(std::uint64_t gen) = 0;
        virtual void                        printHeader(void) = 0;
        
        // FileMonitor functions
        bool                                isFileMonitor( void ) const;
        void                                openStream(bool reopen);
        void                                setAppend(bool tf);   //!< Set if the monitor should append to an existing file
        void                                setPrintVersion(bool tf);  //!< Set flag whether to print the version

        // functions you may want to overwrite
        virtual void                        closeStream(void);
    
    protected:
        std::fstream                        out_stream;  //!< output file stream
        
        // parameters
        path                                filename;  //!< input name of the output file
        path                                working_file_name;  //!< actual output file name, including extension if applicable
        bool                                append;  //!< whether to append to an existing file
        bool                                flatten;  //!< whether vectors should be flattened in the output (i.e each element treated as a separate variable)
        bool                                write_version;  //!< whether to write the version
        
    };
    
}

#endif

