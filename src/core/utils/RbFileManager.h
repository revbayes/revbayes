#ifndef RbFileManager_H
#define RbFileManager_H

#include <fstream>
#include <string>
#include <vector>

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>

namespace RevBayesCore {

    using namespace boost::filesystem;
    
    std::istream&           safeGetline(std::istream& is, std::string& t); //!< Gets one line from a stream
    path                    expandUserDir(std::string path); //!< Get full path to user directory

    bool                    setStringWithNamesOfFilesInDirectory(const path& dirpath, std::vector<path>& sv, bool recursive=true);  //!< Recursively fills in a string vector with the contents of the directory passed in argument
    path                    appendToStem(const path& p, const std::string& s);
    void                    formatError(const path& p, std::string& errorStr);  //!< Format the error string when (mis)reading files

    void                    createDirectoryForFile(const path& p);

    std::stringstream       readFileAsStringStream(const path& fname);
    std::string             readFileAsString(const path& fname);
}

#endif
