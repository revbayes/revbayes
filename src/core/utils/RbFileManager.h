#ifndef RbFileManager_H
#define RbFileManager_H

#include <fstream>
#include <string>
#include <vector>


namespace RevBayesCore {
    

/** @brief Contains functions for files and directories management
 *
 * This class takes advantage of the dirent.h and sys/stat.h libraries. Besides doing
 * some very basic things like opening and closing files for input or output,
 * it also checks for the presence of a directory or file, can recursively
 * list the contents of a directory, and fill in a vector (recursively) with
 * the file names in a directory.
 *
 */
    class RbFileManager {

    public:
                                RbFileManager(void);
                                RbFileManager(const std::string &fn);  //!< Constructor with full file/directory name
                                RbFileManager(const std::string &pn, const std::string &fn);  //!< Constructor with path name and file/directory name

        static std::string      expandUserDir(std::string path); //!< Get full path to user directory

        void                    closeFile(std::ifstream& strm);  //!< Close input file
        void                    closeFile(std::ofstream& strm);  //!< Close output file
        void                    createDirectoryForFile(void);  //!< Create the directories in the path of full_file_name
        void                    formatError(std::string& errorStr);  //!< Format the error string when (mis)reading files
        std::string             getCurrentDirectory(void) const;  //!< Returns the current directory for the process
        std::string             getFileExtension(void) const;  //!< Returns the file extension from file_name
        const std::string&      getFileName(void) const;
        std::string             getFileNameWithoutExtension(void) const;  //!< Returns file_name without extension
        const std::string&      getFilePath(void) const;
        const std::string&      getFullFileName(void) const;
        std::string             getFullFilePath(void) const;  //!< Get absolute file path from file_path
        std::string             getLastPathComponent(void);  //!< Get last component of the full_file_name
        const std::string&      getNewLine(void) const;
        const std::string&      getPathSeparator(void) const;
        std::istream&           safeGetline(std::istream& is, std::string& t); //!< Gets one line from a stream
        std::string             getStringByDeletingLastPathComponent(const std::string& s);  //!< Get path by removing last component
        bool                    isDirectory(void) const;  //!< Is full_file_name an existing directory ?
        bool                    isDirectoryPresent(const std::string &mp) const;  //!< Checks for presence of a directory
        bool                    isFile(void) const;  //!< Is file_path + file_name an existing file ?
        bool                    isFileNamePresent(void) const;  //!< Checks whether file_name is non-empty
        bool                    listDirectoryContents(void);  //!< Recursively lists the contents of the directory given by file_path
        bool                    openFile(std::ifstream& strm);  //!< Open file for input
        bool                    openFile(std::ofstream& strm);  //!< Open file for output
        void                    setFileName(const std::string &s);
        void                    setFilePath(const std::string &s);
        bool                    setStringWithNamesOfFilesInDirectory(std::vector<std::string>& sv, bool recursive=true);  //!< Recursively fills in a string vector with the contents of the directory given by file_path
        bool                    setStringWithNamesOfFilesInDirectory(const std::string& dirpath, std::vector<std::string>& sv, bool recursive=true);  //!< Recursively fills in a string vector with the contents of the directory passed in argument
        bool                    testDirectory(void);  //!< Tests whether the directory given by file_path exists
        bool                    testFile(void);  //!< Tests whether the file given by file_path + file_name exists

    private:

        std::string             getLastPathComponent(const std::string& s);  //!< Get last component of given path
        bool                    isFilePresent(const std::string &fn) const;  //!< Checks for the presence of a file
        bool                    isFilePresent(const std::string &mp, const std::string &mf) const;  //!< Checks for the presence of a file
        bool                    listDirectoryContents(const std::string& dirpath);  //!< Recursively lists the contents of a directory
        bool                    makeDirectory(const std::string &dn); //!< Creates directory from path
        bool                    parsePathFileNames(const std::string &s);  //!< Divides a string into the file path and file name
    
        std::string             file_name; //!< file name
        std::string             file_path; //!< file path
        std::string             full_file_name; //!< full file path (i.e file_path + file_name)
        std::string             path_separator; //!< path separator (based on OS)
        std::string             new_line; //!< newline character (based on OS)
    };
    
}

#endif
