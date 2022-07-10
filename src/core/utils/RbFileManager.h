#ifndef RbFileManager_H
#define RbFileManager_H

#include <fstream>
#include <string>
#include <vector>


namespace RevBayesCore {
    
    std::istream&           safeGetline(std::istream& is, std::string& t); //!< Gets one line from a stream
    std::string             expandUserDir(std::string path); //!< Get full path to user directory

    std::string             getPathSeparator(void);

    bool                    isDirectoryPresent(const std::string &mp);  //!< Checks for presence of a directory
    bool                    isFilePresent(const std::string &fn);  //!< Checks for the presence of a file
    bool                    isFilePresent(const std::string &mp, const std::string &mf);  //!< Checks for the presence of a file

    std::string             getLastPathComponent(const std::string& s);  //!< Get last component of given path
    std::string             getStringByDeletingLastPathComponent(const std::string& s);  //!< Get path by removing last component
    bool                    setStringWithNamesOfFilesInDirectory(const std::string& dirpath, std::vector<std::string>& sv, bool recursive=true);  //!< Recursively fills in a string vector with the contents of the directory passed in argument

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


        void                    createDirectoryForFile(void);  //!< Create the directories in the path of full_file_name
        void                    formatError(std::string& errorStr);  //!< Format the error string when (mis)reading files
        std::string             getFileExtension(void) const;  //!< Returns the file extension from file_name
        const std::string&      getFileName(void) const;
        std::string             getFileNameWithoutExtension(void) const;  //!< Returns file_name without extension
        const std::string&      getFilePath(void) const;
        const std::string&      getFullFileName(void) const;
        std::string             getFullFilePath(void) const;  //!< Get absolute file path from file_path
        std::string             getLastPathComponent(void);  //!< Get last component of the full_file_name
        bool                    isDirectory(void) const;  //!< Is full_file_name an existing directory ?
        bool                    isFile(void) const;  //!< Is file_path + file_name an existing file ?
        bool                    isFileNamePresent(void) const;  //!< Checks whether file_name is non-empty
        void                    setFileName(const std::string &s);
        void                    setFilePath(const std::string &s);
        bool                    setStringWithNamesOfFilesInDirectory(std::vector<std::string>& sv, bool recursive=true);  //!< Recursively fills in a string vector with the contents of the directory given by file_path
        bool                    testDirectory(void);  //!< Tests whether the directory given by file_path exists
        bool                    testFile(void);  //!< Tests whether the file given by file_path + file_name exists

    private:

        std::string             file_name; //!< file name
        std::string             file_path; //!< file path
        std::string             full_file_name; //!< full file path (i.e file_path + file_name)
    };
    
}

#endif
