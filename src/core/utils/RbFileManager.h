#ifndef RbFileManager_H
#define RbFileManager_H

#include <fstream>
#include <string>
#include <vector>


namespace RevBayesCore {
    

/**
 * This is the interface for a class that manages files and directories. 
 * It takes advantage of the dirent.h and sys/stat.h libraries. Besides doing
 * some very basic things like opening and closing files for input or output,
 * it also checks for the presence of a directory or file, can recursively
 * list the contents of a directory, and fill in a vector (recursively) with
 * the file names in a directory.
 *
 */
    class RbFileManager {

    public:
                                RbFileManager(void);                                                                                                    //!< Default constructor
                                RbFileManager(const std::string &fn);                                                                                   //!< Constructor with file/directory name
                                RbFileManager(const std::string &pn, const std::string &fn);                                                            //!< Constructor with path name and file/directory name

        static std::string      expandUserDir(std::string path);
        
        void                    closeFile(std::ifstream& strm);                                                                                         //!< Close input file
        void                    closeFile(std::ofstream& strm);                                                                                         //!< Close output file
        void                    createDirectoryForFile(void);                                                                                           //!< Create the directories in which this files exist.
        void                    formatError(std::string& errorStr);                                                                                     //!< Format the error string when (mis)reading files
        const std::string&      getCurrentDirectory(void) const;                                                                                        //!< Returns the default directory for the process
        std::string             getFileExtension(void) const;                                                                                           //!< Returns the name of the file (could be empty)
        const std::string&      getFileName(void) const;                                                                                                //!< Returns the name of the file (could be empty)
        std::string             getFileNameWithoutExtension(void) const;                                                                                //!< Returns the name of the file (could be empty)
        const std::string&      getFilePath(void) const;                                                                                                //!< Returns the name of the path
        const std::string&      getFullFileName(void) const;
        std::string             getFullFilePath(void) const;                                                                                            //!< Returns the name of the path
        std::string             getLastPathComponent(void);
        const std::string&      getNewLine(void) const;                                                                                           //!< Returns the path separator
        const std::string&      getPathSeparator(void) const;                                                                                           //!< Returns the path separator
        std::istream&           safeGetline(std::istream& is, std::string& t);
        std::string             getStringByDeletingLastPathComponent(const std::string& s);
        bool                    isDirectory(void) const;                                                                                                //!< Is this a directory
        bool                    isDirectoryPresent(const std::string &mp) const;                                                                        //!< Checks for presence of a directory
        bool                    isFile(void) const;                                                                                                     //!< Is this a file
        bool                    isFileNamePresent(void) const;                                                                                          //!< Checks whether the file name is present (true) or empty (false)
        bool                    listDirectoryContents(void);                                                                                            //!< Recursively lists the contents of the directory filePath
        bool                    openFile(std::ifstream& strm);                                                                                          //!< Open file for input
        bool                    openFile(std::ofstream& strm);                                                                                          //!< Open file for output
        void                    setFileName(const std::string &s);                                                                                      //!< Setter function for the fileName
        void                    setFilePath(const std::string &s);                                                                                      //!< Setter function for the filePath
        bool                    setStringWithNamesOfFilesInDirectory(std::vector<std::string>& sv, bool recursive=true);                                //!< Recursively fills in a string vector with the contents of the directory filePath
        bool                    setStringWithNamesOfFilesInDirectory(const std::string& dirpath, std::vector<std::string>& sv, bool recursive=true);    //!< Recursively fills in a string vector with the contents of the directory passed in as the argument dirpath
        bool                    testDirectory(void);                                                                                                    //!< Tests whether the directory filePath is present
        bool                    testFile(void);                                                                                                         //!< Tests whether the file fileName is present

    private:

        std::string             getLastPathComponent(const std::string& s);
        bool                    isFilePresent(const std::string &fn) const;                                                                             //!< Checks for the presence of a file
        bool                    isFilePresent(const std::string &mp, const std::string &mf) const;                                                      //!< Checks for the presence of a file
        bool                    listDirectoryContents(const std::string& dirpath);                                                                      //!< Recursively lists the contents of the directory passed in as the argument dirpath
        bool                    makeDirectory(const std::string &dn);
        bool                    parsePathFileNames(const std::string &s);                                                                               //!< Divides a string into the file path and file name
    
        std::string             file_name;                                                                                                               //!< string with file name
        std::string             file_path;                                                                                                               //!< string with file path
        std::string             full_file_name;
        std::string             path_separator;
        std::string             new_line;
    };
    
}



#endif
