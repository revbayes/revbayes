#include <sys/stat.h>
#include <stdio.h>
#include <sys/types.h> // IWYU pragma: keep
#include <iostream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <algorithm>
#include <vector>

#include "RbFileManager.h"
#include "RbSettings.h"
#include "StringUtilities.h"
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>

#ifdef _WIN32
#	include <dirent.h>
#   include <unistd.h>
#   include <windows.h>
#   include <shlwapi.h>
#else
#	include <dirent.h>
#   include <unistd.h>
#endif


using namespace RevBayesCore;
namespace fs = boost::filesystem;

/** Default constructor, creating a file manager object with the file
 path equal to the current (default) directory and an empty file name */
RbFileManager::RbFileManager( void ) :
    file_name( "" ),
    file_path( "" ),
    full_file_name( "" ),
    path_separator( "" ),
    new_line( "" )
{
    
#	ifdef _WIN32
    path_separator = "\\";
    new_line = "\r\n";
#   else
    path_separator = "/";
    new_line = "\n";
#   endif
    
    file_path = RbSettings::userSettings().getWorkingDirectory();
    full_file_name = file_path;
    if ( full_file_name != "")
    {
        full_file_name += path_separator;
    }
    
    full_file_name += file_name;
    
    
#    ifdef _WIN32
    StringUtilities::replaceSubstring(full_file_name,"/","\\");
#   endif
    
}


/** Constructor with full file/directory name */
RbFileManager::RbFileManager(const std::string &fn) :
    file_name( fn ),
    file_path( "" ),
    full_file_name( "" ),
    path_separator( "" )
{
    
#	ifdef _WIN32
    path_separator = "\\";
#	else
    path_separator = "/";
#   endif
    
    parsePathFileNames( fn );
    
#   ifdef _WIN32
    StringUtilities::replaceSubstring(file_path,"/","\\");
#   endif
    
    full_file_name = file_path;
    if ( full_file_name != "")
    {
        full_file_name += path_separator;
    }
    
    full_file_name += file_name;
    
#   ifdef _WIN32
    StringUtilities::replaceSubstring(full_file_name,"/","\\");
#   endif
    
}


/** Constructor with path name and file/directory name */
RbFileManager::RbFileManager(const std::string &pn, const std::string &fn) :
    file_name( fn ),
    file_path( pn ),
    full_file_name( "" ),
    path_separator( "" )
{
    
#	ifdef _WIN32
    path_separator = "\\";
    new_line = "\r\n";
#	else
    path_separator = "/";
    new_line = "\n";
#   endif
    
    // set the path and file for the string
    std::string tmp = pn + path_separator + fn;
    parsePathFileNames( tmp );

#   ifdef _WIN32
    StringUtilities::replaceSubstring(file_path,"/","\\");
#   endif
    
    full_file_name = file_path;
    if ( full_file_name != "")
    {
        full_file_name += path_separator;
    }
    
    full_file_name += file_name;
    
#   ifdef _WIN32
    StringUtilities::replaceSubstring(full_file_name,"/","\\");
#   endif
    
}


/** Closes input file
 * @param strm stream to close
 * */
void RbFileManager::closeFile(std::ifstream& strm)
{   
    strm.close();
}


/** Closes output file
 * @param strm stream to close
 * */
void RbFileManager::closeFile(std::ofstream& strm)
{    
    strm.close();
}


/**
 * Recursively create the directories which are present in the full file name.
 */
void RbFileManager::createDirectoryForFile( void )
{
    
    std::string dir_path = getStringByDeletingLastPathComponent( full_file_name );
    
    std::vector<std::string> pathComponents;
    StringUtilities::stringSplit(dir_path, path_separator, pathComponents);
    
    std::string directoryName = "";
    for ( std::vector<std::string>::const_iterator it=pathComponents.begin(); it != pathComponents.end(); ++it)
    {
        directoryName += *it;
        
        if ( isDirectoryPresent( directoryName ) == false )
        {
            bool success = makeDirectory( directoryName );
            if ( success == false )
            {
                std::cerr << "Failed to build directory with name \"" << directoryName << "\"." << std::endl;
            }
        }
        
        directoryName += path_separator;
    }
    
}


/**
 * Portable code to get a full path by expanding the user home directory.
 *
 * @param path relative path
 * @return full path including user home directory
 */
std::string RbFileManager::expandUserDir(std::string path)
{
    if ( !path.empty() && path[0] == '~')
    {
        char const* home = getenv("HOME");
        
        if (home or ((home = getenv("USERPROFILE"))))
        {
            path.replace(0, 1, home);
        }
    }
    else if ( path.empty() == false )
    {
        char const *hdrive = getenv("HOMEDRIVE"), *hpath = getenv("HOMEPATH");
        if ( hdrive != NULL )
        {
# ifdef _WIN32
            path = std::string(hdrive) + hpath + "\\" + path;
# else
            path.replace(0, 1, std::string(hdrive) + hpath);
# endif
        }
        
    }
    
    return path;
}


/** Format an error exception string for problems specifying the file/path name
 * @param[out] errorStr string to store the formatted error
*/
void RbFileManager::formatError(std::string& errorStr)
{
    
    bool file_nameProvided    = isFileNamePresent();
    bool isfile_nameGood      = testFile();
    bool isDirectoryNameGood = testDirectory();
    
    if ( file_nameProvided == false && isDirectoryNameGood == false )
    {
        errorStr += "Could not read contents of directory \"" + getFilePath() + "\" because the directory does not exist";
    }
    else if (file_nameProvided == true && (isfile_nameGood == false || isDirectoryNameGood == false))
    {
        errorStr += "Could not read file named \"" + getFileName() + "\" in directory named \"" + getFilePath() + "\" ";
        if (isfile_nameGood == false && isDirectoryNameGood == true)
        {
            errorStr += "because the file does not exist";
        }
        else if (isfile_nameGood == true && isDirectoryNameGood == false)
        {
            errorStr += "because the directory does not exist";
        }
        else
        {
            errorStr += "because neither the directory nor the file exist";
        }
    }
    
}


std::string RbFileManager::getCurrentDirectory( void ) const
{
    return fs::current_path().make_preferred().string();
}


std::string RbFileManager::getFileExtension( void ) const
{
    std::vector<std::string> tokens;
    StringUtilities::stringSplit(file_name,".",tokens);
    return tokens[tokens.size()-1];
}


const std::string& RbFileManager::getFileName( void ) const
{
    return file_name;
}


std::string RbFileManager::getFileNameWithoutExtension( void ) const
{
    std::vector<std::string> tokens;
    StringUtilities::stringSplit(file_name,".",tokens);
    std::string name = "";
    
    for (size_t i = 0; i < tokens.size()-1; ++i)
    {
        if(i > 0) name += ".";
        name += tokens[i];
    }
    
    return name;
}


const std::string& RbFileManager::getFilePath( void ) const
{
    return file_path;
}


const std::string& RbFileManager::getFullFileName( void ) const
{
    return full_file_name;
}


/** Get absolute file path from file_path
 * @return absolute path
 */
std::string RbFileManager::getFullFilePath( void ) const
{
    
    std::string fullfile_path = file_path;
    
    // check if file_path is relative or absolute
    // add current working path only if relative
#   ifdef RB_XCODE
    if ( file_path.size() > 0 && path_separator[0] != file_path[0] )
#    else
    boost::filesystem::path tmp_file = boost::filesystem::path(file_path);
    if ( tmp_file.is_absolute() == false )
#    endif
    {        
        fullfile_path = RbSettings::userSettings().getWorkingDirectory() + path_separator + file_path;        
    }
    
    return fullfile_path;
    
}


/** Get the last path component of full_file_name
 * @note any trailing path separator is removed, so x/y/z/ will return z
 * @return last path component
 */
std::string RbFileManager::getLastPathComponent( void )
{
    
    std::string tmp = full_file_name;
    if ( tmp[tmp.size()-1] == path_separator[0] )
    {
        tmp = tmp.substr(0,tmp.size()-1);
    }
    return getLastPathComponent( tmp );
}


/** Get the last path component of a path
 * @note any trailing path separator is NOT removed, so x/y/z/ will return an empty string
 * @param s input path
 * @return last path component
 */
std::string RbFileManager::getLastPathComponent(const std::string& s)
{
    
    std::string tempS = s;
    size_t location = tempS.find_last_of( path_separator );
    if ( location == std::string::npos )
    {
        // There is no path in this string. We
        // must have only the file name.
        return tempS;
    }
    else if ( location == tempS.length() - 1 )
    {
        // It looks like the last character is "/", which
        // means that no file name has been provided.
        return "";
    }
    else
    {
        /* We can divide the path into the path and the file. */
        std::string lpc = tempS.substr( location+1, tempS.length()-location-1 );
        return lpc;
    }
    
    return "";
}


const std::string& RbFileManager::getNewLine(void) const
{
    return new_line;
}


const std::string& RbFileManager::getPathSeparator( void ) const
{
    return path_separator;
}


/** Removes the last path component from a path
 * @note any trailing path separator is NOT removed, so x/y/z/ will return x/y/z
 * @return string without the last path component
 */
std::string RbFileManager::getStringByDeletingLastPathComponent(const std::string& s)
{
    
    std::string tempS = s;
    size_t location = tempS.find_last_of( path_separator );
    
    if ( location == std::string::npos )
    {
        /* There is no path in this string. We
         must have only the file name. */
        return "";
    }
    else if ( location == tempS.length() - 1 )
    {
        /* It looks like the last character is "/", which
         means that no file name has been provided. */
        return tempS;
    }
    else
    {
        /* We can divide the path into the path and the file. */
        tempS.erase( location );
        return tempS;
    }
    
    return "";
}


/** Checks whether full_file_name is a path to an existing directory */
bool RbFileManager::isDirectory( void ) const
{
    
    bool tf = isDirectoryPresent(full_file_name);
    return tf;
}


/** Tests whether a directory is present (and is a directory)
 * @param mp path to check
 * @return result of the test
 */
bool RbFileManager::isDirectoryPresent(const std::string &mp) const
{
    
    if ( mp == "" )
    {
        return true;
    }
    
    struct stat info;
    
    if ( stat( mp.c_str(), &info ) != 0)
    {
        return false;
    }
    else if (info.st_mode & S_IFDIR)
    {
        return true;
    }
    else
    {
        return false;
    }

    // Sebastian (20150416): This was temporary code to solve the problem that sometimes this function errornously did see a directory.
    // I keep it here for a bit as a reference if the problem re-occurs.
    
    //    DIR* d = opendir( mp.c_str() );
    //
    //	if ( d == NULL )
    //    {
    //		return false;
    //    }
    //
    //    struct stat fInfo;
    //    if ( !stat(mp.c_str(), &fInfo) )
    //    {
    //        if ( S_ISDIR(fInfo.st_mode) )
    //        {
    //            closedir(d);
    //
    //            return true;
    //        }
    //        else
    //        {
    //            closedir(d);
    //
    //            return false;
    //        }
    //
    //    }
    //
    //    closedir(d);
    //
    //    return false;
}


/** Checks whether the path given by file_path + file_name is a path to an existing file */
bool RbFileManager::isFile( void ) const
{   
    return isFilePresent(file_path, file_name);
}


/** Checks whether the file name is non-empty */
bool RbFileManager::isFileNamePresent(void) const
{
    
    if ( file_name == "" )
    {
        return false;
    }
    
    return true;
}


/** Checks whether a file passed in as its path and file name components is present (and a file)
 * @param mp file path - if empty, set to the working directory
 * @param mf file name
 * @return whether the file exists
*/
bool RbFileManager::isFilePresent(const std::string &mp, const std::string &mf) const
{ 
    
    std::string f = mp;
    if ( mp == "" )
    {
        f = RbSettings::userSettings().getWorkingDirectory();
    }
    
    f += path_separator + mf;
    
    return isFilePresent(f);
    
}

/** Checks whether a file passed in as its full path is present (and is a file)
 * @param fn full file path
 * @return whether the file exists
*/
bool RbFileManager::isFilePresent(const std::string &fn) const
{
    
    struct stat fInfo;
    if ( !stat(fn.c_str(), &fInfo) )
    {
        if ( S_ISDIR(fInfo.st_mode) )
        {
            return false;
        }
        else
        {
            return true;
        }
        
    }
    
    return false;
    
}


/** Recursively lists the contents of the directory given by the member variable file_path
 * @return true
 */
bool RbFileManager::listDirectoryContents(void)
{    
    return listDirectoryContents(file_path);
}


/** Recursively lists the contents of a directory
 * @param dirpath path to the directory
 * @return true
 */
bool RbFileManager::listDirectoryContents(const std::string& dirpath)
{
    
    DIR* dir = opendir( dirpath.c_str() );
    if (dir)
    {
        struct dirent* entry;
        while ((entry = readdir( dir )))
        {
            struct stat entryinfo;
            std::string entryname = entry->d_name;
            std::string entrypath = dirpath + path_separator + entryname;
            if (!stat( entrypath.c_str(), &entryinfo ))
            {
                if (S_ISDIR( entryinfo.st_mode ))
                {
                    if (entryname == "..")
                    {}
                    else if (entryname == "." )
                    {}
                    else
                    {
                        listDirectoryContents( entrypath );
                    }
                }
                else
                {}
            }
        }
        closedir( dir );
    }
    
    return true;
}


/**
 * Make a directory.
 * @param dn path to the new directory
 * @return true on Windows, whether the operation was successful otherwise
 */
bool RbFileManager::makeDirectory(const std::string &dn)
{
    
#	ifdef _WIN32
    
    CreateDirectory(dn.c_str(), NULL);
    
    return true;
    
#   else
    
    std::string cmd = "mkdir -p \"" + dn + "\"";
    
    return ( system( cmd.c_str() ) == 0 );
    
#   endif
}


/** Opens a file for input
 * @param strm stream to associate with the file
 * @return whether the operation was successful
*/
bool RbFileManager::openFile(std::ifstream& strm)
{
    
    // concatenate path and file name
    std::string file_pathName = file_path + path_separator + file_name;
        
    // here we assume that the presence of the path/file has
    // been checked elsewhere
    strm.open( file_pathName.c_str(), std::ios::in );
    if ( !strm )
    {
        return false;
    }
    
    return true;
}


/** Opens a file for output
 * @param strm stream to associate with the file
 * @return whether the operation was successful
*/
bool RbFileManager::openFile(std::ofstream& strm)
{
    
    // concatenate path and file name
    std::string file_pathName = file_path + path_separator + file_name;
    
    // here we assume that the presence of the path/file has
    // been checked elsewhere
    strm.open( file_pathName.c_str(), std::ios::out );
    if ( !strm )
    {
        return false;
    }
    
    return true;
}


/** Divides a string into the path and file name components
 * and sets member variables to those values
 *
 * @param input_string string to parse
 * @return whether a file name was found
*/
bool RbFileManager::parsePathFileNames(const std::string &input_string)
{
    std::string name = input_string;
    
#	ifdef _WIN32
    StringUtilities::replaceSubstring(name,"/","\\");
#   endif
    
#    ifdef RB_XCODE
    //    std::filesystem::path path(winPathString); // Construct the path from a string.
    if ( name.size() > 0 && name[0] != '/' )
    {
        name = getCurrentDirectory() + path_separator + input_string;
    }
#    else
    boost::filesystem::path tmp_file = boost::filesystem::path(name);
    //    std::filesystem::path path(winPathString); // Construct the path from a string.
    if ( tmp_file.is_absolute() == false )
    {
        name = getCurrentDirectory() + path_separator + input_string;
    }
#    endif

    
    // check if the path is a good one
    bool isDPresent = isDirectoryPresent(name);
    bool isFPresent = isFilePresent(name);
    
    if ( name.length() > 0 && isDPresent == true && isFPresent == false)
    {
        file_name = "";
        size_t location = name.find_last_of( path_separator );
        if ( location == name.length() - 1 )
        {
            name.erase( location );
        }
        file_path = name;
        return true;
    }
    
    // the string that is supposed to hold the path/file information is empty.
    if ( name.length() == 0 )
    {
        file_path = getCurrentDirectory();
        return false;
    }
    
    // Find the location of the last "/".
    //This is where we will divide the path/file string into two.
    size_t location = StringUtilities::findLastOf( name, path_separator[0] );
    
    if ( location == std::string::npos )
    {
        /* There is no path in this string. We
         must have only the file name, and the
         file should be in our current directory. */
        file_name = name;
        file_path = getCurrentDirectory();
    }
    else if ( location == name.length() - 1 )
    {
        // It looks like the last character is "/", which
        // means that no file name has been provided. However,
        // it also means that the directory that has been provided
        // is not valid, otherwise it would have tested as
        // being present (above).
        file_name = "";
        file_path = getCurrentDirectory();
        return false;
    }
    else
    {
        // We can divide the path into the path and the file.
        file_name = name.substr( location+1, name.length()-location-1 );
        name.erase( location );
        file_path = name;
        
    }
    
    return true;
}


/** Get line while safely handling cross-platform line endings.
 *  Modified from: https://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf
 *
 * @param[in] is stream from which to read
 * @param[out] t string to store the line read
 *
 * @return the input stream
 * */
std::istream& RbFileManager::safeGetline(std::istream& is, std::string& t)
{
    t.clear();
    
    // The characters in the stream are read one-by-one using a std::streambuf.
    // That is faster than reading them one-by-one using the std::istream.
    // Code that uses streambuf this way must be guarded by a sentry object.
    // The sentry object performs various tasks,
    // such as thread synchronization and updating the stream state.
    
    std::istream::sentry se(is, true);
    std::streambuf* sb = is.rdbuf();
    
    for(;;) {
        int c = sb->sbumpc();
        switch (c) {
            case '\n':
                return is;
            case '\r':
                if(sb->sgetc() == '\n')
                    sb->sbumpc();
                return is;
            case EOF:
                // Also handle the case when the last line has no line ending
                if(t.empty())
                    is.setstate(std::ios::eofbit);
                return is;
            default:
                t += (char)c;
        }
    }
}


void RbFileManager::setFileName(std::string const &s)
{
    file_name = s;
}


void RbFileManager::setFilePath(std::string const &s)
{
    file_path = s;
#	ifdef _WIN32
    StringUtilities::replaceSubstring(file_path,"/","\\");
#   endif
    
}

/** Fills in a vector with the names of the files in the directory file_path
 *
 * @param[out] sv vector to store the names
 * @param[in] recursive whether to list recursively, default true
 *
 * @return true
*/
bool RbFileManager::setStringWithNamesOfFilesInDirectory(std::vector<std::string>& sv, bool recursive)
{
    
    return setStringWithNamesOfFilesInDirectory(file_path, sv, recursive);
}


/** Fills in a vector with the names of the files in a directory
 *
 * @param[in] dirpath path to the directory to be listed
 * @param[out] sv vector to store the names
 * @param[in] recursive whether to list recursively, default true
 *
 * @return true
*/
bool RbFileManager::setStringWithNamesOfFilesInDirectory(const std::string& dirpath, std::vector<std::string>& sv, bool recursive)
{
        
    DIR* dir = opendir( dirpath.c_str() );
    if (dir)
    {
        struct dirent* entry;
        while ( (entry = readdir( dir )) )
        {
            struct stat entryinfo;
            std::string entryname = entry->d_name;
            std::string entrypath = dirpath + path_separator + entryname;
            
            bool skip_me = false;

#ifdef _WIN32
            if (stat( entrypath.c_str(), &entryinfo )) {
              // if this returned a non-zero value, something is wrong
              skip_me = true;
            }
#else
            if (!lstat( entrypath.c_str(), &entryinfo ))
            {
                
                // avoid recursing into symlinks that point to a directory above us
                if ( S_ISLNK( entryinfo.st_mode ) ) {
                    char *linkname = (char*) malloc(entryinfo.st_size + 1);
                    ssize_t r = readlink(entrypath.c_str(), linkname, entryinfo.st_size + 1);
                    if (r < 0 || r > entryinfo.st_size) {
                        // error
                        skip_me = true;
                    } else {
                        linkname[entryinfo.st_size] = '\0';
                        if (strspn(linkname, "./") == entryinfo.st_size) {
                            // this symlink consists entirely of dots and dashes and is likely a reference to a directory above us
                            skip_me = true;
                        } else {
                            // replace entryinfo with info from stat
                            if ( stat( entrypath.c_str(), &entryinfo ) ) {
                                // if this returned a non-zero value, something is wrong
                                skip_me = true;
                            }
                        }
                    }
                    free(linkname);
                }

            } else {
              // if this returned a non-zero value, something is wrong
              skip_me = true;
            }
#endif

            if (!skip_me) {

                if (entryname == "..")
                {
                    ;
                }
                else if (entryname == "." || entryname[0] == '.')
                {
                    ;
                }
                else if ( recursive == true && S_ISDIR( entryinfo.st_mode ) )
                {
                    setStringWithNamesOfFilesInDirectory( entrypath, sv );
                }
                else
                {
                    sv.push_back( entrypath );
                }
            }
            
        }
        
        closedir( dir );
    }
    
    // make sure that the file names are sorted
    std::sort(sv.begin(), sv.end());
    
    return true;
}


/** Tests whether the directory specified in the object exists
 * @return true if file_path exists
 */
bool RbFileManager::testDirectory(void)
{   
    return isDirectoryPresent(file_path);
}


/** Tests whether the file specified in the object exists
 * @return true if file_path + file_name exists
 */
bool RbFileManager::testFile(void)
{   
    return isFilePresent(file_path, file_name);
}
