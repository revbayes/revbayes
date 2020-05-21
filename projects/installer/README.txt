# Installing RevBayes for TensorPhylo using the installer script

The *install.sh* script is designed to automate most of RevBayes installation on Linux, MacOS and Windows (i.e., cygwin).

## How to use the script on Linux or MacOS

### General requirements
To use the script you will need:
* An internet connection (to download dependencies)
* Basic tools that can be installed with your Linux package manager or homebrew for MacOS. That is:
  * *curl*
  * *zip* and *unzip*
  * *make*
  * *cmake*

The script will warn you if one of these tools is unavailable.

### Installing
1. Open a terminal and go in the *projects/installer*. E.g.,
```bash
$ cd projects/installer
```

2. From the *build/installer* folder, execute the following command:
```bash
$ bash install.sh
```

3. The installation will proceed and may take several minutes. *RevBayes* compilation will use several processor and may slightly impact the responsiveness of your system. Please be patient and do not overload your system (e.g., youtube, etc.).

4. Once done, you will have the *rb* executable installed in the current folder. If the installation was unsuccessful, the *log.txt* will contain more information - if you can't make sense of it: submit an issue with the *log.txt* file on the git repository.


## How to use the script on Windows

### Warning
This script will only install the RevBayes executable without the *gtk*-gui.

### General requirements
To use the script you will need:
* An internet connection (to download dependencies)
* A functional installation of [cygwin](https://www.cygwin.com/)

### Installing
1. Clone or copy RevBayes repository somewhere in your cygwin home folder (e.g., */home/username/RevBayes*). Make sure that you pull/download the *dev_mrm* branch.

2. Copy the cygwin *setup-x86_64.exe* executable in the *build/installer/cygwin/* of the RevBayes Project.

3. Using a cygwin terminal, go into *build/installer/cygwin/* and to install all packages required for the installation, run
```bash
$ cd ~/RevBayes
$ cd projects/installer/cygwin
$ bash installPackages.sh
```

4. Move back to *projects/installer*, and to launch the installation of *RevBayes*, run
```bash
$ cd ..
$ bash install.sh
```

5. The installation will proceed and may take a few hours (hang tight). *RevBayes* compilation will use several processor and may slightly impact the responsiveness of your system. Please be patient and do not overload your systeme.e.g, RevBayes runs, youtube, etc.).

6. Once done, you will have the *rb* executable installed in the current folder. If the installation was unsuccessful, the *log.txt* will contain more information - if you can't make sense of it: submit an issue with the *log.txt* file on the git repository.

7. To use *rb.exe* you will have to update Windows environment variables with the dependencies for RevBayes.

### How to edit the environment variables

1. First figure out which library are required by looking at [RevBayes installation tutorial](https://revbayes.github.io/tutorials/tutorial_structure/) or by following the next steps.

2. In your cygwin terminal run, once in the *build/installer/cygwin* folder, run
```bash
$ ldd.exe rb.exe
```

3. Take note of the different library and paths printed out.

4. Open the editor for environment variables:
    *  push the windows key
    *  enter "Edit the system environment..."
    *  click on environment variables
    *  select the *path* variable and edit it

5. Add all the path from step 3. Be careful:
    *  Slash "/" must be replaced by backslash "\\"
    *  */cygdrive/c/System32* should already be on your path - don't change it.
    *  Other paths must be added, e.g.:
          *  Path */usr/x86_64-w64-mingw32/sys-root/mingw/bin/libwinpthread-1.dll*
          *  Must be added as *C:\cygwin64\usr\x86_64-w64-mingw32\sys-root\mingw\bin*
          *  Where *C:\cygwin64* must match the location of your cygwin installation.
