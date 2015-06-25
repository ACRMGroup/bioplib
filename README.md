
BiopLib
=======

                    (c) 1990-2015 SciTech Software
                       Dr. Andrew C.R. Martin,
                  UCL and The University of Reading

                     EMail: andrew@bioinf.org.uk



Bioplib is a library of routines for the manipulation of protein
structure and sequence using the C programming language. In addition,
the term `Bioplib' refers to routines for more general C programming
purposes.

This library is copyright (c) 1990-2015 and was mostly originally
written by Dr. Andrew C.R. Martin while self-employed. Many
enhancements and some additional routines have been written while at
The University of Reading (2000-2003) and at UCL (1993-1999 and 2004
onwards).

BiopLib is licensed under the GPL Version 3. Commercial licences are
also available - see COPYING.DOC.


INSTALLATION INSTRUCTIONS
-------------------------

####(1) Install libxml2

By default, PDBML (XML) format files are supported. If you wish to do
this, you need to install libxml2

If you do not need PDBML (XML) support, then you can skip this step.

This will normally be already installed and available on Linux
systems. If not then it is installed on Fedora/CentOS systems using
(as root):

        yum -y install libxml2 libxml2-devel

or on Debian/Ubuntu systems using:

        sudo apt-get install libxml2 libxml2-dev

On other systems, you will need to install libxml2 manually from 
http://xmlsoft.org/downloads.html


####(2) Unpack the BiopLib distribution file

If you have downloaded a gzipped tar file, do:

        zcat bioplib-X.Y.tar.gz | tar -xvf -
-or-

        gunzip bioplib-X.Y.tar.gz
        tar -xvf bioplib-X.Y.tar
-or- (if you have Gnu tar)

        tar -zxvf bioplib-X.Y.tar.gz

(where X.Y is the major and minor version numbers - e.g. 3.0)

If you have chosen to download a ZIP file, unpack this using

        unzip bioplib-X.Y.zip

This will create a directory called bioplib-X.Y

Enter this directory and then go into the src sub-directory:

        cd bioplib-X.Y/src



####(3) By default, BiopLib will be installed in sub-directories of your home directory

These directories will be created when you install BiopLib if they do
not exist already:

        ~/include
        ~/include/bioplib
        ~/lib
        ~/data

You can also choose to install the files elsewhere, but need to modify
the Makefile as described below.



####(4) Modify the configuration

If you are using the GNU C compiler and wish to install BiopLib in the
default directories and provide PDBML (XML) support, no configuration
changes should be needed and this section can be skipped.

Otherwise, modify the Makefile as required for your system. 

If you have chosen alternative locations for the include, library and data
directories then you will need to change:

- LIBDEST to the directory where you wish to install the static libraries
- INCDEST to the directory where you wish to install the include files
- DATADEST to the directory where you wish to install the data file

Note that the complete path is required, you cannot do ~/lib.

If you wish to use dynamic libraries (see 'Additional installation
options', below), you may also wish to change their location by
changing:

- SHAREDLIBDEST to the directory where you wish to install the shared libraries

If you do **not** require PDBML (XML) support, comment out the relevant 
COPT line from the Makefile.



####(5) Type the commands:

        make 
        make doxygen
        make install
        make installdata




####(6) Set environment variables

If you are using BiopLib routines that access BiopLib data directories, you must set the environment variable DATADIR to point to the directory in which you have installed the BiopLib data files (default $HOME/data)

        sh/bash:
           export DATADIR=$HOME/data

        csh/tcsh:
           setenv DATADIR $HOME/data
(This command should be placed in your .bashrc, .profile, .tcsh or .cshrc file as appropriate for your shell.)

If you are using the BiopLib interactive help support in your
programs, you must set the environment variable HELPDIR to point to
the directory in which you have installed the BiopLib help files
(default $HOME/help)

        sh/bash:
           export HELPDIR=$HOME/data

        csh/tcsh:
           setenv HELPDIR $HOME/data
(This command should be placed in your .bashrc, .profile, .tcsh or .cshrc file as appropriate for your shell.)




####(7) Additional installation options

You can use BiopLib as a set of shared libraries:

        make shared
        make installshared

You can clean up your compilation directory with:

        make clean




####(8) For more information on using BiopLib, read the file:

        bioplib-X.Y/doc/doxygen/docsrcinput/page_01.dox

or, after doing 'make doxygen', read the formatted version by pointing
a web browser at:

        bioplib-X.Y/doc/html/index.html

