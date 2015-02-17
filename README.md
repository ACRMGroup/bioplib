
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

####(1) Create sub-directories of your home directory called include and
lib:

        mkdir ~/include
        mkdir ~/include/bioplib
        mkdir ~/lib
        mkdir ~/data

####(2) Install libxml2

By default, PDBML (XML) format files are supported. If you wish to do
this, you need to install libxml2 from
http://http://xmlsoft.org/downloads.html

If you do not need PDBML (XML) support, then you can skip this step.

This will normally be already installed and available on Linux
systems. If not then it is installed on Fedora/CentOS systems using
(as root):

        yum -y install libxml2 libxml2-devel

of on Debian/Ubuntu systems using:

        sudo apt-get install libsml2 libxml2-dev

####(3) Unpack the tar file:

        zcat bioplib-X.Y.tar.gz | tar -xvf -
-or-

        gunzip bioplib-X.Y.tar.gz
        tar -xvf bioplib-X.Y.tar
-or- (if you have Gnu tar)

        tar -zxvf bioplib-X.Y.tar.gz

(where X.Y is the major and minor version numbers - e.g. 3.0)

You may also have chosen to download a ZIP file, in which case this is unpacked using

        unzip bioplib-X.Y.zip


####(4) This will create a directory called bioplib-X.Y.  

Enter this directory and then go into the src sub-directory:

        cd bioplib-X.Y/src

Modify the Makefile as required for your system. If you are using the
GNU C compiler and have followed the directions above, no changes
should be needed. 

If you have chosen alternative locations for the include and library
directories then you will need to change LIBDEST to be the library
directory you have created and INCDEST to the include directory you
have created (N.B. the complete path is required, you can't do ~/lib).

If you do not require PDBML (XML) support, comment out the relevant 
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
(this command should be placed in your .bashrc, .profile, .tcsh or .cshrc file as appropriate for your shell).

If you are using the BiopLib interactive help support in your
programs, you must set the environment variable HELPDIR to point to
the directory in which you have installed the BiopLib help files
(default $HOME/help)

        sh/bash:
           export HELPDIR=$HOME/data

        csh/tcsh:
           setenv HELPDIR $HOME/data
(this command should be placed in your .bashrc, .profile, .tcsh or .cshrc file as appropriate for your shell).




####(7) Additional installation options

You can modify the installation directories as desired by editing the first few lines of the Makefile.

You can also use BiopLib as a set of shared libraries:

        make shared
        make installshared

You can clean up your compilation directory with:

        make clean

####(8) For more information on using BiopLib, read the file:

        bioplib-X.Y/doc/doxygen/docsrcinput/page_01.dox

or, after doing 'make doxygen', read the formatted version by pointing
a web browser at:

        bioplib-X.Y/doc/html/index.html

