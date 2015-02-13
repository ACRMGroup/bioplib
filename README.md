
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

####(2) Install libxml2

If you wish to support PDBML (XML) format files, you need to install
libxml2 from http://http://xmlsoft.org/doanloads.html

This will normally be already installed and available on Linux systems. If not then it is installed on Fedora/CentOS systems using (as root):

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

####(6) For more information, read the file:

        bioplib-X.Y/doc/doxygen/docsrcinput/page_01.dox

or, after doing 'make doxygen', read the formatted version by pointing
a web browser at:

        bioplib-X.Y/doc/html/index.html

