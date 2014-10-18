
BiopLib
=======

                    (c) 1990-2014 SciTech Software
                       Dr. Andrew C.R. Martin,
                  UCL and The University of Reading

                     EMail: andrew@bioinf.org.uk



Bioplib is a library of routines for the manipulation of protein
structure and sequence using the C programming language. In addition,
the term `Bioplib' refers to routines for more general C programming
purposes.

This library is copyright (c) 1990-2014 and was mostly originally
written by Dr. Andrew C.R. Martin while self-employed. Many
enhancements and some additional routines have been written while at
The University of Reading (2000-2003) and at UCL (1993-1999 and 2004
onwards).

BiopLib is licensed under the GPL Version 3. Commercial licences are
also available - see COPYING.DOC.


INSTALLATION INSTRUCTIONS
-------------------------

####(1) Create sub-directories of your main directory called include and
lib:

        mkdir ~/include
        mkdir ~/include/bioplib
        mkdir ~/lib

####(2) Unpack the tar file:

        zcat bioplib-X.Y.tar.gz | tar -xvf -
-or-
        gunzip bioplib-X.Y.tar.gz
        tar -xvf bioplib-X.Y.tar
-or- (if you have Gnu tar)
        tar -zxvf bioplib-X.Y.tar.gz

(where X.Y is the major and minor version numbers - e.g. 3.0)

####(3) This will create a directory called bioplib-X.Y.  

Enter this directory and then go into the src sub-directory:

        cd bioplib-X.Y/src

Modify the Makefile as required for your system. If you are using the
GNU C compiler and have followed the directions above, no changes
should be needed. If you have chosen alternative locations for the
include and library directories then you will need to change LIBDEST
to be the library directory you have created and INCDEST to the
include directory you have created (N.B. the complete path is
required, you can't do ~/lib).

####(4) Type the commands:

        make 
        make doxygen
        make install

####(5) For more information, read the file:

        bioplib-X.Y/doc/doxygen/docsrcinput/page_01.dox

or, after doing 'make doxygen', read the formatted version by pointing
a web browser at:

        bioplib-X.Y/doc/html/index.html

