PDBTAGVARS
==========

This is demonstration code for adding to Bioplib the ability to add
XML tags associated with fields of the PDB structure or from
PDB->extras. 

You simply need to create a function that takes a PDB pointer (to an
individual PDB record) and returns a REAL, int or (char *). You then
associate that with an XML tag name using:
   INIT_PDBTAGVAR(&functionName,  PDBTAGVAR_REAL,   "tagname");

Most of the work goes on in pdbtagvars.h and the idea is that the code
in pdbtagvars.c will be replaced with code that will add to the XML
structure rather than just print statements.

- pdbtagvars.h  Header file for using this code
- pdbtagvars.c  Needs to be enhanced for the real version
- main.c        Demo code

