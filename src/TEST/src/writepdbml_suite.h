/* Include file for tests */

#ifndef _WRITEPDBML_H
#define _WRITEPDBML_H

/* Includes for tests */
#include <stdlib.h>
#include <check.h>
#include <unistd.h>

/* Includes from source file */
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <libxml/tree.h>
#include <ctype.h>

#include "../../MathType.h"
#include "../../pdb.h"
#include "../../macros.h"

/* Prototypes */
Suite *writepdbml_suite(void);

#endif
