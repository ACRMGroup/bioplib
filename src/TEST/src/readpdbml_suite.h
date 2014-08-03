/* Include file for tests */

#ifndef _READPDBML_H
#define _READPDBML_H

/* Includes for tests */
#include <stdlib.h>
#include <check.h>

/* Includes from source file */
#include "../../port.h"    /* Required before stdio.h                   */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>
#include <libxml/parser.h>
#include <libxml/tree.h>

#include "../../SysDefs.h"
#include "../../MathType.h"
#include "../../pdb.h"
#include "../../macros.h"
#include "../../fsscanf.h"
#include "../../general.h"

/* Prototypes */
Suite *readpdbml_suite(void);

#endif
