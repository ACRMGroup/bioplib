/* Include file for tests */

#ifndef _WHOLEPDB_H
#define _WHOLEPDB_H

/* Includes for tests */
#include <stdlib.h>
#include <check.h>
#include <unistd.h>

/* Includes from source file */
#include "../../port.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "../../macros.h"
#include "../../general.h"
#include "../../pdb.h"
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <time.h>

/* Prototypes */
Suite *wholepdb_suite(void);

#endif
