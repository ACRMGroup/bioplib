/************************************************************************/
/* Defines and macros
*/
#define MAXTAGNAME 360  /* Maximum length of an XML tag name            */

typedef struct 
{
   REAL (*realFunction)(PDB *);
   int (*intFunction)(PDB *);
   char *(*stringFunction)(PDB *);
   char tag[MAXTAGNAME];
   int  type;
}  PDBTAGVAR;

#define PDBTAGVAR_REAL   0
#define PDBTAGVAR_INT    1
#define PDBTAGVAR_STRING 2

/* Macro to initialize the binding between a tag name and a function that
   can extract a variable from a PDB structure. The function takes a PDB
   pointer (PDB *) as input and returns REAL, int or (char *) as specified
   by the type (PDBTAGVAR_REAL, PDBTAGVAR_INT or PDBTAGVAR_STRING). The
   taglabel is a string specifying an XML tag that should be associated 
   with the value

   Called as:
   INIT_PDBTAGVAR(&myFunction, type, taglabel)
*/
#define INIT_PDBTAGVAR(f, t, l)                                                  \
do {                                                                             \
   if(gPDBTagFunctions == NULL)                                                  \
   {                                                                             \
      gPDBTagFunctions=(PDBTAGVAR *)malloc(sizeof(PDBTAGVAR));                   \
   }                                                                             \
   else                                                                          \
   {                                                                             \
      gPDBTagFunctions=(PDBTAGVAR *)realloc(gPDBTagFunctions,                    \
                    (1+gNPDBTagFunctions)*sizeof(PDBTAGVAR));                    \
   }                                                                             \
   if(gPDBTagFunctions==NULL) break;                                             \
   strncpy(gPDBTagFunctions[gNPDBTagFunctions].tag, l, MAXTAGNAME);              \
   gPDBTagFunctions[gNPDBTagFunctions].realFunction = NULL;                      \
   gPDBTagFunctions[gNPDBTagFunctions].intFunction = NULL;                       \
   gPDBTagFunctions[gNPDBTagFunctions].stringFunction = NULL;                    \
   switch(t)                                                                     \
   {                                                                             \
   case PDBTAGVAR_REAL:                                                          \
      gPDBTagFunctions[gNPDBTagFunctions].realFunction = (REAL (*)(PDB *))f;     \
      break;                                                                     \
   case PDBTAGVAR_INT:                                                           \
      gPDBTagFunctions[gNPDBTagFunctions].intFunction = (int (*)(PDB *))f;       \
      break;                                                                     \
   case PDBTAGVAR_STRING:                                                        \
      gPDBTagFunctions[gNPDBTagFunctions].stringFunction = (char * (*)(PDB *))f; \
      break;                                                                     \
   }                                                                             \
   gPDBTagFunctions[gNPDBTagFunctions].type = t;                                 \
   gNPDBTagFunctions++;                                                          \
} while(0)


/************************************************************************/
/* Globals
*/


/************************************************************************/
/* Prototypes
*/
void blPrintTagVariables(PDB *p);
void blPrintAllTagVariables(PDB *pdb);
#ifdef _PDBTAGVARS_CODE
PDBTAGVAR *gPDBTagFunctions = NULL;
int       gNPDBTagFunctions = 0;
#else
extern PDBTAGVAR *gPDBTagFunctions;
extern int       gNPDBTagFunctions;
#endif
