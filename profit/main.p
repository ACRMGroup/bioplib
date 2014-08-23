void logo(void)
;
int main(int  argc, char **argv)
;
void SetRMSAtoms(char *command)
;
void SetRMSZone(char *command)
;
void Cleanup(void)
;
void Die(char *message)
;
BOOL ReadStructure(int  structure,
                   char *filename,
                   int  strucnum,
                   int  xmasFormat)
;
int InitParser(void)
;
int DoCommandLoop(FILE *script)
;
void SetFitAtoms(char *command)
;
void SetFitZone(char *command, int strucnum)
;
void DelFitZone(char *command, int strucnum)
;
void DelRMSZone(char *command, int strucnum)
;
void SetZoneStatus(char *status)
;
int GetResSpec(char *resspec, int *resnum, char *chain, char *insert)
;
int  ParseZone(char *zonespec,
               int  *start1, 
               int  *stop1,
               char *chain1,
               char *startinsert1,
               char *stopinsert1,
               int  *start2, 
               int  *stop2,
               char *chain2,
               char *startinsert2,
               char *stopinsert2,
               int  strucnum)
;
int FindSeq(char *zonespec,
            char *sequence, 
            int  *start, 
            int  *stop, 
            char *chain)
;
void ShowMatrix(void)
;
void ShowStatus(char *filename)
;
void stdprompt(char *string)
;
void FormatZone(char *zone, char chain, int start, char startinsert, 
                int stop,  char stopinsert)
;
void ReadMulti(char *filename, BOOL xmasFormat)
;
void WriteCoordinates(char *filename, int strucnum)
;
void WriteMulti(char *ext)
;
char *FindDash(char *buffer)
;
PDB *ReadXMAS(FILE *fp, int *natoms)
;
PDB *ReadXMASAtoms(FILE *fp, int *natoms)
;
PDB *doReadXMAS(FILE *fp, int *natoms, int readhet)
;
void SetCentreResidue(char *command)
;
int RunScript(char *command)
;
int ConvertResidueToSequential(ZONE *input_zone, int strucnum)
;
int CheckOverlap(ZONE *inputtest, ZONE *inputlist, int strucnum)
;
void SetZoneFromBValCol(void)
;
int ConvertSequentialToResidue(ZONE *input_zone, int strucnum)
;
ZONE *SortZoneList(ZONE *zonelist)
;
int ConvertZoneList(ZONE *zonelist, int strucnum, int mode)
;
int ConvertAllZones(int mode)
;
int SortAllZones(void)
;
ZONE *ChainList(PDB *pdblist)
;
BOOL SequentialZones(int strucnum)
;
BOOL SequentialZonesWholeSeq(int strucnum)
;
BOOL OneToOneChains(int strucnum)
;
int EnforceOneToOneChains(int strucnum)
;
BOOL *CheckListOverlap(ZONE *zonelist)
;
int CopyPDBListToRef(int strucnum)
;
int SetMobileToReference(int strucnum)
;
int fit_order_cmp(const void *ptr_scoreA, const void *ptr_scoreB)
;
int AllVsAllRMS(char *filename, BOOL print_tab, BOOL set_ref)
;
ZONE *SetOverlappingZones(ZONE *InputA, ZONE *InputB)
;
ZONE *RenumberZone(ZONE *InputZone, ZONE *OverlapZone)
;
int TrimZones(void)
;
int FitStructuresWrapper(void)
;
