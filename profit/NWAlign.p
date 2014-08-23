void NWAlign(int strucnum)
;
void ReadAlignment(char *alnfile)
;
BOOL RemoveDoubleDeletions(char *seqa, char *seqb)
;
void SetNWZones(char *ref_align, char *mob_align, int align_len,
                PDB **RefIndex, PDB **MobIndex, int strucnum)
;
void MergeZones(int strucnum)
;
BOOL VerifySequence(char *seqa, char *seqb)
;
int AlignmentWrapper(int strucnum, char *command, BOOL append)
;
int AlignmentFromZones(char *filename, BOOL fasta)
;
int TruncateSeq(char *outstring, char *instring, int start, int stop)
;
int AlignmentFromZones_PIR(char *filename)
;
int ChainBreakToGap(char *instring)
;
ZONE *AlignToZone(char *ref_align, char *mob_align, 
		  int   ref_start,  int   mob_start)
;
BOOL CommonZones(void)
;
BOOL BuildMultiAlignment(char **seqs, char **alns)
;
void PrintSequencePIR(FILE *fp, char *sequence, int width)
;

