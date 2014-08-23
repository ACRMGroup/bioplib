void FitStructures(void)
;
BOOL DoFitting(int NCoor, int strucnum)
;
int ValidAtom(char *atnam, int mode)
;
REAL CalcRMS(BOOL ByRes, FILE *fp, int strucnum, BOOL UpdateReference,
	     BOOL ByAtm)
;
void ShowNFitted(void)
;
REAL ShowRMS(BOOL ByRes, char *filename, int strucnum, 
             BOOL UpdateReference, BOOL ByAtm)
;
int CheckForConvergence(int NCoor, int strucnum)
;
int UpdateFitArrays(int strucnum)
;
REAL Distance(PDB *p, PDB *q)
;
REAL AlignOnCADistances(PDB **RefIndex, int length1, 
                        PDB **MobIndex, int length2,
                        char *align1, char *align2, int *align_len)
;
REAL TraceBackDistMat(REAL **matrix, 
                      XY   **dirn,
                      int  length1, 
                      int  length2, 
                      PDB  **RefIndex, 
                      PDB  **MobIndex, 
                      char *align1, 
                      char *align2, 
                      int  *align_len)
;
int SearchForBestDistMat(REAL **matrix, 
                         int  length1, 
                         int  length2, 
                         int  *BestI, 
                         int  *BestJ,
                         PDB  **RefIndex, 
                         PDB  **MobIndex, 
                         char *align1, 
                         char *align2)
;
int CreateFitArrays(int strucnum)
;
int CentreOnZone(int strucnum)
;
void SetSymmetricalAtomPAirs(void)
;
void ApplyMatrixCOOR(COOR *incoords,
		     REAL matrix[3][3],
		     int  ncoor)
;
void CalculateRotationMatrix(REAL RotAngle,
			     REAL Matrix[3][3])
;
REAL FitSingleStructure(int  strucnum,
			BOOL single_iteration)
;
void FitStructuresInOrder(REAL sortorder[][2])
;
void NoFitStructures(void)
;
BOOL DoNoFitting(int strucnum)
;
