#ifndef CSYSTEM_H
#define CSYSTEM_H

#include "CStemCell.h"

class CSystem{
private:
    int _NumPoolCell;   // Number of cells in the simulation pool.
    double _NumCell;		// Number of total cells in the animal (in unit 10^6).
    int _MaxNumCell;    // Maximal cell number.
    double _Nrest, _Nprol;  // Number of cells at the resting phase, and proliferating phase. 
    
    double _Prolif;       // Proliferation rate
    int _N0, _N2, _N1, _N3, _N4, _N5;
    
    // double _X1Aver, _X2Aver, _X3Aver;  // the average of three epigenetic states
    double _EpiAver[3]; // the average of three epigenetic states
    double _PhenotypeNum[4]; // the number of phenotype cell, 0: DSCs; 1: DTP; 2: DSCs; 3: Others
    double _Currdrug;  // current drugdose
    
	CStemCell *_cells;	// Cells.

	void OutPutCells(double t);
	void OutPutSys(int step);
	void OutPutSys(FILE *fp);
    void RunOneStep(double t);
    void EpiAver(double X1[], double X2[], double X3[], int Num); //obtain the average of 3 epigenetic states
    void PhenotypeNum(double X1[], double X2[], double X3[], int Num); 
	//obtain the number of different phenotype cell
	void DrugSys(double t);
    
    DCell nextcell;
	
public:
	CSystem();
	~CSystem();
	
	CSystem(int N0);

	CStemCell& operator()(int indx); // Edit to return the type of your own class.

	bool Initialized();
	void Run();
    
    bool SystemUpdate(double t);
};

#endif
