#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "time.h"
#include "string.h"
#include <iostream>
using namespace std;

#include "BCTool.h"
#include "Cell.h"
#include "CStemCell.h"
#include "System.h"
#include "Random.h"

extern struct IMD _MD;
extern CRandom Rand;
extern struct DrugSchedule _Schedules[MAX_Schedules]; 

CSystem::CSystem()
{
	_NumCell = 0;
}

CSystem::CSystem(int N0)
{
    _Prolif=1.0;
    
    _NumPoolCell = N0;      // Number of cells in the simulation pool
    _NumCell = N0;          // Number of total cell numbers
    _MaxNumCell = MAXCELL;
    _Nrest = N0;            // Assume the all the initial cells are in the resting phase
    _Nprol = 0;
    
    _Currdrug=0;  //
    
    _cells = new CStemCell[_MaxNumCell];
}

CSystem::~CSystem()
{
     delete _cells;
}

bool CSystem::Initialized()
{
	int k;
	k=1;
	do
	{
		if((*this)(k).Initialized(k, _MD.cellpar))
		{
			k++;
		}
		else
		{
		}
	}while(k<=_MaxNumCell);

	return(true);
}

bool CSystem::SystemUpdate(double t)
{
    int k;
    int Ntemp;
    double *X1, *X2, *X3;      // new variable X3
	//  X1 stores the state x1, X2 stores the state x2, and X3 stores x3 for all cells
    bool *PQ;     //  store the state of _ProlifQ for all cells
    double *age;    //  store the state of age for all cells
    double Nresttemp, Nproltemp;
    
    X1 = new double[2*_NumPoolCell+1];
    X2 = new double[2*_NumPoolCell+1];
    X3 = new double[2*_NumPoolCell+1];
    PQ = new bool[2*_NumPoolCell+1];
    age = new double[2*_NumPoolCell+1];
    Ntemp = 0;        // Number of cells after a cycle.
    Nresttemp = 0.0;
	Nproltemp = 0.0;
    
    _N0 = 0;          // Number of cells in the pool after cell fate decision.
    _N1 = 0;          // Number of cells removed from resting phase
    _N2 = 0;          // Number of cells undergoing mitosis
    _N3 = 0;          // Number of cells removed from the proliferating phase
    _N4 = 0;          // Number of cells remain unchanged at the proliferating phase
    _N5 = 0;          // Number of cells remain unchanged at the resting phase
    
    DrugSys(t);
    
    for (k=1; k<=_NumPoolCell; k++)
    {
        //nextcell=(*this)(k).CellFateDecision(_NumCell, _MD.dt, t, _Currdrug);
        nextcell=(*this)(k).CellFateDecision(_Nrest, _MD.dt, t, _Currdrug);
        switch(nextcell.type){
            case 3:                 // Cells with ternimal differentiation
                _N1++;
                break;
            case 0:                 // Cells remain unchanged
            case 4:                 // Enter the prolifertaing state
                Ntemp++;
                X1[Ntemp] = nextcell.X1[0];
                X2[Ntemp] = nextcell.X1[1];
                X3[Ntemp] = nextcell.X1[2];
                PQ[Ntemp] = nextcell.ProfQ;
                age[Ntemp] = nextcell.age;
                
                
                _N0++;
                if(nextcell.ProfQ==0)
                {
                    _N5++;
                    Nresttemp = Nresttemp + 1.0;
                }
                else
                {
                    _N4++;
                    Nproltemp = Nproltemp + 1.0;
                }
                
                break;
            case 1:                 // Mitosis
                Ntemp++;
                X1[Ntemp] = nextcell.X1[0];
                X2[Ntemp] = nextcell.X1[1];
                X3[Ntemp] = nextcell.X1[2];
                PQ[Ntemp] = nextcell.ProfQ;
                age[Ntemp] = nextcell.age;

                
                Ntemp++;
                X1[Ntemp] = nextcell.X2[0];
                X2[Ntemp] = nextcell.X2[1];
                X3[Ntemp] = nextcell.X2[2];
                PQ[Ntemp] = nextcell.ProfQ;
                age[Ntemp] = nextcell.age;


                _N0 = _N0+2;
                _N2 = _N2+1;
                Nresttemp = Nresttemp + 2.0;
                
                break;
            case 2:                 // Apoptosis during proliferation
                _N3++;
                break;
        }
        
    }
    
    double _NumCelltemp; // 为了下面求 kinetic rates 
    _NumCelltemp=_NumCell;
    
    _Prolif = 1.0*Ntemp/_NumPoolCell; 
    _NumCell =_NumCell * _Prolif;  // 通过一次迭代后模拟池的数量变化比例来得到整个细胞数量的变化比例 
    
    _Nrest =  _NumCell * (Nresttemp + 0.001)/(Nresttemp + Nproltemp + 0.001);
    _Nprol =  _NumCell - _Nrest;
    
    if(Ntemp==0)
    {
        return(0);
    }
    
    // obtain the average values of three epigenetic states
    EpiAver(X1, X2, X3, Ntemp);
    //obtain the number of different phenotype cell
    PhenotypeNum(X1, X2, X3, Ntemp);
    
    int i;
    double p0;

    p0 = 1.0*_MaxNumCell/Ntemp;
    k=0;
    for (i=1; i<=Ntemp; i++)
    {
        if(Rand()<p0 && k<_MaxNumCell)
        {
            k=k+1;
            (*this)(k)._X[0] = X1[i];
            (*this)(k)._X[1] = X2[i];
            (*this)(k)._X[2] = X3[i];
            (*this)(k)._ProfQ = PQ[i];
            (*this)(k)._age = age[i];
            
			(*this)(k)._Krs[0]=(*this)(k).fbeta(_Nrest, X1[i], X2[i], X3[i],_Currdrug);
            (*this)(k)._Krs[1]=(*this)(k).fkappa(X2[i]);
            (*this)(k)._Krs[2]=(*this)(k).fdeathrate(X1[i], X3[i],_Currdrug);
            (*this)(k)._Krs[3]=(*this)(k).ftau(X1[i], X3[i]);
        }
        if (k==_MaxNumCell)
        {
            break;
        }
    }
    _NumPoolCell = k;
    
    delete X1;
    delete X2;
    delete X3;
    delete PQ;
    delete age;
    return(1);
}

void CSystem::Run()
{
    double t;
    int step;
    int k;
    FILE *fmdsys;
    char fn[ComLength+5];
    sprintf(fn,"%s.dat",_MD.mdcrd);
    if((fmdsys = fopen(fn,"w"))==NULL)
    {
        cout<<"Cannot open the file mdsys."<<endl;
        exit(0);
    }
    
    step=0;
    //OutPutSys(step);

    for (t=0; t<=_MD.T1; t=t+_MD.dt)
    {
        step=step+1;
        for (k=1; k<=_NumPoolCell; k++){
            (*this)(k).Pre_OneStep();
        }
        if(SystemUpdate(t))
        {
            fprintf(fmdsys,"%.2f %.4f %.2f %.2f %.2f %d %d  %.2f %.2f %.2f %.2f %.2f %.2f %.2f  %.1f\n",\
			                t,_Prolif, _NumCell, _Nrest, _Nprol, _NumPoolCell, _N0, _PhenotypeNum[0],  _PhenotypeNum[1],\
							_PhenotypeNum[2], _PhenotypeNum[3],_EpiAver[0], _EpiAver[1], _EpiAver[2],_Currdrug);
			// stop output the cell state at discrete time points for saving memory
            if((_MD.ntpx>0) && (step%_MD.ntpx==0))
            {
                OutPutSys(step);
            }
        }
        else{
        	fprintf(fmdsys,"%.2f %.4f %.2f  %.2f %.2f %d %d  %.2f %.2f %.2f %.2f %.2f %.2f %.2f  %.1f\n",\
			                t,_Prolif, _NumCell,  _Nrest, _Nprol, _NumPoolCell, _N0, _PhenotypeNum[0],  _PhenotypeNum[1],\
							_PhenotypeNum[2], _PhenotypeNum[3],_EpiAver[0], _EpiAver[1], _EpiAver[2],_Currdrug);
            printf("Tumor cells cleaned at day %4.2f\n",t);
            break;
        }
        
    }
    fclose(fmdsys);
}

void CSystem::RunOneStep(double t)
{
}

void CSystem::OutPutSys(int step)
{
    FILE *fp;
    int k;
    char fnc[StrLength];
    
    
    //(*this)(1).OutPutParameters();
    
    sprintf(fnc,"%s-%d.dat",_MD.mdcrd,step);
    if((fp = fopen(fnc,"w"))==NULL)
    {
        cout<<"Cannot open the file fmdc."<<endl;
        exit(0);
    }

	
	for(k = 1; k<= _NumPoolCell; k++)
	{
	// Edit to change you output items.

        fprintf(fp,"%d  %6.3f  %6.3f  %6.3f   %6.3f  %6.3f  %6.3f  %6.3f \n",\
		         k, (*this)(k)._X[0], (*this)(k)._X[1], (*this)(k)._X[2],  (*this)(k)._Krs[0],\
				  (*this)(k)._Krs[1],(*this)(k)._Krs[2], (*this)(k)._Krs[3]);
	}
	 
    fclose(fp);
}

void CSystem::OutPutCells(double t)
{
	int k;
	for(k = 1; k<= _NumPoolCell; k++)
	{
        (*this)(k).OutPut(t);
	}
}


CStemCell& CSystem::operator()(int indx)
{
if(indx>0 && indx<=_MaxNumCell)
{
	return *(_cells+(indx-1));
}
else
{
	cout<<"Err<< CSystem () >>Dimensional error"<<endl;
	exit(0);
}
}

void CSystem::EpiAver(double X1[], double X2[], double X3[], int Num)
{
	if (Num<=0)
	{
		cout << "Cell number is zero, EpiAver cannot be obtained" << endl;
	}
	for (int i=0;i<3;i++)
	{
	    _EpiAver[i]=0.0;
	}

	for(int k=0; k<Num; k++)
    {
        _EpiAver[0]=_EpiAver[0]+X1[k];
        _EpiAver[1]=_EpiAver[1]+X2[k];
        _EpiAver[2]=_EpiAver[2]+X3[k];
	}
	
    _EpiAver[0]=_EpiAver[0]/Num;
    _EpiAver[1]=_EpiAver[1]/Num;
    _EpiAver[2]=_EpiAver[2]/Num;
	
}

void CSystem::PhenotypeNum(double X1[], double X2[], double X3[], int Num)
{   
    double Splitpoints[3]={2.5,2.5,2.5}; 
	if (Num<=0)
	{
		cout << "Cell number is zero, PhenotypeNum cannot be obtained" << endl;
	}
	for (int i=0;i<4;i++)
	{
	    _PhenotypeNum[i]=0.0;
	}
	
	for(int k=0; k<Num; k++)
    {
       if ( (X1[k]>Splitpoints[0]) && (X2[k]<=Splitpoints[1]) && (X3[k]<=Splitpoints[2]) )
       {  
           _PhenotypeNum[0]+=1;
           continue;
	   }
	    if ( (X1[k]<=Splitpoints[0]) && (X2[k]>Splitpoints[1]) && (X3[k]>Splitpoints[2]) )
       {  
           _PhenotypeNum[1]+=1;
           continue;
	   }
	    if ( (X1[k]<=Splitpoints[0]) && (X2[k]<=Splitpoints[1]) && (X3[k]>Splitpoints[2]) )
       {  
           _PhenotypeNum[2]+=1;
           continue;
	   }
	   
	   _PhenotypeNum[3]+=1;
	   
	}
	
	for (int i=0;i<4;i++)
	{
	    _PhenotypeNum[i]=_NumCell*_PhenotypeNum[i]/Num;
	}

}

void CSystem::DrugSys(double t)
{
	
	_Currdrug=0.0; 
	
	int tempMaxcount;
	int * drugtimes; 
	
	if (_MD.schedule!=_Schedules[_MD.schedule].DS_id)
	{
		printf("Check the schedule ID in drug dat and MD.in");
	}
	
	tempMaxcount = _Schedules[_MD.schedule].DS_count;
	drugtimes = _Schedules[_MD.schedule].DS_times;
	if (t==0)
		printf("drug schedule %d:%s\n",_Schedules[_MD.schedule].DS_id,_Schedules[_MD.schedule].DS_name);
    switch (_MD.schedule)
    {
    	case 0:
    		{
    			_Currdrug = 0.0;
    			break;
			}
		case 1:
		case 2:
		case 3:
		case 4:
			{
				for (int i=1;i<=tempMaxcount;i++)
				{
					if (t>=drugtimes[2*i-2] && t<=drugtimes[2*i-1])
						_Currdrug = 1.0; 
				}
				break;
			}

        default:
            printf("Schedule Wrong\n");
            break;

   }
}
