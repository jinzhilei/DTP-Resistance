#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "time.h"
#include "string.h"
#include <iostream>
using namespace std;

#include "Random.h"
#include "BCTool.h"
#include "CStemCell.h"
#include "Cell.h"

extern struct IMD _MD;
extern CRandom Rand;
extern double _ratio19;


CStemCell::CStemCell()
{
}

CStemCell::~CStemCell()
{
}

bool CStemCell::Initialized(int k, char fpar[])
{
    SetDefaultPar();   // blank function--lxy
	ReadPar(fpar);
	_t = 0.0;
	_cellid=k;
	_NumVar=3;
	_NumKrs=4;
	_Krs=new double[_NumKrs];
	_X = new double[_NumVar];
    _X[0] = Rand(_par.x1-0.5,_par.x1+0.5);        // Initial state of the epigenetic state x1
    _X[1] = Rand(_par.x2-0.5,_par.x2+0.5);      // Initial state of the epigenetic state x2
    _X[2] = Rand(_par.x3-0.5,_par.x3+0.5);     // Initial state of the epigenetic state x3
    _ProfQ = 0;               // We assume that all cells at resting phase at the initial state.
    _age = 0;
    
    return true;
}

void CStemCell::ReadPar(char fpar[])
{
	FILE *fp;
	char str[StrLength], *pst;
	if((fp = fopen(fpar,"r"))==NULL)
	{
		cout<<"Cannot open the cell parameter input file."<<endl;
		exit(0);
	}
	rewind(fp);
	while(!feof(fp))
	{
		fgets(str,StrLength,fp);
		if(str[0]=='#'){ continue;}

        if((pst=strstr(str,"beta0="))!=NULL)
        {
            _par.beta0=atof(pst+6);
        }
        if((pst=strstr(str,"kappa0="))!=NULL)
        {
            _par.kappa0=atof(pst+7);
        }
        if((pst=strstr(str,"p_a1="))!=NULL)
        {
            _par.a1=atof(pst+5);
        }
        if((pst=strstr(str,"p_a2="))!=NULL)
        {
            _par.a2=atof(pst+5);
        }
        if((pst=strstr(str,"p_a3="))!=NULL)
        {
            _par.a3=atof(pst+5);
        }
        if((pst=strstr(str,"p_a4="))!=NULL)
        {
            _par.a4=atof(pst+5);
        }
        if((pst=strstr(str,"p_a5="))!=NULL)
        {
            _par.a5=atof(pst+5);
        }
        if((pst=strstr(str,"p_b1="))!=NULL)
        {
            _par.b1=atof(pst+5);
        }
        if((pst=strstr(str,"alpha="))!=NULL)
        {
            _par.alpha=atof(pst+6);
        }
        if((pst=strstr(str,"mu0="))!=NULL)
        {
            _par.mu0=atof(pst+4);
        }
         if((pst=strstr(str,"mu1="))!=NULL)
        {
            _par.mu1=atof(pst+4);
        }
        if((pst=strstr(str,"tau0="))!=NULL)
        {
            _par.tau0=atof(pst+5);
        }
        if((pst=strstr(str,"theta0="))!=NULL)
        {
            _par.theta0=atof(pst+7);
        }
        if((pst=strstr(str,"p_gamma0="))!=NULL)
        {
            _par.gamma0=atof(pst+9);
        }
        if((pst=strstr(str,"a="))!=NULL)
        {
            _par.a=atof(pst+2);
        }
        if((pst=strstr(str,"x1="))!=NULL)
        {
            _par.x1=atof(pst+3);
        }
        if((pst=strstr(str,"x2="))!=NULL)
        {
            _par.x2=atof(pst+3);
        }
        if((pst=strstr(str,"x3="))!=NULL)
        {
            _par.x3=atof(pst+3);
        }
		if((pst=strstr(str,"d="))!=NULL)
        {
            _par.d=atof(pst+2);
        }
 	}
	fclose(fp);
    
//    OutPutParameters();    
}

void CStemCell::OutPutParameters()
{
    printf("beta0=%f, theta0=%f, alpha=%f\n", _par.beta0, _par.theta0, _par.alpha);
    printf("kappa0=%f, mu0=%f, mu1=%f, tau0=%f\n",_par.kappa0, _par.mu0, _par.mu1, _par.tau0);
    printf("a1=%f, a2=%f, a3=%f, a4=%f, b1=%f, gamma0=%f\n", _par.a1, _par.a2, _par.a3, _par.a4, _par.b1, _par.gamma0);
    printf("a=%f, d=%f\n",_par.a, _par.d);
}

void CStemCell::SetDefaultPar()
{
}

void CStemCell::GetRandForces(double b[][NUMRAND], double x[], int k)
{
}

void CStemCell::GetTrends(double a[], double x[])
{
}

DCell CStemCell::CellFateDecision(long N0, double dt, double t)
{
    int i;
    double mu;
    double beta;
    double kappa;
    double tau;  //a new variable since each cell has its unique duration of cell cycle.
    double drugdose; // new variable
    double _rand;
    
    
    mu = fdeathrate(_X[0],_X[2]) * dt;
    beta = fbeta(N0, _X[0], _X[1], _X[2]) * dt;
    kappa = fkappa(_X[1]) * dt;
    tau = ftau(_X[0], _X[2]);
    drugdose = fdrug(t);
    
    DCell nextcell;
    
    nextcell.type = 0;  //Default, cells remain unchanged.
    for (i=0;i<_NumVar;i++)
    {
        nextcell.X1[i] = _X[i];
    }
    nextcell.ProfQ = _ProfQ;
    nextcell.age = _age;
    
    if (_ProfQ == 0)        // If the cell is at the resting phase
    {
        _rand = Rand();
        if (_rand < kappa)  // Terminal differentiation from the resting phase
        {
            nextcell.type = 3;
        }
        else
        {
            if (_rand < kappa + beta)   // Enter the proliferting phase
            {
                nextcell.type = 4;
                for (i=0;i<_NumVar;i++)
                {
                    nextcell.X1[i] = _X[i];
                }
                nextcell.ProfQ = 1;
                nextcell.age = 0;
            }
        }
    }
    
    if (_ProfQ == 1)    // If the cell is at the proliferating phase
    {
        _rand = Rand();

        if (_rand < mu)     // The apoptosis during the proliferating phase
        {
            nextcell.type = 2;
        }
        else
        {
            if (_age < tau)
            {
                nextcell.age = _age + dt;
            }
            else  // Perform mitosis
            {
                nextcell.type = 1;
                for (i=0;i<_NumVar;i++)
                {
                    nextcell.X1[i] = GetnextEpi(i, t, _X[0], _X[1], _X[2], drugdose);
                    nextcell.X2[i] = GetnextEpi(i, t, _X[0], _X[1], _X[2], drugdose);
                }

                nextcell.ProfQ = 0;
                nextcell.age = 0;
            }
        }
    }
    
    return(nextcell);
}

double CStemCell::GetnextEpi(int i,double t, double x1, double x2, double x3, double drug)
{
    double phi;
    double a, b;
    double z;
    double K1, K2, K3;
    double n1,n2,n3;
    switch(i)
    {
        case 0:  // p(x1; y) 
        	n1=2.0;
        	K1 = 4.9- 4.2 * exp(-0.32 * x3);
        	//phi = 0.7 + 4.9 * pow(x1, n1)/(pow(K1, n1)+pow(x1, n1));
            phi = 0.6 + (5.0-2.5*drug) * pow(x1, n1)/(pow(K1, n1)+pow(x1, n1))-0.15*drug;
            break;
        case 1:
            n2=2.2;
			K2 = 3.1;
            phi = 0.65 + 5.0 * pow(x2, n2)/(pow(K2, n2)+pow(x2, n2));
//            n2=2.0;
//			K2 = 3.2;
//            phi = 0.55 + 5.0 * pow(x2, n2)/(pow(K2, n2)+pow(x2, n2));
            break;
        case 2:
        	n3=2.0;
        	K3 = 4.9- 4.0* exp(-0.25 * x1);
			phi =  0.6 + 5.0 * pow(x3, n3)/(pow(K3, n3)+pow(x3, n3));
        	break;
    }
    a = _par.a;
    b = phi/a;
    z = Rand.GammaDistribution(a, b);
    return(z);
}

void CStemCell::Update(int mu)
{
}

void CStemCell::Propensities(double a[])
{
}

double CStemCell::fbeta(long N0, double x1, double x2, double x3)
{
    double beta, beta0;
    double theta;
    double n_beta;
    n_beta=4.0;
    beta0 = _par.beta0;
    theta = _par.theta0 * (1.0 + _par.a5 * pow(_par.a4 * (x1+_par.alpha *x3), 6.0)/(1.0 + pow(_par.a4 * (x1+_par.alpha *x3), 6.0)));
	beta=beta0 * (1.0/(1.0 + N0/theta)) * ((_par.a1 * x2 + pow(_par.a2 * x2,n_beta))/(1+ pow(_par.a3 * x2,n_beta)));  
	 // Proliferation rate of cells
    
	return(beta);
}

double CStemCell::fkappa(double x2)
{
    double kappa;
    
    kappa = _par.kappa0 * 1.0/(1.0 + pow(_par.b1 * x2, 4.0));

    return(kappa);
}

double CStemCell::fdeathrate(double x1, double x3)
{
    double mu;
    
    mu = _par.mu0 * exp(_par.mu1 * (ftau(x1,x3)-_par.tau0) );
    
    return(mu);
}
double CStemCell::ftau(double x1, double x3)  // new function
{
	double tau,tau0,gamma1;
	
	tau0=_par.tau0;
	gamma1=4.0 + _par.alpha;
	
	tau=tau0 * exp(_par.gamma0 * ( gamma1 - (x1 + _par.alpha * x3) ) );
	
	return(tau);
}                   
double CStemCell::fdrug(double t)    // new function
{
	double drugdose;
	drugdose=0.0;
	
	if(t>720.0) //t>30 day
	{
	    drugdose = _par.d;
	}
	
	return(drugdose);
}

void CStemCell::Pre_OneStep()
{
}



