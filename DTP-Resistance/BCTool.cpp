#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "time.h"
#include "string.h"
#include <iostream>
using namespace std;

#include "BCTool.h"
#include "CStemCell.h"
#include "System.h"
#include "Random.h"

#include "Random.cpp"
#include "CStemCell.cpp"
#include "System.cpp"
#include "Cell.cpp"

struct IMD _MD;
CRandom Rand;
struct DrugSchedule _Schedules[MAX_Schedules]; 

void ReadIPF(char *fn);
void SetParVal(char *str, char const *conststr, char val[]);
void help();
void SetDefault();
void OutputParameter();
int Read_DrugSchedule(const char *filename, DrugSchedule schedules[]);
void Display_DrugSchedule(const DrugSchedule *schedules);

int main(int argc, char *argv[])
{
	CSystem *sys;
	
	clock_t start,finish;
    double totaltime;
    start=clock();
	
  	if (argc<2)
  	{
  		help();
		exit(0);
  	}
	
	ReadIPF(argv[1]);
	if(_MD.seed > 0) Rand.Initialized(_MD.seed);
	
	//OutputParameter();
	
    sys = new CSystem(_MD.N0);
    if(sys->Initialized())
    {
        sys->Run();
        
        
        finish=clock();
        totaltime=(double)(finish-start)/CLOCKS_PER_SEC;
        cout<<"\n The running time is "<<totaltime<<" seconds！"<<endl;
        
        
        return(1);
    }
    else
    {
        return(0);
    }
    
    

    
    
}

void help()
{
	cout<<"Usage: "<<endl;
	cout<<"bct_StemCell inputfile"<<endl;
	cout<<"inputfile: The name of input file."<<endl;
	cout<<" For example: bct md.in."<<endl;
}

void ReadIPF(char *fn)
{
	FILE *fp;
	char str[StrLength], *pst;
	if((fp = fopen(fn,"r"))==NULL)
	{
		cout<<"Cannot open the input file."<<endl;
		exit(0);
	}
	SetDefault();
	rewind(fp);
	while(!feof(fp))
	{
		fgets(str,StrLength,fp);
		if(str[0]=='#'){ continue;}
		SetParVal(str, "mdcrd=\"", _MD.mdcrd);
		SetParVal(str, "cellpar=\"",_MD.cellpar);
		SetParVal(str, "drugpar=\"",_MD.drugpar);
        if((pst=strstr(str,"dt="))!=NULL)
		{
			_MD.dt=atof(pst+3);
		}
		if((pst=strstr(str,"T1="))!=NULL)
		{
			_MD.T1=atof(pst+3);
		}
        if((pst=strstr(str,"T0="))!=NULL)
        {
            _MD.T0=atof(pst+3);
        }
        if((pst=strstr(str,"ntpx="))!=NULL)
		{
            _MD.ntpx=atoi(pst+5);
        }
        if((pst=strstr(str,"ntpr="))!=NULL)
		{
            _MD.ntpr=atoi(pst+5);
        }
		if((pst=strstr(str,"seed="))!=NULL)
		{
			_MD.seed=atoi(pst+5);
		}
        if((pst=strstr(str,"N0="))!=NULL)
        {
            _MD.N0=atoi(pst+3);
        }
        if((pst=strstr(str,"schedule="))!=NULL)
        {
            _MD.schedule=atoi(pst+9);
        }
//        if((pst=strstr(str,"Tstart="))!=NULL)
//        {
//            _MD.Tstart=atof(pst+7);
//        }
//        if((pst=strstr(str,"Tend="))!=NULL)
//        {
//            _MD.Tend=atof(pst+5);
//        }
	}
	fclose(fp);
	
	memset(_Schedules, 0, sizeof(_Schedules));
	int DS_count = Read_DrugSchedule(_MD.drugpar, _Schedules);
	if (DS_count < 0) {
        fprintf(stderr, "Drug schedule file is empty, exit\n");
        exit(0);
    }
	
}

void SetParVal(char *str, char const *conststr, char val[])
{
	char *pst;
	if((pst=strstr(str,conststr))!=NULL)
	{
		strcpy(val,pst+strlen(conststr));
		if((pst = strstr(val,"\""))!=NULL)
		{
			val[pst - val] = '\0';
		}			
	}
	return;
}

void OutputParameter()
{
     cout<<"mdcrd="<<_MD.mdcrd<<endl;
	 cout<<"N0="<<_MD.N0<<endl;
	 cout<<"dt="<<_MD.dt<<endl;
     cout<<"T1="<<_MD.T1<<endl;
	 cout<<"ntpr="<<_MD.ntpr<<endl;
	 cout<<"seed="<<_MD.seed<<endl;
}

void SetDefault()
{
     _MD.ntpx=0;
     _MD.ntpr=0;
     _MD.N0=1;
	 _MD.seed=0;
}

int Read_DrugSchedule(const char *filename, DrugSchedule schedules[]) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        perror("Cannot open drug schedule file");
        return -1;
    }
	
    char line[MAX_LINE];
    int schedule_index = 0;
    int line_num = 0;

    // 逐行读取文件
    while (fgets(line, MAX_LINE, file) != NULL && schedule_index <= MAX_Schedules) {
        line_num++;
        // skip the blank line and comment line (begin with #)
        if (line[0] == '\n' || line[0] == '#') {
            continue;
        }

        // Drug schedule： Schedule ID,name, drug on time 1,drug off time 1,...
        // Example：1,control,1,10,20
        char *token = strtok(line, ",");
        if (token == NULL) {
            fprintf(stderr, "Error in the %d line, skip\n", line_num);
            continue;
        }

        schedules[schedule_index].DS_id = atoi(token);
        
        token = strtok(NULL, ",");
        if (token == NULL) {
            fprintf(stderr, "Miss schedule name in line %d，skip\n", line_num);
            continue;
        }
        int tempnamelen = 50;
        strncpy(schedules[schedule_index].DS_name, token, tempnamelen - 1);
        schedules[schedule_index].DS_name[tempnamelen-1] = '\0';
        //printf("schedules name：%s\n", schedules[schedule_index].DS_name);
        
        token = strtok(NULL,",");
        if (token == NULL) {
            fprintf(stderr, "Miss the number of drug administration in line %d，skip\n", line_num);
            continue;
        }
        schedules[schedule_index].DS_count = atoi(token);
        
		// drug on and off time points
        int tempMaxcount = 0, tempcount = 0;
		tempMaxcount = 2 * schedules[schedule_index].DS_count;
        schedules[schedule_index].DS_times = new int[tempMaxcount];
        while ((token = strtok(NULL, ",")) != NULL && tempcount < tempMaxcount) {	
            schedules[schedule_index].DS_times[tempcount] = atoi(token);
            tempcount = tempcount + 1;
        }
        
        schedule_index++;
        
    }

    fclose(file);
    
    if (schedule_index == 0)
	schedule_index = -1;
    
    return schedule_index;
}

void Display_DrugSchedule(const DrugSchedule *schedules) {
    printf("\n========== Drug Schedule %d ==========\n", schedules->DS_id);
    printf("schedules name: %s\n", schedules->DS_name);
    printf("drug on and off time points：");
    for (int i = 0; i < 2*schedules->DS_count; i++) {
        printf("%dh ", schedules->DS_times[i]);
    }

}
