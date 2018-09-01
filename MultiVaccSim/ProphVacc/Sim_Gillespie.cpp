
/*
Build a simple c++ script that will 
simulate a stochastic SIR epidemiological
model. The simulation procedes according
to the Gillespie algorithm, and simulates 
the use of vaccination to prevent the invasion
of a zoonotic pathogen. 

-This sim is TRUE Gillespie
*/

#include <iostream> // input/output std::cout, cin
#include <math.h>   // Math functions
#include <stdlib.h>
#include <time.h>
#include <fstream> //Used to open files to which data are saved
#include <vector> 
#include <string.h>
#include <sys/stat.h> // mkdir
#include <cmath> //fmod

//---------------------------------------------------------START HEADER FILE

//********
//CONSTANTS

const int NTrials = 25;
const int TPathLEN = 1;
const int IpInitLEN = 1;
const int tvLEN = 4;
const int BpLEN = 1;
const int NParSets = 1;
const int NumPars = 12;

//********
//USER-ASSIGNED VARIABLES
char SimName[50] = "Sim_Test";
bool VerboseWriteFlag = true;

std::vector<double> tvVals;

double TPathMIN = 5*365; double TPathMAX = 6*365; 
std::vector<double> TPathInvVals;

double bpvals[] = {0.00005};
std::vector<double> BpVals; 

int ipinitvals[]={100};
std::vector<int> IpInitVals; 

double TMax = 10.0*365.0; double tick = 1.0; //OneSim writes data at time-intervals tick

int SInit = 10000;

//********
//ARRAYS
double ParMat [NParSets][NumPars];
double TExtMat [NParSets][NTrials];

//********
//CRITICAL VARIABLES
int S, Iv, Ip, V, P, NPop, Par, IpInit;
char FileNamePar[50] = "Data/ParMat_"; 
char FileNameTExt[50] = "Data/TExtMat_";
char DirName[50] = "Data/";
char FileSuffix[50], FileNameDat[50];
double b0, tb, T, Nv,tv, b, gamv, R0p, Bp, gamp, d, Event_Rate, Event_Rate_Prod, RandDeath, dTime, t, ti, TPathInv;
std::ofstream out_data;

//OTHER VARIABLES
int nbirths, ndeaths, ninfv, ninfp, nrecv, nrecp, S_death, Iv_death, Ip_death, V_death, P_death;
int ntrial, whichindex;
double whichmin, unif, tmodT;
double Nudge = 0.000001;

//FUNCTION DECLARATIONS
void CheckEventConflict();
void ApplyEvent();
void GetTime();
void Initialize(); //Fills ParMat, initializes TExtMat, returns #Rows
void OneSim (double StartTime, double EndTime, bool StopOnErad);
double Rand () ;
std::vector<double> Seq(double, double, int); //Returns sequence
void VaccFun();
void WriteMat(double *arr, int NRow, int NCol, char* filename); 

//---------------------------------------------------------END HEADER FILE

//MAIN
int main()
{
  srand ( time(NULL) );

  strcat(FileNamePar, SimName);  
  strcat(FileNameTExt, SimName);  
  Initialize(); //Fill in the parameter matrix

  WriteMat((double *) ParMat, NParSets, NumPars, FileNamePar); //Write ParMat
  
  for(int Par = 0; Par < NParSets; Par++) //Loop through parameters
    {
      if(VerboseWriteFlag){
	strcat(DirName, SimName);  
	mkdir(DirName, ACCESSPERMS);
	strcat(DirName, "/");	  
	strcpy(FileNameDat, DirName);
	sprintf(FileSuffix, "Par_%d",Par);
	strcat(FileNameDat, FileSuffix);
	out_data.open(FileNameDat);
	out_data << "time S Iv Ip V P N births deaths ninfv ninfp nrecv nrecp S_death Iv_death Ip_death V_death P_death\n";
      }
      
      //Extract parameters
      //Parmat[0] corresponds to Par;
      b0 = ParMat[Par][1]; //b0
      d = ParMat[Par][2]; //d
      Bp = ParMat[Par][3]; //Bp
      Nv = ParMat[Par][4]; //Nv
      tv = ParMat[Par][5]; //tv
      gamv = ParMat[Par][6]; //gamv
      gamp = ParMat[Par][7]; //gamp
      tb = ParMat[Par][8]; //tb
      T = ParMat[Par][9]; //T
      IpInit = (int) ParMat[Par][10]; // Initial level of pathogen upon invasion
      TPathInv = ParMat[Par][11]; // Time of pathogen invasion
      
      for(ntrial = 0; ntrial < NTrials; ntrial++)
	{
	  S = SInit; //S
	  Iv = 0; //Iv
	  Ip = 0; //Ip
	  V = 0; //V
	  P = 0; //P
	  NPop = S + Iv + Ip + V + P;

	  //Simulate to quasi steady state (rewrites State)
	  OneSim(0.0, TPathInv, false);
	      
	  //Simulate invasion until TMax years, or pathogen extinction
	  Ip = IpInit;
	  OneSim(TPathInv, TMax, true);

	  
	  //Store final value of t in TExtMat
	  TExtMat[Par][ntrial] = t;
	      
	  std::cout << "Sim Trial: " << ntrial << std::endl;
	  
	}//End loop through NTrials

      if(VerboseWriteFlag)
	{
	  out_data.close();
	}
      
      std::cout << "*****************************" << std::endl;	  
      std::cout << "Finished Parameter Set " << Par+1 << " / " << NParSets << std::endl;
      std::cout << "*****************************" << std::endl;
      
    }//End Loop through NParSets

  WriteMat((double *)TExtMat, NParSets, NTrials, FileNameTExt); //Write TExtMat      
}//End Main

      
      
//////////////////////////////////////
///////// DEFINE FUNCTIONS ///////////
//////////////////////////////////////

//************************************
//Function that updates state variables
void ApplyEvent() {
  if(Event_Rate_Prod <= b) //Birth
    { S++; NPop++; nbirths++; } 
  else if(Event_Rate_Prod <= b + d*NPop) //Death
    {
      RandDeath = Rand()*NPop;
      if( RandDeath < S ) //S dies
	{S--; S_death++;}
      else if(RandDeath < S + Iv) //Iv dies
	{Iv--; Iv_death++;}
      else if(RandDeath < S + Iv + Ip) //Ip dies
	{Ip--; Ip_death++;}
      else if(RandDeath < S + Iv + Ip + V) //V dies
	{V--; V_death++;}
      else{P--; P_death++;}  //P dies	
      NPop--; ndeaths++;
    }  
  else if(Event_Rate_Prod <= b + d*NPop + Bp*Ip*S)    //Event: Pathogen infection of S
    {S--; Ip++; ninfp++;}
  else if(Event_Rate_Prod <= b + d*NPop + Bp*Ip*S + Bp*Ip*Iv) //Event: Pathogen infection of V
    {Iv--; Ip++; ninfp++;}
  else if(Event_Rate_Prod <= b + d*NPop + Bp*Ip*S + Bp*Ip*Iv + gamv*Iv) //Event: Iv Recovery
    {Iv--;V++;nrecv++;}
  else if(Event_Rate_Prod <= b + d*NPop + Bp*Ip*S + Bp*Ip*Iv + gamv*Iv + gamp*Ip) //Event: Ip Recovery
    {Ip--; P++; nrecp++;}
}

//*************************************
//function that finds the min of 3 doubles, and min's index CHECK
void CheckEventConflict (){
  whichmin = T - tmodT;
  whichindex = 0;
  if(tmodT < tb && tb-tmodT < whichmin)
    {
      whichmin = tb-tmodT;
      whichindex = 1;
    }
  if(tmodT < tv && tv-tmodT < whichmin)
    {
      whichmin = tv-tmodT;
      whichindex = 2;
    }
}

//************************************
//Function that returns the time for the next event
void GetTime (){
  // Select random time step (doesn't allow for infite time step)
  do
    {
      dTime = (double) rand() / RAND_MAX ; //Use var dTime for storage
      dTime = -log(dTime) / Event_Rate; //Assign next timestep
    } while (dTime == INFINITY); //Avoid infinite timesteps
}

//************************************
//function to Initialize values of 2D array
void Initialize()
{
  tvVals = Seq(1, 365, tvLEN);
  TPathInvVals = Seq(TPathMIN, TPathMAX, TPathLEN);
  BpVals.assign(bpvals, bpvals + BpLEN);
  IpInitVals.assign(ipinitvals, ipinitvals + IpInitLEN);
  
  //Fill in ParMat
  int i = 0;
  for(int i1=0; i1<tvVals.size(); i1++)
    for(int i2=0; i2<BpVals.size(); i2++)
      for(int i3=0; i3<TPathInvVals.size(); i3++)
	for(int i4=0; i4<IpInitVals.size(); i4++)
	  {
	    ParMat[i][0] = i; //Par
	    ParMat[i][1] = 400.0;   //b0
	    ParMat[i][2] = 0.004; //d
	    ParMat[i][3] = BpVals[i2]; //Bp
	    ParMat[i][4] = 5000.0; //Nv
	    ParMat[i][5] = tvVals[i1]; //tv
	    ParMat[i][6] = 0.007; //gamv
	    ParMat[i][7] = 0.007; //gamp
	    ParMat[i][8] = 90.0; //tb
	    ParMat[i][9] = 365.0; //T
	    ParMat[i][10] = (double) IpInitVals[i4]; //IpInit
	    ParMat[i][11] = TPathInvVals[i3]; //TPathInv
	    i++;
	  }
}

void OneSim (double StartTime, double EndTime, bool StopOnErad = false)
{
  //Set initial conditions: No Pathogen, no vaccination
  ti = StartTime;
  t = StartTime;

  tmodT = std::fmod(t,T);
  b = (double) b0*(tmodT < tb);

  //Track births, deaths, recoveries, etc for error checking
  nbirths = 0; ndeaths = 0; ninfv = 0; ninfp = 0; nrecv = 0; 
  nrecv = 0; S_death = 0; Iv_death = 0; Ip_death = 0;  
  V_death = 0; P_death = 0;

  while(t < EndTime && (Ip > 0 || !StopOnErad))
    {
      if( t >= ti && VerboseWriteFlag) //Write old values at intervals "tick"
	{
	  out_data << t << " " << S << " " << Iv << " " << Ip << " " << V << " " << P  << " " << 
	    NPop << " " << nbirths << " " << ndeaths <<  " " << ninfv << " " << ninfp << " " << nrecv << 
	    " " << nrecp << " " << S_death << " " << Iv_death << " " << Ip_death << " " << V_death << " " << P_death << "\n"; 
	  ti += tick;
	  nbirths = 0; ndeaths = 0;
	}  
      Event_Rate = b + d*NPop + Bp*S*Ip + Bp*Iv*Ip + gamv*Iv + gamp*Ip;
      Event_Rate_Prod = Event_Rate*Rand();
      
      GetTime(); //Get time to next event

      CheckEventConflict(); //Finds whichmin, the time of the next conflict
      if(dTime < whichmin) //If no conflicts, proceed with Gillepsie event
	{
	  ApplyEvent();	
	}
      else{ //If conflict, stop at conflict and perform necessary action
	
	dTime = whichmin;
	switch(whichindex)
	  {
	  case 0 : b = b0; break;//Start of birthing season
	  case 1 : b = 0.0; break;//End of birthing season
	  case 2 : VaccFun(); break;//Update S,Iv due to vaccination 
	  }
	if(dTime==0.0)
	  {dTime+=Nudge;}
      }    
      t += dTime;
      tmodT = std::fmod(t,T);
    }//End While
  
  if(VerboseWriteFlag)
    {
      out_data << t << " " << S << " " << Iv << " " << Ip << " " << V << " " << P  << " " << NPop << " " << nbirths << " " << ndeaths <<  " " << ninfv << " " << ninfp << " " << nrecv << " " << nrecp << " " << S_death << " " << Iv_death << " " << Ip_death << " " << V_death << " " << P_death << "\n"; 
    }
}

//************************************
//Function returning a random number between 0 and 1
double Rand (){
  do
    {
      unif = (double) rand() / RAND_MAX ;
    } while(unif==0.0);
  return unif;
}  

//************************************
//function that returns a sequence from minval to maxval, and of length lengthval
std::vector<double> Seq(double minval, double maxval, int lengthval) {
  std::vector<double> vec(lengthval);
  if(lengthval==1)
    {
      vec[0] = minval;
    }else{
  for(int seqindex = 0; seqindex < lengthval; seqindex++)
    {
      vec[seqindex] = minval + (double) seqindex/(lengthval-1)*(maxval-minval);
    }
  }
  return vec;
}

//************************************
//Vaccination function
void VaccFun(){
  int temp = (int) round( (double) Nv*S/(std::max(NPop,1)));
  temp  = std::min( temp,  S);
  S -= temp;
  Iv += temp;
}

//************************************
//function to write values of 2D array
void WriteMat(double *arr, int NRow, int NCol, char*filename)
{
  out_data.open(filename);
  int i,j;
  for(i=0; i<NRow;i++)
    {
      for(j=0; j<NCol; j++)
	{
	  //	  out_data << arr[i][j] << " "; 
	  out_data << *((arr+i*NRow) + j) << " "; 
	}
      out_data << std::endl;
    }
  out_data.close();
}
