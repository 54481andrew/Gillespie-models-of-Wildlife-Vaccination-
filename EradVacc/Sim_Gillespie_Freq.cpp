
/*
Build a simple c++ script that will 
simulate a stochastic SIR epidemiological
model. The simulation procedes according
to the Gillespie algorithm, and simulates 
the use of vaccination to prevent the invasion
of a zoonotic pathogen. 
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

//***********
//CONSTANTS
//***********
const int NTrials = 500;
const int TVaccLEN = 1; //TVacc is the year in which vaccination begins
const int IpInitLEN = 1; int ipinitvals[]={100};
const int tvLEN = 52; //double tvvals[] = {90.0};
const int tbLEN = 2; double tbvals[] = {60.0, 90.0};
//const int BpLEN = 10; double bpvals[] = {0.0105, 0.0140, 0.0246, 0.0280};
const int R0pLEN = 4; double r0pvals[] = {1.1, 1.5, 2, 3, 5};
//const int NvLEN = 3; int nvvals[] = {};
const int RhoLEN = 3; double rhovals[] = {0.5, 1, 1.5};
const int gampLEN = 2; double gampvals[] = {0.01, 0.02};
const int dLEN = 1; double dvals[] = {0.00548};
const int bLEN = 1; double bvals[] = {4};

const int NParSets = tvLEN*tbLEN*R0pLEN*RhoLEN*gampLEN*dLEN;

const int NumPars = 12; //Number of columns in ParMat
const bool VerboseWriteFlag = false;

//************************
//USER-ASSIGNED VARIABLES
//************************
char SimName[50] = "C_Freq";

std::vector<double> tvVals;
double TVaccMIN = 5*365; double TVaccMAX = 6*365; 
std::vector<double> TVaccStartVals;
std::vector<double> BpVals;
std::vector<double> R0pVals;  
std::vector<double> NvVals;
std::vector<double> RhoVals;
std::vector<int> IpInitVals; 
std::vector<double> tbVals;
std::vector<double> gampVals;
std::vector<double> dVals;

double TMax = 11.0*365.0; double tick = 1.0; //OneSim writes data at time-intervals tick

int SInit = 1000;

//********
//ARRAYS
//********
double ParMat [NParSets][NumPars];
double TExtMat [NParSets][NTrials];

//*******************
//CRITICAL VARIABLES
//*******************
int S, Iv, Ip, V, P, NPop, Par, IpInit, NFails;
char FileNamePar[50];
char FileNameTExt[50];
char DirName[50] = "Data/";
char FileSuffix[50], FileNameDat[50];
double b0, tb, T, Nv,tv, b, gamv, R0p, Bp, gamp, d, Event_Rate, Event_Rate_Prod, RandDeath, dTime, t, ti, TVaccStart;
std::ofstream out_data;

//OTHER VARIABLES
int ntrial, whichindex, nbirths, ndeaths, ninfv, ninfp, nrecv, nrecp, S_death, Iv_death, Ip_death, V_death, P_death, svacc, npopvacc, totvacc, totbirthson, totbirthsoff;
double whichmin, unif, tmodT;
double Nudge = 0.0000001;

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
  srand( time(NULL) );

  //Build directory for simulation results
  strcat(DirName, SimName);  
  mkdir(DirName, ACCESSPERMS);
  strcat(DirName, "/");	  

  //Create filenames within directory for parmat, textmat, ipmat
  strcpy(FileNamePar, DirName);  
  strcat(FileNamePar, "ParMat");

  strcpy(FileNameTExt, DirName);
  strcat(FileNameTExt, "TExtMat");

  Initialize(); //Fill in the parameter matrix

  WriteMat((double *)ParMat, NParSets, NumPars, FileNamePar); //Write ParMat

  for(int Par = 0; Par < NParSets; Par++) //Loop through parameters
    {
      if(VerboseWriteFlag){
	strcpy(FileNameDat, DirName);
	sprintf(FileSuffix, "Par_%d",Par);
	strcat(FileNameDat, FileSuffix);
	out_data.open(FileNameDat);
	out_data << "time S Iv Ip V P N births deaths ninfv ninfp nrecv nrecp S_death Iv_death Ip_death V_death P_death svacc npopvacc totvacc totbirthson totbirthsoff\n";
      }
      
      //Extract parameters
      //Parmat[0] corresponds to Par;
      b0 = ParMat[Par][1]; //b0
      d = ParMat[Par][2]; //d
      Bp = ParMat[Par][3]; //Bp
      //----Nv is set within for(ntrial = 0...) loop----
      tv = ParMat[Par][5]; //tv
      gamv = ParMat[Par][6]; //gamv
      gamp = ParMat[Par][7]; //gamp
      tb = ParMat[Par][8]; //tb
      T = ParMat[Par][9]; //T
      IpInit = (int) ParMat[Par][10]; // Initial level of pathogen upon invasion
      TVaccStart = ParMat[Par][11]; // Time at which vaccination begins
 
      for(ntrial = 0; ntrial < NTrials; ntrial++)
	{
	  NFails = 0;
	  do {
	    S = SInit; //S
	    Iv = 0; //Iv
	    Ip = IpInit; //Ip introduced immediately
	    V = 0; //V
	    P = 0; //P
	    NPop = S + Iv + Ip + V + P;
	    
	    //Simulate to quasi steady state (rewrites State)
	    Nv = 0.0; //No vaccination at first
	    OneSim(0.0, TVaccStart, true);
	    NFails++;
	  } while(t < ( TVaccStart) && NFails < 100 ); 
	  //above while: Throw away sims that don't last TVaccstart+tv time, 
	  //quit after ntrials failures

	  std::cout << ntrial << "\n";
	    
	  if(t > TVaccStart){
	    //Simulate invasion until TMax years, or pathogen extinction
	    Nv = ParMat[Par][4]; //Nv 
	    OneSim(TVaccStart, TMax, true);
	  }
	    //Store final value of t in TExtMat
	    TExtMat[Par][ntrial] = t;	      
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
  else if(Event_Rate_Prod <= b + d*NPop + Bp*Ip*S/NPop)    //Event: Pathogen infection of S
    {S--; Ip++; ninfp++;}
  else if(Event_Rate_Prod <= b + d*NPop + (Bp*Ip*S + Bp*Ip*Iv)/NPop ) //Event: Pathogen infection of V
    {Iv--; Ip++; ninfp++;}
  else if(Event_Rate_Prod <= b + d*NPop + (Bp*Ip*S + Bp*Ip*Iv)/NPop + gamv*Iv) //Event: Iv Recovery
    {Iv--;V++;nrecv++;}
  else if(Event_Rate_Prod <= b + d*NPop + (Bp*Ip*S + Bp*Ip*Iv)/NPop + gamv*Iv + gamp*Ip) //Event: Ip Recovery
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

//******************************************
//function to Initialize values of 2D array
//******************************************
void Initialize()
{
  tvVals = Seq(0.1, 364.9, tvLEN);
  //tvVals.assign(tvvals, tvvals + tvLEN);

  tbVals.assign(tbvals, tbvals + tbLEN);

  TVaccStartVals = Seq(TVaccMIN, TVaccMAX, TVaccLEN);

  //BpVals.assign(bpvals, bpvals + BpLEN);
  //BpVals = Seq(0.0105,0.028,BpLEN);

  //R0pVals = Seq(1.1, 5.0, BpLEN);
  R0pVals.assign(r0pvals, r0pvals + R0pLEN);

  IpInitVals.assign(ipinitvals, ipinitvals + IpInitLEN);

  //NvVals = Seq(1.0,500.0,NvLEN);
  //NvVals.assign(nvvals,nvvals+NvLEN);

  RhoVals.assign(rhovals, rhovals + RhoLEN);

  gampVals.assign(gampvals, gampvals + gampLEN);

  dVals.assign(dvals, dvals + dLEN);
  
  //Fill in ParMat
  int i = 0;
  for(int i1=0; i1<tvVals.size(); i1++)
    for(int i2=0; i2<R0pVals.size(); i2++) //Note R0pVals
      for(int i3=0; i3<TVaccStartVals.size(); i3++)
	for(int i4=0; i4<IpInitVals.size(); i4++)
	  for(int i5=0; i5<RhoVals.size(); i5++) //Note RhoVals
	    for(int i6=0; i6<tbVals.size(); i6++)
	      for(int i7=0; i7<gampVals.size(); i7++)
		for(int i8=0; i8<dVals.size(); i8++)
		  {
		    d = dVals[i8];
		    gamp = gampVals[i7];
		    ParMat[i][0] = i; //Par
		    ParMat[i][1] = 4.0;   //b0
		    ParMat[i][2] = d; //d
		    ParMat[i][3] = R0pVals[i2]*(d + gamp); //Bp
		    ParMat[i][4] = RhoVals[i5]*4.0*tbVals[i6]/(365*d); //Nv
		    ParMat[i][5] = tvVals[i1]; //tv
		    ParMat[i][6] = 0.07; //gamv
		    ParMat[i][7] = gampVals[i7]; //gamp
		    ParMat[i][8] = tbVals[i6]; //tb
		    ParMat[i][9] = 365.0; //T
		    ParMat[i][10] = (double) IpInitVals[i4]; //IpInit
		    ParMat[i][11] = TVaccStartVals[i3]; //TVaccStart
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
  V_death = 0; P_death = 0; svacc = 0; npopvacc = 0; 
  totvacc = 0; totbirthson = 0; totbirthsoff = 0;

  while(t < EndTime && (Ip > 0 || !StopOnErad))
    {
      if( t >= ti && VerboseWriteFlag) //Write old values at intervals "tick"
	{
	  out_data << t << " " << S << " " << Iv << " " << Ip << " " << V << " " << P  << " " << NPop << " " << nbirths << " " << ndeaths <<  " " << ninfv << " " << ninfp << " " << nrecv << " " << nrecp << " " << S_death << " " << Iv_death << " " << Ip_death << " " << V_death << " " << P_death << " " << svacc << " " << npopvacc << " " << totvacc << " " << totbirthson << " " <<  totbirthsoff << "\n"; 
	  ti += tick;
	  nbirths = 0; ndeaths = 0; svacc = 0; npopvacc = 0; totvacc = 0; totbirthson = 0; 
	  totbirthsoff = 0; 
	}  
      Event_Rate = b + d*NPop + (Bp*S*Ip + Bp*Iv*Ip)/NPop + gamv*Iv + gamp*Ip;
      Event_Rate_Prod = Event_Rate*Rand();
      
      GetTime(); //Get time to next event
      
      CheckEventConflict(); //Finds whichmin, the time of the next conflict

      if(dTime < whichmin) //If no conflicts, proceed with Gillepsie event
	{
	  ApplyEvent();	
	}
      else{ //If conflict, stop at conflict and perform necessary action

	switch(whichindex)
	  {
	  case 0 : b = b0; totbirthson++; break; //Start of birthing season
	  case 1 : b = 0.0; totbirthsoff++; break; //End of birthing season
	  case 2 : VaccFun(); totvacc++; //Update S,Iv due to vaccination 
	  }//End else
	dTime = whichmin + Nudge; //Nudge ensures that t moves past conflict
      }    

      t += dTime;
      tmodT = std::fmod(t,T);

    }//End While
  
  if(VerboseWriteFlag)
    {
      out_data << t << " " << S << " " << Iv << " " << Ip << " " << V << " " << P  << " " << NPop << " " << nbirths << " " << ndeaths <<  " " << ninfv << " " << ninfp << " " << nrecv << " " << nrecp << " " << S_death << " " << Iv_death << " " << Ip_death << " " << V_death << " " << P_death << " " << svacc << " " << npopvacc << " " << totvacc << " " << totbirthson << " " << totbirthsoff << "\n"; 
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
	  out_data << *((arr+i*NCol) + j) << " "; 
	}
      out_data << std::endl;
    }
  out_data.close();
}
