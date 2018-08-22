/*To do: 
Gillespie.cpp: 
-add pulse vaccination at times tv%Tb==0
*/



/*
Build a simple c++ script that will 
simulate a stochastic SIR epidemiological
model. The simulation procedes according
to the Gillespie algorithm.  
*/

#include <iostream> // input/output cout, cin
#include <math.h> // Math functions
#include <stdlib.h>
#include <time.h>
#include <fstream> //Used to open files to which data are saved
#include <string>

using namespace std;

//VARIABLES
int S_next, S, Iv_next, Iv, Ip_next, Ip, V_next, V, P_next, P, Npop_next, Npop, nbirths, ndeaths, ninfv, ninfp, nrecv, nrecp, S_death, Iv_death, Ip_death, V_death, P_death, NVacc;
char FileName[50];
double b0, tb, Tb, tv, b, gamv, R0p, Bp, gamp, d, Sum_rate, UniSeed, RandDeath, Picker, Sum, dTime, t, Tmax, tick, ti;
const double pi = 3.1415926535;

ofstream out_data;

//FUNCTION DECLARATIONS
double Rand () ;
double BirthRate ( double ) ;
double VaccFun ( double ) ;

//MAIN
int main()
{

  srand ( time(NULL) );
  

  cout << "Enter filename...\n";
  cin >> FileName;
  //FileName = "1";
  
  out_data.open(FileName);

  out_data << "time S Iv Ip V P N births deaths ninfv ninfp nrecv nrecp S_death Iv_death Ip_death V_death P_death\n"; //Place these names in the open file

  //************//
  // PARAMETERS //
  //************//
  
  //Set simulation parameters
  t = 0;
  Tmax = (double) 2*365;

  //Set biological parameters
  b0 = 500;
  tb = 30;
  Tb = 365;
  d = (double) 0.001;

  //VACCINATION
  NVacc = 1000;
  
  //PATHOGEN
  Bp = 0.00001;
  gamp = 0.001;
  R0p = 0;
  
  /*These parameters correspond to the events and rates at which the events occur: 
Birth baseline:     b0 + Ab sin(2*Pi*t/Tb)
Death:     d*Npop
Infection: Bp*S*Ip
Recovery:  gamp*Ip
  */

  //Set initial conditions  
  S = (int) 1000 ; // Set susceptible population at carrying capacity  
  Iv = (int) 0;
  Ip = (int) 100;
  P = (int) 0;
  V = (int) 0;
  Npop = S + Iv + Ip + P + V;

  //  cout << "d: " << d << "\n";
  cout << "S0: " << S << "; Iv0: " << Iv << "; Ip0: " << Ip << "; V0: " << V << "; P0: " << P << "\n" << "\n"; 

  
  tick = 1; //This program writes data after each time interval tick
  ti = 0; // Works with tick

  nbirths = 0;
  ndeaths = 0;

  ninfv = 0;
  ninfp = 0;
  
  nrecv = 0;
  nrecv = 0;
  
  S_death = 0;
  Iv_death = 0;
  Ip_death = 0;
  V_death = 0;
  P_death = 0;
  
  while(t < Tmax) {
        
    /*
Update the state variables S, Ip, P
as functions of the old state variables. 

Event1: Birth:     b
Event2: Death:     d*Npop
Event3: P Infection: Bp*S*Ip
Event4: V Recovery: gamv*Iv
Event5: P Recovery:  gamp*Ip
Event6: Vaccination: when t <= tv and t + dTime > tv

    */    

    //Get reaction time for next step
    b = BirthRate( t ) ;    
    Sum_rate = b + d*Npop + Bp*S*Ip + gamv*Iv + gamp*Ip
    Picker = Rand();
    Sum = 0;

    
    
    //Event: Birth
    if (Picker < b / Sum_rate)
      {
	S_next = S + 1;
	Iv_next = Iv;
	Ip_next = Ip;
	V_next = V;
	P_next = P;
	Npop_next = Npop + 1;
	nbirths += 1;
      }    
    Sum += b;

    //Event: Death
    if (Sum / Sum_rate <= Picker && Picker < (Sum + d*Npop) / Sum_rate)
      {
	//Pick the host that dies
	RandDeath = Rand();
	//Host of class S dies
	if( RandDeath < (double) S/Npop)
	  {
	    S_next = S - 1;
	    Iv_next = Iv;
	    Ip_next = Ip;
	    V_next = V;
	    P_next = P;
	    S_death ++ ;
	  }
	else {//Host of class Iv dies
	  if(RandDeath < (double) (S + Iv)/Npop)
	    {
	      S_next = S;
	      Iv_next = Iv - 1;
	      Ip_next = Ip;
	      V_next = V;
	      P_next = P;
	      Iv_death ++ ;
	    }
	  else {//Host of class Ip dies
	    if(RandDeath < (double) (S + Iv + Ip)/Npop)
	      {
		S_next = S;
		Iv_next = Iv;
		Ip_next = Ip - 1;
		V_next = V;
		P_next = P;		
		Ip_death ++ ;
	      }
	    else{//Host of class V dies
	      if(RandDeath < (double) (S + Iv + Ip + V)/Npop)
		{
		  S_next = S;
		  Iv_next = Iv;
		  Ip_next = Ip;
		  V_next = V - 1;
		  P_next = P;
		  V_death ++ ;		  
		}
	      else{//Host of class P dies
		  S_next = S;
		  Iv_next = Iv;
		  Ip_next = Ip;
		  V_next = V;
		  P_next = P - 1;
		  P_death ++ ;		  
	      }
	    }
	  }
	}
	Npop_next = Npop - 1;
	ndeaths += 1;
      }
    Sum = Sum + d*Npop;

    //Event: Pathogen infection
    if (Sum / Sum_rate <= Picker && Picker < (Sum + Bp*Ip*S) / Sum_rate)
      {
	S_next = S - 1;
	Ip_next = Ip + 1;
	P_next = P;
	Npop_next = Npop;
	ninfp += 1;
      }
    Sum += Bp*Ip*S;

      //Event: Vaccine recovery 
    if ( Sum / Sum_rate <= Picker && Picker <= (Sum + gamv*Iv) / Sum_rate )
      {
	S_next = S;
	Iv_next = Iv - 1;
	V_next = V + 1;
	Npop_next = Npop;
	nrecv += 1;
      }
    Sum += gamv*Iv;
    
      //Event: Pathogen recovery
    if ( Sum / Sum_rate <= Picker && Picker <= (Sum + gamp*Ip) / Sum_rate )
      {
	S_next = S;
	Ip_next = Ip - 1;
	P_next = P + 1;
	Npop_next = Npop;
	nrecp += 1;
      }

    //Vaccination: this is assumed to be a sudden pulse
    if ( t < tv && t + dTime > tv )
      {
	S_next = S - VaccFun(S,Npop);
	Iv_next = Iv + VaccFun(S,Npop);
      }
    
    
    //Write old values
    if( t >= ti )
      {
  out_data << t << " " << S << " " << Iv << " " << Ip << " " << V << " " << P  << " " << Npop << " " << nbirths << " " << ndeaths <<  " " << ninfv << " " << ninfp << " " << nrecv << " " << nrecp << " " << S_death << " " << Iv_death << " " << Ip_death << " " << V_death << " " << P_death << "\n"; //Write current values    
	cout << "Time: " << t << endl; //"; check sum: " << (Sum + gamp*Ip) / Sum_rate << "Check R " << P << endl; 
	ti = ti + tick;
	nbirths = 0;
	ndeaths = 0;
	//ninf = 0;
	//nrec = 0;
      }

    
    //Update current state and time
    S = S_next;
    Iv = Iv_next;
    Ip = Ip_next;
    V = V_next;
    P = P_next;
    Npop = Npop_next;
    t += dTime;

  }//End While
}//End Main

//////////////////////////////////////
///////// DEFINE FUNCTIONS ///////////
//////////////////////////////////////

//Function returning a random number between 0 and 1
double Rand (){
  double unif;
  unif = 0;
  while (unif == 0)
    {
      unif = (double) rand() / RAND_MAX ;
   }
  return unif;
}  

//Birth function
 double BirthRate( double time ) {
   //Local function declaration
   double result;
   double modtime;
   modtime = fmod(time, Tb ) ;
   result = b0*(modtime < tb) ;
   return result;
 }

//Vaccination function
int VaccFun( int S, int Npop) {
  int result;
  double temp;
  numvacc = round(NVacc*S/(Npop)) ;
  result = min( numvacc , S ) ;
  return result;
}
