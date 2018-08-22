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
int S_next, S, Ip_next, Ip, Rp_next, Rp, Npop_next, Npop, nbirths, ndeaths, ninf, nrec, nsdeath, nipdeath, nrpdeath;
char FileName[50];
double b, R0p, Bp, delp, d, Sum_rate, UniSeed, RandDeath, Picker, Sum, TimeOfReaction, t, Tmax, tick, ti;

ofstream out_data;

//FUNCTIONS
double Rand () ;

//MAIN
int main()
{

  srand ( time(NULL) );
  
  cout << "Enter filename...\n";
  cin >> FileName;
  
  out_data.open(FileName);

  out_data << "time S Ip Rp N births deaths inf rec nsdeath nipdeath nrpdeath\n"; //Place these names in the open file

  //Set simulation parameters
  t = 0;
  Tmax = 1000;

  //Set biological parameters
  b = 500;
  R0p = 10;
  delp = 0.1;
  d = 0.5;
  Bp = R0p*d*(d + delp) / b;
  /*These parameters correspond to the events:rates: 
Birth:     b
Death:     d*Npop
Infection: Bp*S*Ip
Recovery:  delp*Ip
  */

  //Set initial conditions  
  S = (int) (b / d * 1 / R0p) ; // Set susceptible population at carrying capacity  
  Ip = (int) (b / (d + delp) * (1 - 1 / R0p )) ;
  Rp = (int) (delp / d * b / (d + delp) * (1 - 1 / R0p )) ;
  Npop = S + Ip + Rp;

  cout << "S0: " << S << "; Ip0: " << Ip << "; Rp0: " << Rp << "; Rp0alt:" << delp / d * b / (d + delp) * (1 - 1 / R0p ) << "\n" << "\n"; 
  
  tick = 1; //Writes data after each time interval tick
  ti = 0; // Works with tick
  nbirths = 0;
  ndeaths = 0;
  ninf = 0;
  nrec = 0;
  nsdeath = 0;
  nipdeath = 0;
  nrpdeath = 0;
  
  while(t < Tmax) {
        
    /*
Update the state variables S, Ip, Rp
as functions of the old state variables. 

Event1: Birth:     b
Event2: Death:     d*Npop
Event3: Infection: Bp*S*Ip
Event4: Recovery:  delp*Ip
    */    

    //Get reaction time for next step
    Sum_rate = b + d*Npop + Bp*S*Ip + delp*Ip;
    TimeOfReaction = -log( Rand() ) / Sum_rate;
    Picker = Rand();
    Sum = 0;
    
    //Event: Birth
    if (Picker < b / Sum_rate)
      {
	S_next = S + 1;
	Ip_next = Ip;
	Rp_next = Rp;
	Npop_next = Npop + 1;
	nbirths += 1;
      }    
    Sum += b;

    //Event: Death
    if (Sum / Sum_rate <= Picker && Picker < (Sum + d*Npop) / Sum_rate)
      {
	//Pick the host that dies
	RandDeath = Rand(); 
	if( RandDeath < (double) S/Npop)
	  {
	    S_next = S - 1; 
	    Ip_next = Ip;
	    Rp_next = Rp;
	    nsdeath ++ ;
	  }
	else {//Host of class Ip or Rp dies
	  if(RandDeath < (double) (S + Ip)/Npop)
	    {
	      S_next = S;
	      Ip_next = Ip - 1;
	      Rp_next = Rp;
	      nipdeath ++ ;
	    }
	  else {//Host of class Rp dies
	      S_next = S;
	      Ip_next = Ip;
	      Rp_next = Rp - 1;
	      nrpdeath ++ ;
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
	Rp_next = Rp;
	Npop_next = Npop;
	ninf += 1;
      }
    Sum += Bp*Ip*S;

      //Event: Recovery
    if ( Sum / Sum_rate <= Picker && Picker <= (Sum + delp*Ip) / Sum_rate )
      {
	S_next = S;
	Ip_next = Ip - 1;
	Rp_next = Rp + 1;
	Npop_next = Npop;
	nrec += 1;
      }

    //Write old values
    if( t >= ti )
      {
	out_data  << t << " " << S << " " << Ip << " " << Rp  << " " << Npop << " " << nbirths << " " << ndeaths << " " << ninf << " " << nrec << " " << nsdeath << " " << nipdeath << " " << nrpdeath << "\n"; //Write current values    
	cout << "Time: " << t << endl; //"; check sum: " << (Sum + delp*Ip) / Sum_rate << "Check R " << Rp << endl; 
	ti = ti + tick;
	nbirths = 0;
	ndeaths = 0;
	ninf = 0;
	nrec = 0;
      }

    
    //Update current state and time
    S = S_next;
    Ip = Ip_next;
    Rp = Rp_next;
    Npop = Npop_next;
    t += TimeOfReaction;

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

