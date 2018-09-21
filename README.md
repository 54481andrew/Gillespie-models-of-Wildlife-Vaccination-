This repository contains two directories, ProphVacc and EradVacc, that
contain code to simulate a vaccination campaign applied to a
fluctuating wildlife population. The three directories simulate three 
scenarios: 

******************
    SCENARIOS 
******************

*NewPropVacc: the pathogen is introduced into a completely susceptible, fluctuating population. A single premptive or reactive pulse vaccination occurs. Prior to pathogen invasion, the susceptible population is simulated for 8 years to reach a stable limit cycle. Output figures depict the probability that the pathogen remains in the population, one year after vaccination.  

*ProphVacc: the pathogen is introduced into a population that has been annually vaccinated for 8 years prior. The simulation records the length of time that the pathogen remains in the population. Output figures depict the probability that the pathogen remains in the population, one year after vaccination.  

*EradVacc: the pathogen has been circulating in a fluctuating population for 8 years. On year 8, annual pulse vaccinations begin. Output figures show the proportion of simulations in which the pathogen remains after 1 year.   


******************
       CODE
******************

Each directory contains c++ scripts Sim_Gillespie and Sim_Gillespie_Freq that use a modified Gillespie algorithm to simulate the scenario forward in time. The Gillespie algorithm is modified so that simulation stops for the following fixed time events: start of breeding season (birth rate set nonzero), end of breeding season (birth rate set to zero), and a pulse vaccination. 

In addition, each directory contains a Data directory with R graphing scripts. The graphing scripts plot 2 types of figures: the first shows the probability that the pathogen remains in the population for a fixed time; the 2nd type of figure is a line graph of the same probability, averaged over all of the possible times of pathogen invasion throughout a year. 
