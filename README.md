This repository contains two directories, ProphVacc and EradVacc, that
contain code to simulate a vaccination campaign applied to a
fluctuating wildlife population. The two directories simulate two
scenarios: in ProphVacc, the pathogen is introduced into a
pathogen-naive population that is annually vaccinated. The goal is to
prevent the pathogen from becoming established in the population. In
EradVacc, the pathogen is introduced into the population at the start
of simulations, and annual vaccination is initiated 5 years
after. Here, the goal is to drive the pathogen to extinction.

Each directory contains Sim_Gillespie, a c++ script that uses the Gillespie algorithm to simulate each scenario forward in time.  