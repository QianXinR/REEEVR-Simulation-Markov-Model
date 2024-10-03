# REEEVR-Simulation-Markov-Model
This is a simulated Markov model that will randomly generate base case inputs and uncertainties around each input based on the number of health states, treatment arms, and specific assumptions.
We generate state costs values in increasing order, i.e., the more severer the health sate, the higher costs correlated. The final health state(death) is assumed to have no associated costs.
We generate state utility values in decreasing order, i.e., the more severer the health sate, the lower utility correlated. The final health state(death) is assumed to have no associated utility.
We assumed: 

- patient cohort all start from health state 1, and the the cohort changes with each model cycle.
- the cycle length is one year, so no further calculation is conducted transferring utility to QALYs.
- no treatment discontinuation is considered, patients taking the treatments within the time horizon, except when they moved into the last health states.

The randomly generated inputs, along with the associated uncertainty for each input, are then used to test the reliability of the value of information (VOI) analysis, particularly the expected value of perfect information (EVPPI), using various methods.
