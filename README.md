# PIGSS

Plants (Individuals) and Gametes Stochastic Simulation. Code for simulating plant populations using a Cleave and Rescue gamete killer

## How to use

To run the simulations, the files Simulation.py and PlantClasses.py are required, as well as an empty folder called "final_data"

This module was written using Python 3.7, and requires the following python packages:
* pandas
* time
* tqdm
* numpy

The code is run from within the working directory, using the command

    python Simulation.py function_name

Many functions are already written, and will automatically output the data to a folder named "final_data." Here is a list of function names, the data they generate, and the names of the files they will create. All file names here are listed as "onePartner," however for other number of partners, the file name matches such that a function called "fivePartner_..." will have file name "fivePartner..."

Function Names | Data Simulated | Output file names 
--- | --- | --- 
explore_fitness_costs_diploid | 301 | 283 
onePartner_explore_fitness_costs_haploid, fivePartner_explore_fitness_costs_haploid, twentyPartner_explore_fitness_costs_haploid| explores haploid fitness costs [blah blah blah] for given number of partners | onePartner_haploidFitCosts_noSterile_adults.csv, onePartner_haploidFitCosts_noSterile_allele.csv, onePartner_haploidFitCosts_noSterile_total.csv
intro20_onePartner_explore_fitness_costs_haploid, intro20_fivePartner_explore_fitness_costs_haploid, intro20_twentyPartner_explore_fitness_costs_haploid| the same as twentyPartner_explore_fitness_costs_haploid, but with an introduction frequency of 20% instead of 10% | int20_onePartner_haploidFitCosts_noSterile_adults.csv, int20_onePartner_haploidFitCosts_noSterile_allele.csv, int20_onePartner_haploidFitCosts_noSterile_total.csv. | 269
onePartner_female_sterility_MC, fivePartner_female_sterility_MC, twentyPartner_female_sterility_MC | blah | blah 
intro20_onePartner_female_sterility_MC, intro20_fivePartner_female_sterility_MC, intro20_onePartner_female_sterility_MC | blah | blah 
onePartner_male_sterility_MC, fivePartner_male_sterility_MC, twentyPartner_male_sterility_MC | blah | blah 
intro20_onePartner_male_sterility_MC, intro20_fivePartner_male_sterility_MC, intro20_onePartner_male_sterility_MC | blah | blah 

## How it works

This model uses three classes, stored in PlantClasses.py, and two main functions, stored in Simulation.py, in order to simulate the classes. The lowest level class is Diploid, which tracks the genotypes and alleles of a single individual. The next highest class is Haploid, which tracks its own alleles and also its Diploid parent. The highest class, in which most of the data is stored, is StochasticSim. This class object contains all the parameters necessary to run a simulation, and can store the individuals that are created over the course of a simulation. The parameters stored by StochasticSim include the alleles and genotypes possible, the haploid and diploid fitnesses associated with each genotype, the fecundity of each individual, and other parameters used over the course of the simulation. The function stochastic_sim calls on a StochasticSim object to perform the simulation, and for each generation it takes the previous generation, generates gametes, mates gametes, and then produces adult diploid individuals for the next generation. The function stochastic_sim is called by run_stochastic_sim, which handles doing multiple runs and writing the data out to files. More on the specifics of stochastic_sim below.

![flowchart showing the process or run_stochastic_sim. The steps are in a circle, with adults of a generation being represented by purple rectangles, the eggs and sperm represented by blue, and diploid offspring which have not yet reached maturity are colored red](./flowchart_3.png)

Each generation starts with a pool of adult individuals, shown at the top in purple. This information is stored as a list where each index represents a possible genotype, and the associated value represents the number of individuals of that genotype. 
1. For each mating, a single mother produces a pool of ovules and multiple males produce a single pool of pollen.
2. The number of ovules and pollen are reduced, based on their haploid fitness costs Fh such that a pollen with no fitness costs will have a 100% chance of survival, and a pollen with fitness cost 0.1 will only have a 90% chance of surviving.
3. From the pool of possible pollen and ovules, each ovule is mated with a single pollen to produce a possible offspring.
4. The possible offspring are randomly grouped as male or female
5. The number of possible offspring is reduced to only those that survive. For each mating, the base chance of survival is 2 over the number of expected offspring o, as we expect each mating between a female and a male to produce 2 offspring when the population is at carrying capacity. This probability of survival is further modified by Fd, the fitness of the individual, and S. S denotes the density dependence function $S(P) = \frac{g}{1+(g-1)\frac{P}{K}}$ where P is the parent generation’s population size, and K is the carrying capacity. S ranges from some growth factor g > 1 at low densities, to 1 at densities near carrying capacity. The offspring that survive become the parents of the next generation.