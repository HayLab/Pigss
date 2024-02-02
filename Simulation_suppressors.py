"""
Author(s)     : Michelle Johnson and Tobin Ivy
Lab           : Hay Lab
Description   : Performing stochastic simulations on gene drive, where individual gametes
are tracked. Fitness costs for males and females must be roughly equal, as the simulation
depends on having equal/similar numbers of males and females for correct mating

NOTE: sex determination is not currently handled, assumes 50/50 chance of being male/female
"""
import sys
import random
import time
import copy
import pandas as pd
from tqdm import tqdm
from SuppressorMutationClasses import *

# import cProfile
# import pstats

ALLELES = [['C', 'A'], ['V', 'W']]
POP_MAX = 10000
NUM_GENS = 50
CLEAVAGE_PROB = 0.95
NUM_REPS = 10
NUM_OVULES = 30
GROWTH_FACTOR = 6
R_D = [[50], [50]]

BASIC_ADD_INTRO = [[]]
BASIC_RUN_LABEL = "test"
BASIC_MC_PROB = 0
BASIC_FILE_NAME = "blankTest"

def run_stochastic_sim(alleles, num_reps,
                   num_gens, intro, f_c, hf_c, s_c, num_partners, mut_flag,
                   r_d = R_D, cleavage_prob = CLEAVAGE_PROB, add_intro = BASIC_ADD_INTRO,
                   run_label = "test", k = POP_MAX, n_o = NUM_OVULES, g_f = GROWTH_FACTOR,
                   mc_prob = BASIC_MC_PROB, file_name = BASIC_FILE_NAME, FC_flag = None):
    """ a function that performs a stochastic simulation, given all starting parameters
    params:
        alleles     -- list of locis, with each loci being a list of possible alleles at that loci.
                        it is assumed that wt goes last
        num_reps    -- int number of times the stochastic simulation should repeat
        num_gens    -- int number of generations the simulation should run for
        intro       -- list containing information for the introduction of gene drive,
                        takes the form [sex, genotype, frequency] where genotype is an int,
                        derived from Simulation.genotypes. Assume the 0 position refers to the 
                        completely non-wildtype. In this case, [CC, VV]
        f_c         -- list of fitness costs associated with diploids. Each takes the form
                        [sex, list of alleles required for fitness cost,
                        fitness cost, list of alleles that rescue]
        hf_c        -- list of fitness costs associated with haploids. Each takes the form
                        [sex, list of alleles required for fitness cost, 
                        fitness cost, list of PARENT alleles that rescue],
                        with the parent alleles rescuing via maternal carryover
        s_c         -- list of sterility costs, each sterility cost takes the form
                        [sex, list of alleles required, sterility cost]
        num_partners -- integer number for how many male individuals' pollen are
                        pooled together for each female

    optional params:
        r_d         -- list of recombination distances, 
                        values range from 0 (co-inherited) to 50 (separate chromosomes)
        cleavage_prob -- float probability of cleavage
        add_intro   -- list of additional releases, same form as intro
        run_label   -- label for the runs, can be used to identify runs and under what parameters
        k           -- int carrying capacity
        n_o         -- int number of ovules, is the average number of offspring each female has
        g_f         -- int growth factor. The low density growth rate
        mc_prob     -- float probability that maternal carryover occurs
        file_name   -- start of file name to store data in
    """

    cross_dict = {}
    gametogenesis_dict = {}

    adult_file_name = 'arabidopsis_data/' + file_name + '_adults.csv'
    allele_file_name = 'arabidopsis_data/' + file_name + '_NEWallele.csv'
    total_file_name = 'arabidopsis_data/' + file_name + '_total.csv'

    sim_object = StochasticSim(num_gens, alleles, intro, f_c, hf_c, mc_prob, s_c, cross_dict,
                                gametogenesis_dict, r_d,
                                add_intro, cleavage_prob, k, n_o, g_f, num_partners, mut_flag,
                                FC_flag)

    for rep in range(1,num_reps+1):
        rep_label = run_label + "_" + str(rep)

        start_time = time.time()
        data = stochastic_sim(sim_object, rep_label)
        end_time = time.time()
        print(f'time taken on run {rep}: {end_time - start_time}')

        # write information about genotypes to file
        df_adults = data[0]
        df_adults.to_csv(adult_file_name, mode = 'a')

        # write information about alleles to file
        df_alleles = data[1]
        df_alleles.to_csv(allele_file_name, mode = 'a')

        # write information about population to file
        df_total = data[2]
        df_total.to_csv(total_file_name, mode = 'a')

        # transfer over dictionary
        sim_object = data[3]

    print(f'appended {run_label} to file {file_name}')
    return None

def stochastic_sim(Simulation: StochasticSim, label):
    """ function that performs a stochastic simulation: Just the loop through the generations!
     Everything else should be handled prior to that.
    Params:
        Simulation  -- a simulation object
        label       -- the "run" label 
    Returns:
        df_adults   -- a dataframe containing each genotype possible, and the number of adults
                        of that genotype, for each generation simulated
        df_alleles  -- a dataframe containing each allele, and the count of those alleles
                        in the population for each generation simulated
        df_total    -- a dataframe containing total number of females and males for each generation
        Simulation  -- the input simulation object, with possible edits to dictionaries
    """

    working_adults = copy.deepcopy(Simulation.adults)
    for gen in tqdm(range(Simulation.num_gens)):
        # add an empty list for the next generation
        working_adults[0].append([])
        working_adults[1].append([])

        # add individuals
        # TODO: this is PROBABLY broken, I haven't looked at it
        # for a_i, add_release_list in enumerate(Simulation.additonal_release_list):
        #     if gen in add_release_list:
        #         for individual in range(int(Simulation.add_intro[a_i][2]*Simulation.K)):
        #             working_adults[Simulation.add_intro[a_i][0]][gen].append(
        #                 Simulation.genotypes[Simulation.add_intro[a_i][0]][Simulation.add_intro[a_i][1]])
      
        # check that we are not at zero populations
        if sum(working_adults[1][gen]) != 0:
            # initialize temporary adults, they are sexless
            temp_adults = np.zeros(len(Simulation.genotypes))
            # determine survival modifier
            # survival modifier must range from ~growth factor at population near zero,
            #   to 1 at carrying capacity
            survival_modifier = Simulation.growth_factor/(1+(Simulation.growth_factor-1)*
                                                (sum(working_adults[0][gen])+sum(working_adults[1][gen]))
                                                /Simulation.k)

            # to generate gametes, get mother and father information
            # format is [# of individuals with genotype 0, # of individuals with genotype 1 ...]
            mother_counts = working_adults[0][gen]
            father_counts = working_adults[1][gen]

            if sum(mother_counts) > sum(father_counts):
                print("\nERROR: less fathers than mothers")

            # get ovules and pollens, grouped by parent
            grouped_ovules, grouped_pollens = Simulation.get_gametes(mother_counts, father_counts)
            random.shuffle(grouped_ovules)
            # TODO: grouped ovules take the form 
            random.shuffle(grouped_pollens)

            # store a random integer, for later random choosing of fathers.
            # this is done OUTSIDE of the loop, as rng.choice inside the loop is time-intensive
            start_pos = len(grouped_pollens)
            for ovules in grouped_ovules:
                # remove empty fathers
                while [] in grouped_pollens:
                    grouped_pollens.remove([])

                # when females are the introduced population, may lead to 
                #   more females than males, leading to not enough matings
                if grouped_pollens != []:
                    # if we have RUN OUT of fathers, re-use buckets
                    if len(grouped_pollens) < Simulation.num_partners:
                        grouped_pollens = [[bucket] for father in grouped_pollens for bucket in father]

                    # we generate pseudo-random numbers, as rng.choice is time-intensive
                    # max_pos has +1 so that we can use every number

                    # note: we randomize start position for pulling our fathers, as we do not want
                    #   the same fathers to be our sperm pool each time
                    max_pos = len(grouped_pollens) + 1 - Simulation.num_partners
                    start_pos = int(len(grouped_ovules) * start_pos 
                                    + sum(working_adults[0][gen])*Simulation.k) % max_pos

                    # choose our pseudo-random fathers
                    father_indexes = range(start_pos, start_pos + Simulation.num_partners)
                    pollens = []

                    # loop through 'fathers' to pool all pollens
                    for index in father_indexes:
                        father = grouped_pollens[index]
                        bucket = father.pop(0)
                        pollens.extend(bucket)

                    ###### apply fitness costs ##############
                    # get random integers
                    rand_floats = rng.random(len(pollens))
                    # apply fitness costs to each pollen
                    surviving_pollens = [pollen for index, pollen in enumerate(pollens)
                                         if rand_floats[index] < Simulation.haplo_fitness[(pollen, 1)]]

                    # apply fitness costs to each ovule
                    rand_floats = rng.random(len(ovules))
                    # apply fitness costs to each pollen
                    surviving_ovules = [ovule for index, ovule in enumerate(ovules)
                                        if rand_floats[index] < Simulation.haplo_fitness[(ovule, 0)]]

                    # select mating pollens
                    sperm = []
                    # note: there should ALWAYS be an excess of pollen
                    if len(surviving_pollens) >= len(surviving_ovules):
                        sperm = list(rng.choice(surviving_pollens, len(surviving_ovules),
                                                replace = False))
                    else:
                        # the only case where this should not be true, is if males infertile.
                        sperm = surviving_pollens

                    for index, single_sperm in enumerate(sperm):
                        mother_egg = surviving_ovules[index]
                        father_sperm = single_sperm
                        individual_index = Simulation.cross_dict[(mother_egg, father_sperm)]
                        temp_adults[individual_index] += 1

            # given individuals, get sex
            temp_females = [rng.binomial(adult_count, 0.5) for adult_count in temp_adults]
            # we MUST ENSURE that there are more males than females, by at least 1, for mating
            while sum(temp_females) > int(sum(temp_adults)/2):
                temp_females = [rng.binomial(adult_count, 0.5) for adult_count in temp_adults]
            temp_males_1 = [adults_count - temp_females[a_index]
                            for a_index, adults_count in enumerate(temp_adults)]


            # given individuals, get survival
            # probability of survival is (2/ number of offspring) because each mating should
            #   produce two offspring at carrying capacity
            # survival is then modified by fitness and the survival_modifier
            temp_females = [rng.binomial(female_count, (2/Simulation.num_ovules) 
                                         * Simulation.fitness[0][f_index]*survival_modifier)
                            for f_index, female_count in enumerate(temp_females)]
            # number of offspring is num_ovules for both females and males, because num_ovules
            #   is the number of individuals produced; pollen should be in excess
            temp_males = [rng.binomial(male_count, (2/Simulation.num_ovules) 
                                       * Simulation.fitness[1][m_index]*survival_modifier)
                          for m_index, male_count in enumerate(temp_males_1)]
            # if number of males is < females, repeat
            while sum(temp_males) < sum(temp_females):
                temp_males = [rng.binomial(male_count, (2/Simulation.num_ovules) 
                                           * Simulation.fitness[1][m_index]*survival_modifier)
                          for m_index, male_count in enumerate(temp_males_1)]

            # add counts to adults
            working_adults[0][gen+1] = temp_females
            working_adults[1][gen+1] = temp_males

        # must maintain zero vectors if we are at zero populations
        else:
            temp_adults = np.zeros(len(Simulation.genotypes))
            working_adults[0][gen+1] = temp_adults
            working_adults[1][gen+1] = temp_adults

    # create a temporary variable to hold array versions of adults
    adults_temp = [[], []]
    # turn lists into arrays
    for sex in [0, 1]:
        adults_temp[sex] = np.array([np.array(x) for x in working_adults[sex]])

    # generate dataframe of Population Counts (total # females & males in each generation)
    # summarize female adults data in a dataframe
    df_total_females = pd.DataFrame(np.sum(adults_temp[0], axis = 1), columns = ['Count'])
    df_total_females['Generation'] = range(Simulation.num_gens+1)
    df_total_females['Sex'] = 'Female'
    # summarize male adults data in a dataframe
    df_total_males = pd.DataFrame(np.sum(adults_temp[1], axis = 1), columns = ['Count'])
    df_total_males['Generation'] = range(Simulation.num_gens+1)
    df_total_males['Sex'] = 'Male'
    # combine dataframes
    df_total = pd.concat([df_total_females, df_total_males])
    df_total = df_total.reset_index(drop=True)
    df_total['Run'] = label

    # generate df of genotype counts
    # summarize female by genotype data into dataframe
    df_females = pd.DataFrame(adults_temp[0], columns = [str(diploid.genotype) 
                                                         for diploid in Simulation.genotypes])
    df_females['Generation'] = range(Simulation.num_gens+1)
    df_females = df_females.melt(id_vars = 'Generation', var_name = 'Genotype',
                                 value_name = 'Count')
    df_females['Sex'] = 'Female'
    # repeat for males
    df_males = pd.DataFrame(adults_temp[1], columns = [str(diploid.genotype) 
                                                       for diploid in Simulation.genotypes])
    df_males['Generation'] = range(Simulation.num_gens+1)
    df_males = df_males.melt(id_vars = 'Generation', var_name = 'Genotype', value_name = 'Count')
    df_males['Sex'] = 'Male'
    # combine dataframes
    df_adults = pd.concat([df_females, df_males])
    df_adults = df_adults.reset_index(drop=True)
    df_adults['Run'] = label

    # generate df of allele counts
    alleles_temp = [[], []]
    for sex in [0, 1]:
        for allele in Simulation.all_alleles:
            a_count = [0]*(Simulation.num_gens+1)
            for d_index, diploid in enumerate(Simulation.genotypes):
                allele_list = [allele for locus in diploid.alleles for allele in locus]
                # no idea if this works
                a_count += adults_temp[sex][:,d_index] * allele_list.count(allele)
            alleles_temp[sex].append(a_count)
        alleles_temp[sex] = np.array([np.array(x) for x in alleles_temp[sex]])

    df_females_a = pd.DataFrame(np.transpose(alleles_temp[0]), 
                                columns = list(Simulation.all_alleles))
    df_females_a['Generation'] = range(Simulation.num_gens+1)
    df_females_a = df_females_a.melt(id_vars = 'Generation', var_name = 'Allele',
                                     value_name = 'Count')
    df_females_a['Sex'] = 'Female'

    df_males_a = pd.DataFrame(np.transpose(alleles_temp[1]), columns = list(Simulation.all_alleles))
    df_males_a['Generation'] = range(Simulation.num_gens+1)
    df_males_a = df_males_a.melt(id_vars = 'Generation', var_name = 'Allele', value_name = 'Count')
    df_males_a['Sex'] = 'Male'
    df_alleles = pd.concat([df_females_a, df_males_a])
    df_alleles = df_alleles.reset_index(drop=True)
    df_alleles.drop(df_alleles.index[(df_alleles["Sex"] == 'Female') 
                                     & (df_alleles["Allele"] == 'Y')].tolist(), 
                    inplace = True)
    df_alleles['Run'] = label

    return(df_adults, df_alleles, df_total, Simulation)


##########################################
######### population modification ########
##########################################

def resistant_modification_onePartner_haploid_fc():
    """runs an arabidopsis cleaver stochastic simulation, with varying fitness costs.
    Data is all stored in three files:
    final_data/resistant_modification_onePartner_haploidFitCosts_noSterile_total.csv
    final_data/resistant_modification_onePartner_haploidFitCosts_noSterile_NEWallele.csv
    final_data/resistant_modification_onePartner_haploidFitCosts_noSterile_adults.csv
    """
    num_partners = 1
    alleles = [['C', 'R', 'A'], ['V', 'W']] # a resistance allele is uncleavable
    intro = [[1, 0, 0.1], [0, 23, 0.025], [1, 23, 0.025]] # sex, genotype, frequency
    # genotype 0 = cc vv, genotype 23 = ra, ww (wt resistant)
    s_c = []
    f_c = []

    file_name = "mutation_data/resistant_modification_onePartner_haploidFitCosts_noSterile"

    for clvr_cost in [0, 0.05, 0.1, 0.15, 0.2, 0.4]:
        # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
        hf_c = [[0, ['C', 'W'], 1.0, [['C']]], #savable via maternal carryover
                # NOTE. 'C' in maternal carryover does NOT matter here, because MC = 0
                [1, ['C', 'W'], 1.0, []],
                [0, ['V'], clvr_cost, []], # haploid fitness cost
                [1, ['V'], clvr_cost, []]] # haploid fitness cost

        run_label = f'clvr_cost_{clvr_cost}_intro_{intro[0][2]}'
        run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro, f_c, hf_c, s_c, num_partners,
                           mut_flag = "NA",
                           run_label= run_label,
                           file_name= file_name, k=POP_MAX)
    return None

def resistant_modification_fivePartner_haploid_fc():
    """runs an arabidopsis cleaver stochastic simulation, with varying fitness costs.
    Data is all stored in three files:
    final_data/resistant_modification_fivePartner_haploidFitCosts_noSterile_total.csv
    final_data/resistant_modification_fivePartner_haploidFitCosts_noSterile_NEWallele.csv
    final_data/resistant_modification_fivePartner_haploidFitCosts_noSterile_adults.csv
    """
    num_partners = 5
    alleles = [['C', 'R', 'A'], ['V', 'W']] # a resistance allele is uncleavable
    intro = [[1, 0, 0.1], [0, 23, 0.025], [1, 23, 0.025]] # sex, genotype, frequency
    # genotype 0 = cc vv, genotype 23 = ra, ww (wt resistant)
    s_c = []
    f_c = []

    file_name = "mutation_data/resistant_modification_fivePartner_haploidFitCosts_noSterile"

    for clvr_cost in [0, 0.05, 0.1, 0.15, 0.2, 0.4]:
        # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
        hf_c = [[0, ['C', 'W'], 1.0, [['C']]], #savable via maternal carryover
                # NOTE. 'C' in maternal carryover does NOT matter here, because MC = 0
                [1, ['C', 'W'], 1.0, []],
                [0, ['V'], clvr_cost, []], # haploid fitness cost
                [1, ['V'], clvr_cost, []]] # haploid fitness cost

        run_label = f'clvr_cost_{clvr_cost}_intro_{intro[0][2]}'
        run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro, f_c, hf_c, s_c, num_partners,
                           mut_flag = "NA",
                           run_label= run_label,
                           file_name= file_name, k=POP_MAX)
    return None

def resistant_modification_twentyPartner_haploid_fc():
    """runs an arabidopsis cleaver stochastic simulation, with varying fitness costs.
    Data is all stored in three files:
    final_data/resistant_modification_fivePartner_haploidFitCosts_noSterile_total.csv
    final_data/resistant_modification_fivePartner_haploidFitCosts_noSterile_NEWallele.csv
    final_data/resistant_modification_fivePartner_haploidFitCosts_noSterile_adults.csv
    """
    num_partners = 20
    alleles = [['C', 'R', 'A'], ['V', 'W']] # a resistance allele is uncleavable
    intro = [[1, 0, 0.1], [0, 23, 0.025], [1, 23, 0.025]] # sex, genotype, frequency
    # genotype 0 = cc vv, genotype 23 = ra, ww (wt resistant)
    s_c = []
    f_c = []

    file_name = "mutation_data/resistant_modification_twentyPartner_haploidFitCosts_noSterile"

    for clvr_cost in [0, 0.05, 0.1, 0.15, 0.2, 0.4]:
        # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
        hf_c = [[0, ['C', 'W'], 1.0, [['C']]], #savable via maternal carryover
                # NOTE. 'C' in maternal carryover does NOT matter here, because MC = 0
                [1, ['C', 'W'], 1.0, []],
                [0, ['V'], clvr_cost, []], # haploid fitness cost
                [1, ['V'], clvr_cost, []]] # haploid fitness cost

        run_label = f'clvr_cost_{clvr_cost}_intro_{intro[0][2]}'
        run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro, f_c, hf_c, s_c, num_partners,
                           mut_flag = "NA",
                           run_label= run_label,
                           file_name= file_name, k=POP_MAX)
    return None

# region RMmod, haploid v. diploid v. haplo-diploid

def RM_manypercents_onePartner_haploid_fc():
    """runs an arabidopsis cleaver stochastic simulation, with varying fitness costs.
    Data is all stored in three files:
    final_data/resistant_modification_onePartner_haploidFitCosts_noSterile_total.csv
    final_data/resistant_modification_onePartner_haploidFitCosts_noSterile_NEWallele.csv
    final_data/resistant_modification_onePartner_haploidFitCosts_noSterile_adults.csv
    """
    num_partners = 1
    alleles = [['C', 'R', 'A'], ['V', 'W']] # a resistance allele is uncleavable
    
    # genotype 0 = cc vv, genotype 23 = ra, ww (wt resistant)
    s_c = []
    f_c = []

    file_name = "mutation_data/RM_manypercents_onePartner_haploidFitCosts_noSterile"

    for clvr_cost in [0, 0.05, 0.1, 0.15]:
        for introFreq in [0.1]: #, 0.2]:
            intro = [[1, 0, 0.1], [0, 23, introFreq/2], [1, 23, introFreq/2]] # sex, genotype, frequency
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            hf_c = [[0, ['C', 'W'], 1.0, [['C']]], #savable via maternal carryover
                    # NOTE. 'C' in maternal carryover does NOT matter here, because MC = 0
                    [1, ['C', 'W'], 1.0, []],
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]] # haploid fitness cost

        run_label = f'clvr_cost_{clvr_cost}_intro_{introFreq}'
        run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro, f_c, hf_c, s_c, num_partners,
                           mut_flag = "NA",
                           run_label= run_label,
                           file_name= file_name, k=POP_MAX)
    return None

def RM_manypercents_onePartner_haploid_fc2():
    """runs an arabidopsis cleaver stochastic simulation, with varying fitness costs.
    Data is all stored in three files:
    final_data/resistant_modification_onePartner_haploidFitCosts_noSterile_total.csv
    final_data/resistant_modification_onePartner_haploidFitCosts_noSterile_NEWallele.csv
    final_data/resistant_modification_onePartner_haploidFitCosts_noSterile_adults.csv
    """
    num_partners = 1
    alleles = [['C', 'R', 'A'], ['V', 'W']] # a resistance allele is uncleavable
    
    # genotype 0 = cc vv, genotype 23 = ra, ww (wt resistant)
    s_c = []
    f_c = []

    file_name = "mutation_data/RM_manypercents_onePartner_haploidFitCosts_noSterile2"

    for clvr_cost in [0, 0.05, 0.1, 0.15]:
        for introFreq in [0.2]:
            intro = [[1, 0, 0.1], [0, 19, introFreq/2], [1, 19, introFreq/2]] # sex, genotype, frequency
            # 23 = ra ww, 19 = rr ww
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            hf_c = [[0, ['C', 'W'], 1.0, [['C']]], #savable via maternal carryover
                    # NOTE. 'C' in maternal carryover does NOT matter here, because MC = 0
                    [1, ['C', 'W'], 1.0, []],
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]] # haploid fitness cost

        run_label = f'clvr_cost_{clvr_cost}_intro_{introFreq}'
        run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro, f_c, hf_c, s_c, num_partners,
                           mut_flag = "NA",
                           run_label= run_label,
                           file_name= file_name, k=POP_MAX)
    return None

def RM_manypercents_onePartner_diploid_fc():
    """runs an arabidopsis cleaver stochastic simulation, with varying fitness costs.
    """
    num_partners = 1
    alleles = [['C', 'R', 'A'], ['V', 'W']] # a resistance allele is uncleavable
    
    # genotype 0 = cc vv, genotype 23 = ra, ww (wt resistant)
    s_c = []
    hf_c = [[0, ['C', 'W'], 1.0, [['C']]], #savable via maternal carryover
                    # NOTE. 'C' in maternal carryover does NOT matter here, because MC = 0
                    [1, ['C', 'W'], 1.0, []]]

    file_name = "mutation_data/RM_manypercents_onePartner_diploidFitCosts_noSterile"

    for clvr_cost in [0, 0.05, 0.1, 0.15]:
        for introFreq in [0.1, 0.2]:
            intro = [[1, 0, 0.1], [0, 23, introFreq/2], [1, 23, introFreq/2]] # sex, genotype, frequency
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            f_c = [[0, ['V'], clvr_cost, []], # haploid fitness cost
                   [1, ['V'], clvr_cost, []]] # haploid fitness cost

        run_label = f'clvr_cost_{clvr_cost}_intro_{introFreq}'
        run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro, f_c, hf_c, s_c, num_partners,
                           mut_flag = "NA",
                           run_label= run_label,
                           file_name= file_name, k=POP_MAX)
    return None

def RM_manypercents_onePartner_diploidHaploid_fc():
    """runs an arabidopsis cleaver stochastic simulation, with varying fitness costs.
    """
    num_partners = 1
    alleles = [['C', 'R', 'A'], ['V', 'W']] # a resistance allele is uncleavable
    
    # genotype 0 = cc vv, genotype 23 = ra, ww (wt resistant)
    s_c = []

    file_name = "mutation_data/RM_manypercents_onePartner_diploidHaploidFitCosts_noSterile"

    for clvr_cost in [0, 0.05, 0.1, 0.15]:
        for introFreq in [0.1, 0.2]:
            intro = [[1, 0, 0.1], [0, 23, introFreq/2], [1, 23, introFreq/2]] # sex, genotype, frequency
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            f_c = [[0, ['V'], clvr_cost, []], # haploid fitness cost
                   [1, ['V'], clvr_cost, []]] # haploid fitness cost

            hf_c = [[0, ['C', 'W'], 1.0, [['C']]], #savable via maternal carryover
                    # NOTE. 'C' in maternal carryover does NOT matter here, because MC = 0
                    [1, ['C', 'W'], 1.0, []],
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]]

        run_label = f'clvr_cost_{clvr_cost}_intro_{introFreq}'
        run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro, f_c, hf_c, s_c, num_partners,
                           mut_flag = "NA",
                           run_label= run_label,
                           file_name= file_name, k=POP_MAX)
    return None

def RM_manypercents_onePartner_haploid_fc_homo():
    """runs an arabidopsis cleaver stochastic simulation, with varying fitness costs.
    Data is all stored in three files:
    final_data/resistant_modification_onePartner_haploidFitCosts_noSterile_total.csv
    final_data/resistant_modification_onePartner_haploidFitCosts_noSterile_NEWallele.csv
    final_data/resistant_modification_onePartner_haploidFitCosts_noSterile_adults.csv
    """
    num_partners = 1
    alleles = [['C', 'R', 'A'], ['V', 'W']] # a resistance allele is uncleavable
    
    # genotype 0 = cc vv, genotype 23 = ra, ww (wt resistant)
    s_c = []
    f_c = []

    file_name = "mutation_data/RM_manypercents_onePartner_haploidFitCosts_noSterile_homoRelease"

    for clvr_cost in [0, 0.05, 0.1, 0.15]:
        for introFreq in [0.05, 0.1, 0.2]:
            intro = [[1, 0, 0.1], [0, 19, introFreq/2], [1, 19, introFreq/2]] # sex, genotype, frequency
            # 23 = ra ww, 19 = rr ww
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            hf_c = [[0, ['C', 'W'], 1.0, [['C']]], #savable via maternal carryover
                    # NOTE. 'C' in maternal carryover does NOT matter here, because MC = 0
                    [1, ['C', 'W'], 1.0, []],
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]] # haploid fitness cost

            run_label = f'clvr_cost_{clvr_cost}_intro_{introFreq}'
            run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro, f_c, hf_c, s_c, num_partners,
                            mut_flag = "NA",
                            run_label= run_label,
                            file_name= file_name, k=POP_MAX)
    return None

def RM_manypercents_onePartner_diploid_fc_homo():
    """runs an arabidopsis cleaver stochastic simulation, with varying fitness costs.
    """
    num_partners = 1
    alleles = [['C', 'R', 'A'], ['V', 'W']] # a resistance allele is uncleavable
    
    # genotype 0 = cc vv, genotype 23 = ra, ww (wt resistant)
    s_c = []
    hf_c = [[0, ['C', 'W'], 1.0, [['C']]], #savable via maternal carryover
                    # NOTE. 'C' in maternal carryover does NOT matter here, because MC = 0
                    [1, ['C', 'W'], 1.0, []]]

    file_name = "mutation_data/RM_manypercents_onePartner_diploidFitCosts_noSterile_homoRelease"

    for clvr_cost in [0, 0.05, 0.1, 0.15]:
        for introFreq in [0.05, 0.1, 0.2]:
            intro = [[1, 0, 0.1], [0, 19, introFreq/2], [1, 19, introFreq/2]] # sex, genotype, frequency
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            f_c = [[0, ['V'], clvr_cost, []], # haploid fitness cost
                   [1, ['V'], clvr_cost, []]] # haploid fitness cost

            run_label = f'clvr_cost_{clvr_cost}_intro_{introFreq}'
            run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro, f_c, hf_c, s_c, num_partners,
                            mut_flag = "NA",
                            run_label= run_label,
                            file_name= file_name, k=POP_MAX)
    return None

def RM_manypercents_onePartner_diploidHaploid_fc_homo():
    """runs an arabidopsis cleaver stochastic simulation, with varying fitness costs.
    """
    num_partners = 1
    alleles = [['C', 'R', 'A'], ['V', 'W']] # a resistance allele is uncleavable
    
    # genotype 0 = cc vv, genotype 23 = ra, ww (wt resistant)
    s_c = []

    file_name = "mutation_data/RM_manypercents_onePartner_diploidHaploidFitCosts_noSterile_homoRelease"

    for clvr_cost in [0, 0.05, 0.1, 0.15]:
        for introFreq in [0.05, 0.1, 0.2]:
            intro = [[1, 0, 0.1], [0, 19, introFreq/2], [1, 19, introFreq/2]] # sex, genotype, frequency
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            f_c = [[0, ['V'], clvr_cost, []], # haploid fitness cost
                   [1, ['V'], clvr_cost, []]] # haploid fitness cost

            hf_c = [[0, ['C', 'W'], 1.0, [['C']]], #savable via maternal carryover
                    # NOTE. 'C' in maternal carryover does NOT matter here, because MC = 0
                    [1, ['C', 'W'], 1.0, []],
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]]

            run_label = f'clvr_cost_{clvr_cost}_intro_{introFreq}'
            run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro, f_c, hf_c, s_c, num_partners,
                            mut_flag = "NA",
                            run_label= run_label,
                            file_name= file_name, k=POP_MAX)
    return None


def RM_manypercents_fivePartner_haploid_fc():
    """runs an arabidopsis cleaver stochastic simulation, with varying fitness costs.
    Data is all stored in three files:
    final_data/resistant_modification_onePartner_haploidFitCosts_noSterile_total.csv
    final_data/resistant_modification_onePartner_haploidFitCosts_noSterile_NEWallele.csv
    final_data/resistant_modification_onePartner_haploidFitCosts_noSterile_adults.csv
    """
    num_partners = 5
    alleles = [['C', 'R', 'A'], ['V', 'W']] # a resistance allele is uncleavable
    
    # genotype 0 = cc vv, genotype 23 = ra, ww (wt resistant)
    s_c = []
    f_c = []

    file_name = "mutation_data/RM_manypercents_fivePartner_haploidFitCosts_noSterile2"

    for clvr_cost in [0, 0.05, 0.1, 0.15]:
        for introFreq in [0.1, 0.2]:
            intro = [[1, 0, 0.1], [0, 19, introFreq/2], [1, 19, introFreq/2]] # sex, genotype, frequency
            # 23 = ra ww, 19 = rr ww
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            hf_c = [[0, ['C', 'W'], 1.0, [['C']]], #savable via maternal carryover
                    # NOTE. 'C' in maternal carryover does NOT matter here, because MC = 0
                    [1, ['C', 'W'], 1.0, []],
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]] # haploid fitness cost

            run_label = f'clvr_cost_{clvr_cost}_intro_{introFreq}'
            run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro, f_c, hf_c, s_c, num_partners,
                            mut_flag = "NA",
                            run_label= run_label,
                            file_name= file_name, k=POP_MAX)
    return None

def RM_manypercents_fivePartner_diploid_fc():
    """runs an arabidopsis cleaver stochastic simulation, with varying fitness costs.
    """
    num_partners = 5
    alleles = [['C', 'R', 'A'], ['V', 'W']] # a resistance allele is uncleavable
    
    # genotype 0 = cc vv, genotype 23 = ra, ww (wt resistant)
    s_c = []
    hf_c = [[0, ['C', 'W'], 1.0, [['C']]], #savable via maternal carryover
                    # NOTE. 'C' in maternal carryover does NOT matter here, because MC = 0
                    [1, ['C', 'W'], 1.0, []]]

    file_name = "mutation_data/RM_manypercents_fivePartner_diploidFitCosts_noSterile"

    for clvr_cost in [0, 0.05, 0.1, 0.15]:
        for introFreq in [0.1, 0.2]:
            intro = [[1, 0, 0.1], [0, 23, introFreq/2], [1, 23, introFreq/2]] # sex, genotype, frequency
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            f_c = [[0, ['V'], clvr_cost, []], # haploid fitness cost
                   [1, ['V'], clvr_cost, []]] # haploid fitness cost

            run_label = f'clvr_cost_{clvr_cost}_intro_{introFreq}'
            run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro, f_c, hf_c, s_c, num_partners,
                            mut_flag = "NA",
                            run_label= run_label,
                            file_name= file_name, k=POP_MAX)
    return None

def RM_manypercents_fivePartner_diploidHaploid_fc():
    """runs an arabidopsis cleaver stochastic simulation, with varying fitness costs.
    """
    num_partners = 5
    alleles = [['C', 'R', 'A'], ['V', 'W']] # a resistance allele is uncleavable
    
    # genotype 0 = cc vv, genotype 23 = ra, ww (wt resistant)
    s_c = []

    file_name = "mutation_data/RM_manypercents_fivePartner_diploidHaploidFitCosts_noSterile"

    for clvr_cost in [0, 0.05, 0.1, 0.15]:
        for introFreq in [0.1, 0.2]:
            intro = [[1, 0, 0.1], [0, 23, introFreq/2], [1, 23, introFreq/2]] # sex, genotype, frequency
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            f_c = [[0, ['V'], clvr_cost, []], # haploid fitness cost
                   [1, ['V'], clvr_cost, []]] # haploid fitness cost

            hf_c = [[0, ['C', 'W'], 1.0, [['C']]], #savable via maternal carryover
                    # NOTE. 'C' in maternal carryover does NOT matter here, because MC = 0
                    [1, ['C', 'W'], 1.0, []],
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]]

            run_label = f'clvr_cost_{clvr_cost}_intro_{introFreq}'
            run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro, f_c, hf_c, s_c, num_partners,
                            mut_flag = "NA",
                            run_label= run_label,
                            file_name= file_name, k=POP_MAX)
    return None

def RM_manypercents_fivePartner_haploid_fc_homo():
    """runs an arabidopsis cleaver stochastic simulation, with varying fitness costs.
    Data is all stored in three files:
    final_data/resistant_modification_onePartner_haploidFitCosts_noSterile_total.csv
    final_data/resistant_modification_onePartner_haploidFitCosts_noSterile_NEWallele.csv
    final_data/resistant_modification_onePartner_haploidFitCosts_noSterile_adults.csv
    """
    num_partners = 5
    alleles = [['C', 'R', 'A'], ['V', 'W']] # a resistance allele is uncleavable
    
    # genotype 0 = cc vv, genotype 23 = ra, ww (wt resistant)
    s_c = []
    f_c = []

    file_name = "mutation_data/RM_manypercents_fivePartner_haploidFitCosts_noSterile_homoRelease"

    for clvr_cost in [0, 0.05, 0.1, 0.15]:
        for introFreq in [0.05, 0.1, 0.2]:
            intro = [[1, 0, 0.1], [0, 19, introFreq/2], [1, 19, introFreq/2]] # sex, genotype, frequency
            # 23 = ra ww, 19 = rr ww
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            hf_c = [[0, ['C', 'W'], 1.0, [['C']]], #savable via maternal carryover
                    # NOTE. 'C' in maternal carryover does NOT matter here, because MC = 0
                    [1, ['C', 'W'], 1.0, []],
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]] # haploid fitness cost

            run_label = f'clvr_cost_{clvr_cost}_intro_{introFreq}'
            run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro, f_c, hf_c, s_c, num_partners,
                            mut_flag = "NA",
                            run_label= run_label,
                            file_name= file_name, k=POP_MAX)
    return None

def RM_manypercents_fivePartner_diploid_fc_homo():
    """runs an arabidopsis cleaver stochastic simulation, with varying fitness costs.
    """
    num_partners = 5
    alleles = [['C', 'R', 'A'], ['V', 'W']] # a resistance allele is uncleavable
    
    # genotype 0 = cc vv, genotype 23 = ra, ww (wt resistant)
    s_c = []
    hf_c = [[0, ['C', 'W'], 1.0, [['C']]], #savable via maternal carryover
                    # NOTE. 'C' in maternal carryover does NOT matter here, because MC = 0
                    [1, ['C', 'W'], 1.0, []]]

    file_name = "mutation_data/RM_manypercents_fivePartner_diploidFitCosts_noSterile_homoRelease"

    for clvr_cost in [0, 0.05, 0.1, 0.15]:
        for introFreq in [0.05, 0.1, 0.2]:
            intro = [[1, 0, 0.1], [0, 19, introFreq/2], [1, 19, introFreq/2]] # sex, genotype, frequency
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            f_c = [[0, ['V'], clvr_cost, []], # haploid fitness cost
                   [1, ['V'], clvr_cost, []]] # haploid fitness cost

            run_label = f'clvr_cost_{clvr_cost}_intro_{introFreq}'
            run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro, f_c, hf_c, s_c, num_partners,
                            mut_flag = "NA",
                            run_label= run_label,
                            file_name= file_name, k=POP_MAX)
    return None

def RM_manypercents_fivePartner_diploidHaploid_fc_homo():
    """runs an arabidopsis cleaver stochastic simulation, with varying fitness costs.
    """
    num_partners = 5
    alleles = [['C', 'R', 'A'], ['V', 'W']] # a resistance allele is uncleavable
    
    # genotype 0 = cc vv, genotype 23 = ra, ww (wt resistant)
    s_c = []

    file_name = "mutation_data/RM_manypercents_fivePartner_diploidHaploidFitCosts_noSterile_homoRelease"

    for clvr_cost in [0, 0.05, 0.1, 0.15]:
        for introFreq in [0.05, 0.1, 0.2]:
            intro = [[1, 0, 0.1], [0, 19, introFreq/2], [1, 19, introFreq/2]] # sex, genotype, frequency
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            f_c = [[0, ['V'], clvr_cost, []], # haploid fitness cost
                   [1, ['V'], clvr_cost, []]] # haploid fitness cost

            hf_c = [[0, ['C', 'W'], 1.0, [['C']]], #savable via maternal carryover
                    # NOTE. 'C' in maternal carryover does NOT matter here, because MC = 0
                    [1, ['C', 'W'], 1.0, []],
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]]

            run_label = f'clvr_cost_{clvr_cost}_intro_{introFreq}'
            run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro, f_c, hf_c, s_c, num_partners,
                            mut_flag = "NA",
                            run_label= run_label,
                            file_name= file_name, k=POP_MAX)
    return None


def RM_manypercents_twentyPartner_haploid_fc():
    """runs an arabidopsis cleaver stochastic simulation, with varying fitness costs.
    Data is all stored in three files:
    final_data/resistant_modification_onePartner_haploidFitCosts_noSterile_total.csv
    final_data/resistant_modification_onePartner_haploidFitCosts_noSterile_NEWallele.csv
    final_data/resistant_modification_onePartner_haploidFitCosts_noSterile_adults.csv
    """
    num_partners = 20
    alleles = [['C', 'R', 'A'], ['V', 'W']] # a resistance allele is uncleavable
    
    # genotype 0 = cc vv, genotype 23 = ra, ww (wt resistant)
    s_c = []
    f_c = []

    file_name = "mutation_data/RM_manypercents_twentyPartner_haploidFitCosts_noSterile2"

    for clvr_cost in [0, 0.05, 0.1, 0.15]:
        for introFreq in [0.1, 0.2]:
            intro = [[1, 0, 0.1], [0, 19, introFreq/2], [1, 19, introFreq/2]] # sex, genotype, frequency
            # 23 = ra ww, 19 = rr ww
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            hf_c = [[0, ['C', 'W'], 1.0, [['C']]], #savable via maternal carryover
                    # NOTE. 'C' in maternal carryover does NOT matter here, because MC = 0
                    [1, ['C', 'W'], 1.0, []],
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]] # haploid fitness cost

            run_label = f'clvr_cost_{clvr_cost}_intro_{introFreq}'
            run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro, f_c, hf_c, s_c, num_partners,
                            mut_flag = "NA",
                            run_label= run_label,
                            file_name= file_name, k=POP_MAX)
    return None

def RM_manypercents_twentyPartner_diploid_fc():
    """runs an arabidopsis cleaver stochastic simulation, with varying fitness costs.
    """
    num_partners = 20
    alleles = [['C', 'R', 'A'], ['V', 'W']] # a resistance allele is uncleavable
    
    # genotype 0 = cc vv, genotype 23 = ra, ww (wt resistant)
    s_c = []
    hf_c = [[0, ['C', 'W'], 1.0, [['C']]], #savable via maternal carryover
                    # NOTE. 'C' in maternal carryover does NOT matter here, because MC = 0
                    [1, ['C', 'W'], 1.0, []]]

    file_name = "mutation_data/RM_manypercents_twentyPartner_diploidFitCosts_noSterile"

    for clvr_cost in [0, 0.05, 0.1, 0.15]:
        for introFreq in [0.1, 0.2]:
            intro = [[1, 0, 0.1], [0, 23, introFreq/2], [1, 23, introFreq/2]] # sex, genotype, frequency
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            f_c = [[0, ['V'], clvr_cost, []], # haploid fitness cost
                   [1, ['V'], clvr_cost, []]] # haploid fitness cost

            run_label = f'clvr_cost_{clvr_cost}_intro_{introFreq}'
            run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro, f_c, hf_c, s_c, num_partners,
                            mut_flag = "NA",
                            run_label= run_label,
                            file_name= file_name, k=POP_MAX)
    return None

def RM_manypercents_twentyPartner_diploidHaploid_fc():
    """runs an arabidopsis cleaver stochastic simulation, with varying fitness costs.
    """
    num_partners = 20
    alleles = [['C', 'R', 'A'], ['V', 'W']] # a resistance allele is uncleavable
    
    # genotype 0 = cc vv, genotype 23 = ra, ww (wt resistant)
    s_c = []

    file_name = "mutation_data/RM_manypercents_twentyPartner_diploidHaploidFitCosts_noSterile"

    for clvr_cost in [0, 0.05, 0.1, 0.15]:
        for introFreq in [0.1, 0.2]:
            intro = [[1, 0, 0.1], [0, 23, introFreq/2], [1, 23, introFreq/2]] # sex, genotype, frequency
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            f_c = [[0, ['V'], clvr_cost, []], # haploid fitness cost
                   [1, ['V'], clvr_cost, []]] # haploid fitness cost

            hf_c = [[0, ['C', 'W'], 1.0, [['C']]], #savable via maternal carryover
                    # NOTE. 'C' in maternal carryover does NOT matter here, because MC = 0
                    [1, ['C', 'W'], 1.0, []],
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]]

            run_label = f'clvr_cost_{clvr_cost}_intro_{introFreq}'
            run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro, f_c, hf_c, s_c, num_partners,
                            mut_flag = "NA",
                            run_label= run_label,
                            file_name= file_name, k=POP_MAX)
    return None

def RM_manypercents_twentyPartner_haploid_fc_homo():
    """runs an arabidopsis cleaver stochastic simulation, with varying fitness costs.
    Data is all stored in three files:
    final_data/resistant_modification_onePartner_haploidFitCosts_noSterile_total.csv
    final_data/resistant_modification_onePartner_haploidFitCosts_noSterile_NEWallele.csv
    final_data/resistant_modification_onePartner_haploidFitCosts_noSterile_adults.csv
    """
    num_partners = 20
    alleles = [['C', 'R', 'A'], ['V', 'W']] # a resistance allele is uncleavable
    
    # genotype 0 = cc vv, genotype 23 = ra, ww (wt resistant)
    s_c = []
    f_c = []

    file_name = "mutation_data/RM_manypercents_twentyPartner_haploidFitCosts_noSterile_homoRelease"

    for clvr_cost in [0, 0.05, 0.1, 0.15]:
        for introFreq in [0.05, 0.1, 0.2]:
            intro = [[1, 0, 0.1], [0, 19, introFreq/2], [1, 19, introFreq/2]] # sex, genotype, frequency
            # 23 = ra ww, 19 = rr ww
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            hf_c = [[0, ['C', 'W'], 1.0, [['C']]], #savable via maternal carryover
                    # NOTE. 'C' in maternal carryover does NOT matter here, because MC = 0
                    [1, ['C', 'W'], 1.0, []],
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]] # haploid fitness cost

            run_label = f'clvr_cost_{clvr_cost}_intro_{introFreq}'
            run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro, f_c, hf_c, s_c, num_partners,
                            mut_flag = "NA",
                            run_label= run_label,
                            file_name= file_name, k=POP_MAX)
    return None

def RM_manypercents_twentyPartner_diploid_fc_homo():
    """runs an arabidopsis cleaver stochastic simulation, with varying fitness costs.
    """
    num_partners = 20
    alleles = [['C', 'R', 'A'], ['V', 'W']] # a resistance allele is uncleavable
    
    # genotype 0 = cc vv, genotype 23 = ra, ww (wt resistant)
    s_c = []
    hf_c = [[0, ['C', 'W'], 1.0, [['C']]], #savable via maternal carryover
                    # NOTE. 'C' in maternal carryover does NOT matter here, because MC = 0
                    [1, ['C', 'W'], 1.0, []]]

    file_name = "mutation_data/RM_manypercents_twentyPartner_diploidFitCosts_noSterile_homoRelease"

    for clvr_cost in [0, 0.05, 0.1, 0.15]:
        for introFreq in [0.05, 0.1, 0.2]:
            intro = [[1, 0, 0.1], [0, 19, introFreq/2], [1, 19, introFreq/2]] # sex, genotype, frequency
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            f_c = [[0, ['V'], clvr_cost, []], # haploid fitness cost
                   [1, ['V'], clvr_cost, []]] # haploid fitness cost

            run_label = f'clvr_cost_{clvr_cost}_intro_{introFreq}'
            run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro, f_c, hf_c, s_c, num_partners,
                            mut_flag = "NA",
                            run_label= run_label,
                            file_name= file_name, k=POP_MAX)
    return None

def RM_manypercents_twentyPartner_diploidHaploid_fc_homo():
    """runs an arabidopsis cleaver stochastic simulation, with varying fitness costs.
    """
    num_partners = 20
    alleles = [['C', 'R', 'A'], ['V', 'W']] # a resistance allele is uncleavable
    
    # genotype 0 = cc vv, genotype 23 = ra, ww (wt resistant)
    s_c = []

    file_name = "mutation_data/RM_manypercents_twentyPartner_diploidHaploidFitCosts_noSterile_homoRelease"

    for clvr_cost in [0, 0.05, 0.1, 0.15]:
        for introFreq in [0.05, 0.1, 0.2]:
            intro = [[1, 0, 0.1], [0, 19, introFreq/2], [1, 19, introFreq/2]] # sex, genotype, frequency
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            f_c = [[0, ['V'], clvr_cost, []], # haploid fitness cost
                   [1, ['V'], clvr_cost, []]] # haploid fitness cost

            hf_c = [[0, ['C', 'W'], 1.0, [['C']]], #savable via maternal carryover
                    # NOTE. 'C' in maternal carryover does NOT matter here, because MC = 0
                    [1, ['C', 'W'], 1.0, []],
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]]

            run_label = f'clvr_cost_{clvr_cost}_intro_{introFreq}'
            run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro, f_c, hf_c, s_c, num_partners,
                            mut_flag = "NA",
                            run_label= run_label,
                            file_name= file_name, k=POP_MAX)
    return None

# endregion RMmod, haploid v. diploid v. haplo-diploid

# region haplo

# region modification NO MUT., FC types

def explore_fitness_costs_diploid_dom():
    """runs an arabidopsis cleaver stochastic simulation, with varying fitness costs.
    Data is all stored in three files:
    gamete_tracking/fitCosts_noSterile_adults.csv,
    gamete_tracking/fitCosts_noSterile_alleles.csv, and 
    gamete_tracking/fitCosts_noSterile_total.csv"""
    num_partners = 1
    num_gens = 100
    num_reps = 10
    intro = [[1, 0, 0.1]] # sex, genotype, frequency
    hf_c = [[0, ['C', 'W'], 1.0, ['V']], # haploid fitness costs are such that CW is lethal,
                                         # but can be rescued by maternal carryover
            [1, ['C', 'W'], 1.0, []]]
    s_c = []
    pop_max = 10000

    file_name = "mutation_data2/fitCosts_noSterile_dom"

    for clvr_cost in [0, 0.05, 0.1, 0.15]:
        # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
        f_c = [[0, ['V'], clvr_cost, []],
               [1, ['V'], clvr_cost, []]]

        run_label = f'clvr_cost_{clvr_cost}_intro_{intro[0][2]}'
        run_stochastic_sim(ALLELES, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=POP_MAX)

def explore_fitness_costs_diploid_rec():
    """runs an arabidopsis cleaver stochastic simulation, with varying fitness costs.
    Data is all stored in three files:
    gamete_tracking/fitCosts_noSterile_adults.csv,
    gamete_tracking/fitCosts_noSterile_alleles.csv, and 
    gamete_tracking/fitCosts_noSterile_total.csv"""
    num_partners = 1
    num_gens = 100
    num_reps = 10
    intro = [[1, 0, 0.1]] # sex, genotype, frequency
    hf_c = [[0, ['C', 'W'], 1.0, ['V']], # haploid fitness costs are such that CW is lethal,
                                         # but can be rescued by maternal carryover
            [1, ['C', 'W'], 1.0, []]]
    s_c = []
    pop_max = 10000

    file_name = "mutation_data2/fitCosts_noSterile_rec"

    for clvr_cost in [0, 0.05, 0.1, 0.15]:
        # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
        f_c = [[0, ['V', 'V'], clvr_cost, []],
               [1, ['V', 'V'], clvr_cost, []]]

        run_label = f'clvr_cost_{clvr_cost}_intro_{intro[0][2]}'
        run_stochastic_sim(ALLELES, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=POP_MAX)

def explore_fitness_costs_diploid_add():
    """runs an arabidopsis cleaver stochastic simulation, with varying fitness costs.
    Data is all stored in three files:
    gamete_tracking/fitCosts_noSterile_adults.csv,
    gamete_tracking/fitCosts_noSterile_alleles.csv, and 
    gamete_tracking/fitCosts_noSterile_total.csv"""
    num_partners = 1
    num_gens = 100
    num_reps = 10
    intro = [[1, 0, 0.1]] # sex, genotype, frequency
    hf_c = [[0, ['C', 'W'], 1.0, ['V']], # haploid fitness costs are such that CW is lethal,
                                         # but can be rescued by maternal carryover
            [1, ['C', 'W'], 1.0, []]]
    s_c = []
    pop_max = 10000

    file_name = "mutation_data2/fitCosts_noSterile_rec"

    for clvr_cost in [0, 0.05, 0.1, 0.15]:
        # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
        f_c = [[0, ['V'], clvr_cost, []],
               [1, ['V'], clvr_cost, []]]

        run_label = f'clvr_cost_{clvr_cost}_intro_{intro[0][2]}'
        run_stochastic_sim(ALLELES, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=POP_MAX, FC_flag="additive")


############## intro 20 ##################

##########################################
############ female sterility ############
##########################################

#region female sterility

# region resistant suppression / resistant to cleavage

def resistant_suppression_femalesterile_onePartner_MC():
    """runs simulations for multiple maternal carryovers and various haploid 
    fitness costs, for mating 1 female to 1 male"""
    num_partners = 1
    alleles = [['C', 'R', 'A'], ['V', 'W']] # a resistance allele is uncleavable
    intro = [[1, 0, 0.1], [0, 23, 0.025], [1, 23, 0.025]] # sex, genotype, frequency
    # genotype 0 = cc vv, genotype 23 = ra, ww (wt resistant)
    s_c = [[0, ['V', 'V'], 1.0]] #sex, alleles, fert_cost - females homozygous sterile
    f_c = []

    file_name = "mutation_data/resistant_suppression_onePartner_femSterile"

    for maternal_carryover in [0, 0.1, 0.2, 0.3]:
        for clvr_cost in [0, 0.05, 0.1, 0.15, 0.2]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            hf_c = [[0, ['C', 'W'], 1.0, [['V']]],
                    [1, ['C', 'W'], 1.0, []],
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]] # haploid fitness cost

            run_label = f'mc_prob_{maternal_carryover}_FC_{clvr_cost}'
            run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=POP_MAX, mc_prob=maternal_carryover)

    return None

def RS_1percent_femalesterile_onePartner_MC():
    """runs simulations for multiple maternal carryovers and various haploid 
    fitness costs, for mating 1 female to 1 male"""
    num_partners = 1
    alleles = [['C', 'R', 'A'], ['V', 'W']] # a resistance allele is uncleavable
    intro = [[1, 0, 0.1], [0, 23, 0.005], [1, 23, 0.005]] # sex, genotype, frequency
    # genotype 0 = cc vv, genotype 23 = ra, ww (wt resistant)
    s_c = [[0, ['V', 'V'], 1.0]] #sex, alleles, fert_cost - females homozygous sterile
    f_c = []

    file_name = "mutation_data/RS_1percent_onePartner_femSterile"

    for maternal_carryover in [0, 0.1, 0.2, 0.3]:
        for clvr_cost in [0, 0.05, 0.1, 0.15, 0.2]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            hf_c = [[0, ['C', 'W'], 1.0, [['V']]],
                    [1, ['C', 'W'], 1.0, []],
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]] # haploid fitness cost

            run_label = f'mc_prob_{maternal_carryover}_FC_{clvr_cost}'
            run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=POP_MAX, mc_prob=maternal_carryover)

    return None

def RS_05percent_femalesterile_onePartner_MC():
    """runs simulations for multiple maternal carryovers and various haploid 
    fitness costs, for mating 1 female to 1 male"""
    num_partners = 1
    alleles = [['C', 'R', 'A'], ['V', 'W']] # a resistance allele is uncleavable
    intro = [[1, 0, 0.1], [0, 23, 0.0025], [1, 23, 0.0025]] # sex, genotype, frequency
    # genotype 0 = cc vv, genotype 23 = ra, ww (wt resistant)
    s_c = [[0, ['V', 'V'], 1.0]] #sex, alleles, fert_cost - females homozygous sterile
    f_c = []

    file_name = "mutation_data/RS_05percent_onePartner_femSterile"

    for maternal_carryover in [0, 0.1, 0.2, 0.3]:
        for clvr_cost in [0, 0.05, 0.1, 0.15, 0.2]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            hf_c = [[0, ['C', 'W'], 1.0, [['V']]],
                    [1, ['C', 'W'], 1.0, []],
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]] # haploid fitness cost

            run_label = f'mc_prob_{maternal_carryover}_FC_{clvr_cost}'
            run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=POP_MAX, mc_prob=maternal_carryover)

    return None

def RS_01percent_femalesterile_onePartner_MC():
    """runs simulations for multiple maternal carryovers and various haploid 
    fitness costs, for mating 1 female to 1 male
    
    for .1%, or 10^-3"""
    num_partners = 1
    alleles = [['C', 'R', 'A'], ['V', 'W']] # a resistance allele is uncleavable
    intro = [[1, 0, 0.1], [0, 23, 0.0005], [1, 23, 0.0005]] # sex, genotype, frequency
    # genotype 0 = cc vv, genotype 23 = ra, ww (wt resistant)
    s_c = [[0, ['V', 'V'], 1.0]] #sex, alleles, fert_cost - females homozygous sterile
    f_c = []

    file_name = "mutation_data/RS_01percent_onePartner_femSterile"

    for maternal_carryover in [0, 0.1, 0.2, 0.3]:
        for clvr_cost in [0, 0.05, 0.1, 0.15, 0.2]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            hf_c = [[0, ['C', 'W'], 1.0, [['V']]],
                    [1, ['C', 'W'], 1.0, []],
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]] # haploid fitness cost

            run_label = f'mc_prob_{maternal_carryover}_FC_{clvr_cost}'
            run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=POP_MAX, mc_prob=maternal_carryover)

    return None

def RS_o01percent_femalesterile_onePartner_MC():
    """runs simulations for multiple maternal carryovers and various haploid 
    fitness costs, for mating 1 female to 1 male
    
    for 10^-4"""
    num_partners = 1
    alleles = [['C', 'R', 'A'], ['V', 'W']] # a resistance allele is uncleavable
    intro = [[1, 0, 0.1], [0, 23, 0.00005], [1, 23, 0.00005]] # sex, genotype, frequency
    # genotype 0 = cc vv, genotype 23 = ra, ww (wt resistant)
    s_c = [[0, ['V', 'V'], 1.0]] #sex, alleles, fert_cost - females homozygous sterile
    f_c = []

    file_name = "mutation_data/RS_o01percent_onePartner_femSterile"

    K = 1000000 # 1 000 000

    for maternal_carryover in [0, 0.3]:
        for clvr_cost in [0, 0.05, 0.1, 0.15]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            hf_c = [[0, ['C', 'W'], 1.0, [['V']]],
                    [1, ['C', 'W'], 1.0, []],
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]] # haploid fitness cost

            run_label = f'mc_prob_{maternal_carryover}_FC_{clvr_cost}'
            run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=K, mc_prob=maternal_carryover)

    return None

def RS_o001percent_femalesterile_onePartner_MC():
    """runs simulations for multiple maternal carryovers and various haploid 
    fitness costs, for mating 1 female to 1 male"""
    num_partners = 1
    alleles = [['C', 'R', 'A'], ['V', 'W']] # a resistance allele is uncleavable
    intro = [[1, 0, 0.1], [0, 23, 0.000005], [1, 23, 0.000005]] # sex, genotype, frequency
    # genotype 0 = cc vv, genotype 23 = ra, ww (wt resistant)
    s_c = [[0, ['V', 'V'], 1.0]] #sex, alleles, fert_cost - females homozygous sterile
    f_c = []

    file_name = "mutation_data/RS_o001percent_onePartner_femSterile"

    K = 10000000 # 10 000 000

    for maternal_carryover in [0, 0.3]:
        for clvr_cost in [0, 0.05, 0.1, 0.15]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            hf_c = [[0, ['C', 'W'], 1.0, [['V']]],
                    [1, ['C', 'W'], 1.0, []],
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]] # haploid fitness cost

            run_label = f'mc_prob_{maternal_carryover}_FC_{clvr_cost}'
            run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=K, mc_prob=maternal_carryover)

    return None

def RS_o0001percent_femalesterile_onePartner_MC():
    """runs simulations for multiple maternal carryovers and various haploid 
    fitness costs, for mating 1 female to 1 male"""
    num_partners = 1
    alleles = [['C', 'R', 'A'], ['V', 'W']] # a resistance allele is uncleavable
    intro = [[1, 0, 0.1], [0, 23, 0.0000005], [1, 23, 0.0000005]] # sex, genotype, frequency
    # genotype 0 = cc vv, genotype 23 = ra, ww (wt resistant)
    s_c = [[0, ['V', 'V'], 1.0]] #sex, alleles, fert_cost - females homozygous sterile
    f_c = []

    file_name = "mutation_data/RS_o0001percent_onePartner_femSterile"

    K = 100000000 # 100 000 000

    for maternal_carryover in [0, 0.3]:
        for clvr_cost in [0, 0.05, 0.1, 0.15]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            hf_c = [[0, ['C', 'W'], 1.0, [['V']]],
                    [1, ['C', 'W'], 1.0, []],
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]] # haploid fitness cost

            run_label = f'mc_prob_{maternal_carryover}_FC_{clvr_cost}'
            run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=K, mc_prob=maternal_carryover)

    return None


def int20_resistant_suppression_femalesterile_onePartner_MC():
    """runs simulations for multiple maternal carryovers and various haploid 
    fitness costs, for mating 1 female to 1 male"""
    num_partners = 1
    alleles = [['C', 'R', 'A'], ['V', 'W']] # a resistance allele is uncleavable
    intro = [[1, 0, 0.2], [0, 23, 0.025], [1, 23, 0.025]] # sex, genotype, frequency
    # genotype 0 = cc vv, genotype 23 = ra, ww (wt resistant)
    s_c = [[0, ['V', 'V'], 1.0]] #sex, alleles, fert_cost - females homozygous sterile
    f_c = []

    file_name = "mutation_data/int20_resistant_suppression_onePartner_femSterile"

    for maternal_carryover in [0, 0.1, 0.2, 0.3]:
        for clvr_cost in [0, 0.05, 0.1, 0.15, 0.2]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            hf_c = [[0, ['C', 'W'], 1.0, [['V']]],
                    [1, ['C', 'W'], 1.0, []],
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]] # haploid fitness cost

            run_label = f'mc_prob_{maternal_carryover}_FC_{clvr_cost}'
            run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=POP_MAX, mc_prob=maternal_carryover)

    return None

# endregion resistant suppression / resistant to cleavage

# region clvr loss of function

def lof_clvr_suppression_femalesterile_onePartner_MC():
    """runs simulations for multiple maternal carryovers and various haploid 
    fitness costs, for mating 1 female to 1 male"""
    num_partners = 1
    alleles = [['C', 'A'], ['V', 'R', 'W']] # a resistance allele is uncleavable
    intro = [[1, 0, 0.08], [1, 13, 0.02]] # sex, genotype, frequency
    # genotype 0 = cc vv, 13 = [['C', 'A'], ['R', 'R']]
    s_c = [[0, ['V', 'V'], 1.0],
           [0, ['V', 'R'], 1.0],
           [0, ['R', 'V'], 1.0],
           [0, ['R', 'R'], 1.0]] #sex, alleles, fert_cost - females homozygous sterile
    f_c = []

    file_name = "mutation_data/lof_clvr_suppression_onePartner_femSterile"

    for maternal_carryover in [0, 0.1, 0.2, 0.3]:
        for clvr_cost in [0, 0.05, 0.1, 0.15, 0.2]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            hf_c = [[0, ['C', 'W'], 1.0, [['V']]],
                    [1, ['C', 'W'], 1.0, []],
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]] # haploid fitness cost

            run_label = f'mc_prob_{maternal_carryover}_FC_{clvr_cost}'
            run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=POP_MAX, mc_prob=maternal_carryover)

    return None

def lof_clvr_suppression_femalesterile_fivePartner_MC():
    """runs simulations for multiple maternal carryovers and various haploid 
    fitness costs, for mating 1 female to 5 males"""
    num_partners = 5
    alleles = [['C', 'A'], ['V', 'R', 'W']] # a resistance allele is uncleavable
    intro = [[1, 0, 0.08], [1, 13, 0.02]] # sex, genotype, frequency
    # genotype 0 = cc vv, 13 = [['C', 'A'], ['R', 'R']]
    s_c = [[0, ['V', 'V'], 1.0],
           [0, ['V', 'R'], 1.0],
           [0, ['R', 'V'], 1.0],
           [0, ['R', 'R'], 1.0]] #sex, alleles, fert_cost - females homozygous sterile
    f_c = []

    file_name = "mutation_data/lof_clvr_suppression_fivePartner_femSterile"

    for maternal_carryover in [0, 0.1, 0.2, 0.3]:
        for clvr_cost in [0, 0.05, 0.1, 0.15, 0.2]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            hf_c = [[0, ['C', 'W'], 1.0, [['V']]],
                    [1, ['C', 'W'], 1.0, []],
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]] # haploid fitness cost

            run_label = f'mc_prob_{maternal_carryover}_FC_{clvr_cost}'
            run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=POP_MAX, mc_prob=maternal_carryover)

    return None

def lof_clvr_suppression_femalesterile_twentyPartner_MC():
    """runs simulations for multiple maternal carryovers and various haploid 
    fitness costs, for mating 1 female to 20 males"""
    num_partners = 20
    alleles = [['C', 'A'], ['V', 'R', 'W']] # a resistance allele is uncleavable
    intro = [[1, 0, 0.08], [1, 13, 0.02]] # sex, genotype, frequency
    # genotype 0 = cc vv, 13 = [['C', 'A'], ['R', 'R']]
    s_c = [[0, ['V', 'V'], 1.0],
           [0, ['V', 'R'], 1.0],
           [0, ['R', 'V'], 1.0],
           [0, ['R', 'R'], 1.0]] #sex, alleles, fert_cost - females homozygous sterile
    f_c = []

    file_name = "mutation_data/lof_clvr_suppression_twentyPartner_femSterile"

    for maternal_carryover in [0, 0.1, 0.2, 0.3]:
        for clvr_cost in [0, 0.05, 0.1, 0.15, 0.2]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            hf_c = [[0, ['C', 'W'], 1.0, [['V']]],
                    [1, ['C', 'W'], 1.0, []],
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]] # haploid fitness cost

            run_label = f'mc_prob_{maternal_carryover}_FC_{clvr_cost}'
            run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=POP_MAX, mc_prob=maternal_carryover)

    return None

# endregion clvr loss of function

# region clvr LOF - new IF's

def lof_clvr_if1_suppression_femalesterile_onePartner_MC():
    """runs simulations for multiple maternal carryovers and various haploid 
    fitness costs, for mating 1 female to 1 male"""
    num_partners = 1
    alleles = [['C', 'A'], ['V', 'R', 'W']] # a resistance allele is uncleavable
    intro = [[1, 0, 0.099], [1, 13, 0.001]] # sex, genotype, frequency
    # genotype 0 = cc vv, 13 = [['C', 'A'], ['R', 'R']]
    s_c = [[0, ['V', 'V'], 1.0],
           [0, ['V', 'R'], 1.0],
           [0, ['R', 'V'], 1.0],
           [0, ['R', 'R'], 1.0]] #sex, alleles, fert_cost - females homozygous sterile
    f_c = []

    file_name = "mutation_data2/lof_clvr_if1_suppression_onePartner_femSterile"

    for maternal_carryover in [0, 0.3]:
        for clvr_cost in [0, 0.05, 0.1, 0.15, 0.2]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            hf_c = [[0, ['C', 'W'], 1.0, [['V']]],
                    [1, ['C', 'W'], 1.0, []],
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]] # haploid fitness cost

            run_label = f'mc_prob_{maternal_carryover}_FC_{clvr_cost}'
            run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=POP_MAX, mc_prob=maternal_carryover)

    return None

def lof_clvr_if2_suppression_femalesterile_onePartner_MC():
    """runs simulations for multiple maternal carryovers and various haploid 
    fitness costs, for mating 1 female to 1 male"""
    num_partners = 1
    alleles = [['C', 'A'], ['V', 'R', 'W']] # a resistance allele is uncleavable
    intro = [[1, 0, 0.098], [1, 13, 0.002]] # sex, genotype, frequency
    # genotype 0 = cc vv, 13 = [['C', 'A'], ['R', 'R']]
    s_c = [[0, ['V', 'V'], 1.0],
           [0, ['V', 'R'], 1.0],
           [0, ['R', 'V'], 1.0],
           [0, ['R', 'R'], 1.0]] #sex, alleles, fert_cost - females homozygous sterile
    f_c = []

    file_name = "mutation_data2/lof_clvr_if2_suppression_onePartner_femSterile"

    for maternal_carryover in [0, 0.3]:
        for clvr_cost in [0, 0.05, 0.1, 0.15, 0.2]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            hf_c = [[0, ['C', 'W'], 1.0, [['V']]],
                    [1, ['C', 'W'], 1.0, []],
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]] # haploid fitness cost

            run_label = f'mc_prob_{maternal_carryover}_FC_{clvr_cost}'
            run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=POP_MAX, mc_prob=maternal_carryover)

    return None

def lof_clvr_if5_suppression_femalesterile_onePartner_MC():
    """runs simulations for multiple maternal carryovers and various haploid 
    fitness costs, for mating 1 female to 1 male"""
    num_partners = 1
    alleles = [['C', 'A'], ['V', 'R', 'W']] # a resistance allele is uncleavable
    intro = [[1, 0, 0.095], [1, 13, 0.005]] # sex, genotype, frequency
    # genotype 0 = cc vv, 13 = [['C', 'A'], ['R', 'R']]
    s_c = [[0, ['V', 'V'], 1.0],
           [0, ['V', 'R'], 1.0],
           [0, ['R', 'V'], 1.0],
           [0, ['R', 'R'], 1.0]] #sex, alleles, fert_cost - females homozygous sterile
    f_c = []

    file_name = "mutation_data2/lof_clvr_if5_suppression_onePartner_femSterile"

    for maternal_carryover in [0, 0.3]:
        for clvr_cost in [0, 0.05, 0.1, 0.15, 0.2]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            hf_c = [[0, ['C', 'W'], 1.0, [['V']]],
                    [1, ['C', 'W'], 1.0, []],
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]] # haploid fitness cost

            run_label = f'mc_prob_{maternal_carryover}_FC_{clvr_cost}'
            run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=POP_MAX, mc_prob=maternal_carryover)

    return None

def lof_clvr_if1_suppression_femalesterile_fivePartner_MC():
    """runs simulations for multiple maternal carryovers and various haploid 
    fitness costs, for mating 1 female to 1 male"""
    num_partners = 5
    alleles = [['C', 'A'], ['V', 'R', 'W']] # a resistance allele is uncleavable
    intro = [[1, 0, 0.099], [1, 13, 0.001]] # sex, genotype, frequency
    # genotype 0 = cc vv, 13 = [['C', 'A'], ['R', 'R']]
    s_c = [[0, ['V', 'V'], 1.0],
           [0, ['V', 'R'], 1.0],
           [0, ['R', 'V'], 1.0],
           [0, ['R', 'R'], 1.0]] #sex, alleles, fert_cost - females homozygous sterile
    f_c = []

    file_name = "mutation_data2/lof_clvr_if1_suppression_fivePartner_femSterile"

    for maternal_carryover in [0, 0.3]:
        for clvr_cost in [0, 0.05, 0.1, 0.15, 0.2]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            hf_c = [[0, ['C', 'W'], 1.0, [['V']]],
                    [1, ['C', 'W'], 1.0, []],
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]] # haploid fitness cost

            run_label = f'mc_prob_{maternal_carryover}_FC_{clvr_cost}'
            run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=POP_MAX, mc_prob=maternal_carryover)

    return None

def lof_clvr_if2_suppression_femalesterile_fivePartner_MC():
    """runs simulations for multiple maternal carryovers and various haploid 
    fitness costs, for mating 1 female to 1 male"""
    num_partners = 5
    alleles = [['C', 'A'], ['V', 'R', 'W']] # a resistance allele is uncleavable
    intro = [[1, 0, 0.098], [1, 13, 0.002]] # sex, genotype, frequency
    # genotype 0 = cc vv, 13 = [['C', 'A'], ['R', 'R']]
    s_c = [[0, ['V', 'V'], 1.0],
           [0, ['V', 'R'], 1.0],
           [0, ['R', 'V'], 1.0],
           [0, ['R', 'R'], 1.0]] #sex, alleles, fert_cost - females homozygous sterile
    f_c = []

    file_name = "mutation_data2/lof_clvr_if2_suppression_fivePartner_femSterile"

    for maternal_carryover in [0, 0.3]:
        for clvr_cost in [0, 0.05, 0.1, 0.15, 0.2]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            hf_c = [[0, ['C', 'W'], 1.0, [['V']]],
                    [1, ['C', 'W'], 1.0, []],
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]] # haploid fitness cost

            run_label = f'mc_prob_{maternal_carryover}_FC_{clvr_cost}'
            run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=POP_MAX, mc_prob=maternal_carryover)

    return None

def lof_clvr_if5_suppression_femalesterile_fivePartner_MC():
    """runs simulations for multiple maternal carryovers and various haploid 
    fitness costs, for mating 1 female to 1 male"""
    num_partners = 5
    alleles = [['C', 'A'], ['V', 'R', 'W']] # a resistance allele is uncleavable
    intro = [[1, 0, 0.095], [1, 13, 0.005]] # sex, genotype, frequency
    # genotype 0 = cc vv, 13 = [['C', 'A'], ['R', 'R']]
    s_c = [[0, ['V', 'V'], 1.0],
           [0, ['V', 'R'], 1.0],
           [0, ['R', 'V'], 1.0],
           [0, ['R', 'R'], 1.0]] #sex, alleles, fert_cost - females homozygous sterile
    f_c = []

    file_name = "mutation_data2/lof_clvr_if5_suppression_fivePartner_femSterile"

    for maternal_carryover in [0, 0.3]:
        for clvr_cost in [0, 0.05, 0.1, 0.15, 0.2]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            hf_c = [[0, ['C', 'W'], 1.0, [['V']]],
                    [1, ['C', 'W'], 1.0, []],
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]] # haploid fitness cost

            run_label = f'mc_prob_{maternal_carryover}_FC_{clvr_cost}'
            run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=POP_MAX, mc_prob=maternal_carryover)

    return None

def lof_clvr_if1_suppression_femalesterile_twentyPartner_MC():
    """runs simulations for multiple maternal carryovers and various haploid 
    fitness costs, for mating 1 female to 1 male"""
    num_partners = 20
    alleles = [['C', 'A'], ['V', 'R', 'W']] # a resistance allele is uncleavable
    intro = [[1, 0, 0.099], [1, 13, 0.001]] # sex, genotype, frequency
    # genotype 0 = cc vv, 13 = [['C', 'A'], ['R', 'R']]
    s_c = [[0, ['V', 'V'], 1.0],
           [0, ['V', 'R'], 1.0],
           [0, ['R', 'V'], 1.0],
           [0, ['R', 'R'], 1.0]] #sex, alleles, fert_cost - females homozygous sterile
    f_c = []

    file_name = "mutation_data2/lof_clvr_if1_suppression_twentyPartner_femSterile"

    for maternal_carryover in [0, 0.3]:
        for clvr_cost in [0, 0.05, 0.1, 0.15, 0.2]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            hf_c = [[0, ['C', 'W'], 1.0, [['V']]],
                    [1, ['C', 'W'], 1.0, []],
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]] # haploid fitness cost

            run_label = f'mc_prob_{maternal_carryover}_FC_{clvr_cost}'
            run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=POP_MAX, mc_prob=maternal_carryover)

    return None

def lof_clvr_if2_suppression_femalesterile_twentyPartner_MC():
    """runs simulations for multiple maternal carryovers and various haploid 
    fitness costs, for mating 1 female to 1 male"""
    num_partners = 20
    alleles = [['C', 'A'], ['V', 'R', 'W']] # a resistance allele is uncleavable
    intro = [[1, 0, 0.098], [1, 13, 0.002]] # sex, genotype, frequency
    # genotype 0 = cc vv, 13 = [['C', 'A'], ['R', 'R']]
    s_c = [[0, ['V', 'V'], 1.0],
           [0, ['V', 'R'], 1.0],
           [0, ['R', 'V'], 1.0],
           [0, ['R', 'R'], 1.0]] #sex, alleles, fert_cost - females homozygous sterile
    f_c = []

    file_name = "mutation_data2/lof_clvr_if2_suppression_twentyPartner_femSterile"

    for maternal_carryover in [0, 0.3]:
        for clvr_cost in [0, 0.05, 0.1, 0.15, 0.2]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            hf_c = [[0, ['C', 'W'], 1.0, [['V']]],
                    [1, ['C', 'W'], 1.0, []],
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]] # haploid fitness cost

            run_label = f'mc_prob_{maternal_carryover}_FC_{clvr_cost}'
            run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=POP_MAX, mc_prob=maternal_carryover)

    return None

def lof_clvr_if5_suppression_femalesterile_twentyPartner_MC():
    """runs simulations for multiple maternal carryovers and various haploid 
    fitness costs, for mating 1 female to 1 male"""
    num_partners = 20
    alleles = [['C', 'A'], ['V', 'R', 'W']] # a resistance allele is uncleavable
    intro = [[1, 0, 0.095], [1, 13, 0.005]] # sex, genotype, frequency
    # genotype 0 = cc vv, 13 = [['C', 'A'], ['R', 'R']]
    s_c = [[0, ['V', 'V'], 1.0],
           [0, ['V', 'R'], 1.0],
           [0, ['R', 'V'], 1.0],
           [0, ['R', 'R'], 1.0]] #sex, alleles, fert_cost - females homozygous sterile
    f_c = []

    file_name = "mutation_data2/lof_clvr_if5_suppression_twentyPartner_femSterile"

    for maternal_carryover in [0, 0.3]:
        for clvr_cost in [0, 0.05, 0.1, 0.15, 0.2]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            hf_c = [[0, ['C', 'W'], 1.0, [['V']]],
                    [1, ['C', 'W'], 1.0, []],
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]] # haploid fitness cost

            run_label = f'mc_prob_{maternal_carryover}_FC_{clvr_cost}'
            run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=POP_MAX, mc_prob=maternal_carryover)

    return None


# endregion clvr LOF - new IF's... the only thing different is maternal carryover. So.

def lof_clvr_if1_suppression_femalesterile_onePartner_MC_2():
    """runs simulations for multiple maternal carryovers and various haploid 
    fitness costs, for mating 1 female to 1 male"""
    num_partners = 1
    alleles = [['C', 'A'], ['V', 'R', 'W']] # a resistance allele is uncleavable
    intro = [[1, 0, 0.099], [1, 13, 0.001]] # sex, genotype, frequency
    # genotype 0 = cc vv, 13 = [['C', 'A'], ['R', 'R']]
    s_c = [[0, ['V', 'V'], 1.0],
           [0, ['V', 'R'], 1.0],
           [0, ['R', 'V'], 1.0],
           [0, ['R', 'R'], 1.0]] #sex, alleles, fert_cost - females homozygous sterile
    f_c = []

    file_name = "mutation_data2/lof_clvr_if1_suppression_onePartner_femSterile_2"

    for maternal_carryover in [0, 0.3]:
        for clvr_cost in [0, 0.05, 0.1, 0.15, 0.2]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            hf_c = [[0, ['C', 'W'], 1.0, [['V'], ['R']]],
                    [1, ['C', 'W'], 1.0, []],
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]] # haploid fitness cost

            run_label = f'mc_prob_{maternal_carryover}_FC_{clvr_cost}'
            run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=POP_MAX, mc_prob=maternal_carryover)

    return None

def lof_clvr_if1_suppression_femalesterile_fivePartner_MC_2():
    """runs simulations for multiple maternal carryovers and various haploid 
    fitness costs, for mating 1 female to 1 male"""
    num_partners = 5
    alleles = [['C', 'A'], ['V', 'R', 'W']] # a resistance allele is uncleavable
    intro = [[1, 0, 0.099], [1, 13, 0.001]] # sex, genotype, frequency
    # genotype 0 = cc vv, 13 = [['C', 'A'], ['R', 'R']]
    s_c = [[0, ['V', 'V'], 1.0],
           [0, ['V', 'R'], 1.0],
           [0, ['R', 'V'], 1.0],
           [0, ['R', 'R'], 1.0]] #sex, alleles, fert_cost - females homozygous sterile
    f_c = []

    file_name = "mutation_data2/lof_clvr_if1_suppression_fivePartner_femSterile_2"

    for maternal_carryover in [0, 0.3]:
        for clvr_cost in [0, 0.05, 0.1, 0.15, 0.2]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            hf_c = [[0, ['C', 'W'], 1.0, [['V'], ['R']]],
                    [1, ['C', 'W'], 1.0, []],
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]] # haploid fitness cost

            run_label = f'mc_prob_{maternal_carryover}_FC_{clvr_cost}'
            run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=POP_MAX, mc_prob=maternal_carryover)

    return None

def lof_clvr_if1_suppression_femalesterile_twentyPartner_MC_2():
    """runs simulations for multiple maternal carryovers and various haploid 
    fitness costs, for mating 1 female to 1 male"""
    num_partners = 20
    alleles = [['C', 'A'], ['V', 'R', 'W']] # a resistance allele is uncleavable
    intro = [[1, 0, 0.099], [1, 13, 0.001]] # sex, genotype, frequency
    # genotype 0 = cc vv, 13 = [['C', 'A'], ['R', 'R']]
    s_c = [[0, ['V', 'V'], 1.0],
           [0, ['V', 'R'], 1.0],
           [0, ['R', 'V'], 1.0],
           [0, ['R', 'R'], 1.0]] #sex, alleles, fert_cost - females homozygous sterile
    f_c = []

    file_name = "mutation_data2/lof_clvr_if1_suppression_twentyPartner_femSterile_2"

    for maternal_carryover in [0, 0.3]:
        for clvr_cost in [0, 0.05, 0.1, 0.15, 0.2]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            hf_c = [[0, ['C', 'W'], 1.0, [['V'], ['R']]],
                    [1, ['C', 'W'], 1.0, []],
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]] # haploid fitness cost

            run_label = f'mc_prob_{maternal_carryover}_FC_{clvr_cost}'
            run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=POP_MAX, mc_prob=maternal_carryover)

    return None

def lof_clvr_if1_suppression_femalesterile_onePartner_MC_3():
    """runs simulations for multiple maternal carryovers and various haploid 
    fitness costs, for mating 1 female to 1 male"""
    num_partners = 1
    alleles = [['C', 'A'], ['V', 'R', 'W']] # a resistance allele is uncleavable
    intro = [[1, 0, 0.08], [1, 13, 0.02]] # sex, genotype, frequency
    # genotype 0 = cc vv, 13 = [['C', 'A'], ['R', 'R']]
    s_c = [[0, ['V', 'V'], 1.0],
           [0, ['V', 'R'], 1.0],
           [0, ['R', 'V'], 1.0],
           [0, ['R', 'R'], 1.0]] #sex, alleles, fert_cost - females homozygous sterile
    f_c = []

    file_name = "mutation_data2/lof_clvr_if1_suppression_onePartner_femSterile_3"

    for maternal_carryover in [0, 0.3]:
        for clvr_cost in [0, 0.05, 0.1, 0.15, 0.2]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            hf_c = [[0, ['C', 'W'], 1.0, [['V'], ['R']]],
                    [1, ['C', 'W'], 1.0, []],
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]] # haploid fitness cost

            run_label = f'mc_prob_{maternal_carryover}_FC_{clvr_cost}'
            run_stochastic_sim(alleles, 2*NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=POP_MAX, mc_prob=maternal_carryover)

    return None

def lof_clvr_if1_suppression_femalesterile_fivePartner_MC_3():
    """runs simulations for multiple maternal carryovers and various haploid 
    fitness costs, for mating 1 female to 1 male"""
    num_partners = 5
    alleles = [['C', 'A'], ['V', 'R', 'W']] # a resistance allele is uncleavable
    intro = [[1, 0, 0.08], [1, 13, 0.02]] # sex, genotype, frequency
    # genotype 0 = cc vv, 13 = [['C', 'A'], ['R', 'R']]
    s_c = [[0, ['V', 'V'], 1.0],
           [0, ['V', 'R'], 1.0],
           [0, ['R', 'V'], 1.0],
           [0, ['R', 'R'], 1.0]] #sex, alleles, fert_cost - females homozygous sterile
    f_c = []

    file_name = "mutation_data2/lof_clvr_if1_suppression_fivePartner_femSterile_3"

    for maternal_carryover in [0, 0.3]:
        for clvr_cost in [0, 0.05, 0.1, 0.15, 0.2]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            hf_c = [[0, ['C', 'W'], 1.0, [['V'], ['R']]],
                    [1, ['C', 'W'], 1.0, []],
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]] # haploid fitness cost

            run_label = f'mc_prob_{maternal_carryover}_FC_{clvr_cost}'
            run_stochastic_sim(alleles, 2*NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=POP_MAX, mc_prob=maternal_carryover)

    return None

def lof_clvr_if1_suppression_femalesterile_twentyPartner_MC_3():
    """runs simulations for multiple maternal carryovers and various haploid 
    fitness costs, for mating 1 female to 1 male"""
    num_partners = 20
    alleles = [['C', 'A'], ['V', 'R', 'W']] # a resistance allele is uncleavable
    intro = [[1, 0, 0.08], [1, 13, 0.02]] # sex, genotype, frequency
    # genotype 0 = cc vv, 13 = [['C', 'A'], ['R', 'R']]
    s_c = [[0, ['V', 'V'], 1.0],
           [0, ['V', 'R'], 1.0],
           [0, ['R', 'V'], 1.0],
           [0, ['R', 'R'], 1.0]] #sex, alleles, fert_cost - females homozygous sterile
    f_c = []

    file_name = "mutation_data2/lof_clvr_if1_suppression_twentyPartner_femSterile_3"

    for maternal_carryover in [0, 0.3]:
        for clvr_cost in [0, 0.05, 0.1, 0.15, 0.2]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            hf_c = [[0, ['C', 'W'], 1.0, [['V'], ['R']]],
                    [1, ['C', 'W'], 1.0, []],
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]] # haploid fitness cost

            run_label = f'mc_prob_{maternal_carryover}_FC_{clvr_cost}'
            run_stochastic_sim(alleles, 2*NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=POP_MAX, mc_prob=maternal_carryover)

    return None


# region recombination distance

def rd1_clvr_suppression_fs_onePartner():
    """runs simulations for multiple maternal carryovers and various haploid 
    fitness costs, for mating 1 female to 1 male"""
    num_partners = 1
    alleles = [['C', 'A'], ['V', 'W'], ['R', 'X']] # cleavable allele, clvr, grna/cargo 
    intro = [[1, 0, 0.1]] # release of homozygous males at 10% population frequency
    s_c = [[0, ['V', 'V'], 1.0],
           [0, ['R', 'R'], 1.0]] #sex, alleles, fert_cost - females homozygous sterile
    
    recomb_d = [[50, 1], [50, 1]]

    ## note: here, an individual that has a V on one chromosome and an R on the second is NOT sterile

    f_c = []

    file_name = "mutation_data/rd_1_suppression_onePartner_femSterile"

    for maternal_carryover in [0, 0.3]:
        for clvr_cost in [0, 0.05, 0.1, 0.15]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            hf_c = [[0, ['C', 'W', 'X'], 1.0, [['R']]],
                    [1, ['C', 'W', 'X'], 1.0, []],
                    [0, ['C', 'V', 'X'], 1.0, [['R']]],
                    [1, ['C', 'V', 'X'], 1.0, []], # maternal resuce happens via 'R'
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]] # haploid fitness cost

            run_label = f'mc_prob_{maternal_carryover}_FC_{clvr_cost}'
            run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            r_d= recomb_d,
                            mut_flag= "recomb_dist_1", run_label= run_label,
                            file_name= file_name, k=POP_MAX, mc_prob=maternal_carryover)

    return None

def rd5_clvr_suppression_fs_onePartner():
    """runs simulations for multiple maternal carryovers and various haploid 
    fitness costs, for mating 1 female to 1 male"""
    num_partners = 1
    alleles = [['C', 'A'], ['V', 'W'], ['R', 'X']] # cleavable allele, clvr, grna/cargo 
    intro = [[1, 0, 0.1]] # release of homozygous males at 10% population frequency
    s_c = [[0, ['V', 'V'], 1.0],
           [0, ['R', 'R'], 1.0]] #sex, alleles, fert_cost - females homozygous sterile
    
    recomb_d = [[50, 5], [50, 5]]

    ## note: here, an individual that has a V on one chromosome and an R on the second is NOT sterile

    f_c = []

    file_name = "mutation_data/rd_5_suppression_onePartner_femSterile"

    for maternal_carryover in [0, 0.3]:
        for clvr_cost in [0, 0.05, 0.1, 0.15]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            hf_c = [[0, ['C', 'W', 'X'], 1.0, [['R']]],
                    [1, ['C', 'W', 'X'], 1.0, []],
                    [0, ['C', 'V', 'X'], 1.0, [['R']]],
                    [1, ['C', 'V', 'X'], 1.0, []], # maternal resuce happens via 'R'
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]] # haploid fitness cost

            run_label = f'mc_prob_{maternal_carryover}_FC_{clvr_cost}'
            run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            r_d= recomb_d,
                            mut_flag= "recomb_dist_1", run_label= run_label,
                            file_name= file_name, k=POP_MAX, mc_prob=maternal_carryover)

    return None

def rd1_clvr_suppression_fs_fivePartner():
    """runs simulations for multiple maternal carryovers and various haploid 
    fitness costs, for mating 1 female to 1 male"""
    num_partners = 5
    alleles = [['C', 'A'], ['V', 'W'], ['R', 'X']] # cleavable allele, clvr, grna/cargo 
    intro = [[1, 0, 0.1]] # release of homozygous males at 10% population frequency
    s_c = [[0, ['V', 'V'], 1.0],
           [0, ['R', 'R'], 1.0]] #sex, alleles, fert_cost - females homozygous sterile
    
    recomb_d = [[50, 1], [50, 1]]

    ## note: here, an individual that has a V on one chromosome and an R on the second is NOT sterile

    f_c = []

    file_name = "mutation_data/rd_1_suppression_fivePartner_femSterile"

    for maternal_carryover in [0, 0.3]:
        for clvr_cost in [0, 0.05, 0.1, 0.15]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            hf_c = [[0, ['C', 'W', 'X'], 1.0, [['R']]],
                    [1, ['C', 'W', 'X'], 1.0, []],
                    [0, ['C', 'V', 'X'], 1.0, [['R']]],
                    [1, ['C', 'V', 'X'], 1.0, []], # maternal resuce happens via 'R'
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]] # haploid fitness cost

            run_label = f'mc_prob_{maternal_carryover}_FC_{clvr_cost}'
            run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            r_d= recomb_d,
                            mut_flag= "recomb_dist_1", run_label= run_label,
                            file_name= file_name, k=POP_MAX, mc_prob=maternal_carryover)

    return None

def rd1_clvr_suppression_fs_twentyPartner():
    """runs simulations for multiple maternal carryovers and various haploid 
    fitness costs, for mating 1 female to 1 male"""
    num_partners = 20
    alleles = [['C', 'A'], ['V', 'W'], ['R', 'X']] # cleavable allele, clvr, grna/cargo 
    intro = [[1, 0, 0.1]] # release of homozygous males at 10% population frequency
    s_c = [[0, ['V', 'V'], 1.0],
           [0, ['R', 'R'], 1.0]] #sex, alleles, fert_cost - females homozygous sterile
    
    recomb_d = [[50, 1], [50, 1]]

    ## note: here, an individual that has a V on one chromosome and an R on the second is NOT sterile

    f_c = []

    file_name = "mutation_data/rd_1_suppression_twentyPartner_femSterile"

    for maternal_carryover in [0, 0.3]:
        for clvr_cost in [0, 0.05, 0.1, 0.15]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            hf_c = [[0, ['C', 'W', 'X'], 1.0, [['R']]],
                    [1, ['C', 'W', 'X'], 1.0, []],
                    [0, ['C', 'V', 'X'], 1.0, [['R']]],
                    [1, ['C', 'V', 'X'], 1.0, []], # maternal resuce happens via 'R'
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]] # haploid fitness cost

            run_label = f'mc_prob_{maternal_carryover}_FC_{clvr_cost}'
            run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            r_d= recomb_d,
                            mut_flag= "recomb_dist_1", run_label= run_label,
                            file_name= file_name, k=POP_MAX, mc_prob=maternal_carryover)

    return None

def rd1_clvr_suppression_fs_onePartner_2():
    """runs simulations for multiple maternal carryovers and various haploid 
    fitness costs, for mating 1 female to 1 male"""
    num_partners = 1
    alleles = [['C', 'A'], ['V', 'W'], ['R', 'X']] # cleavable allele, clvr, grna/cargo 
    intro = [[1, 0, 0.1]] # release of homozygous males at 10% population frequency
    s_c = [[0, ['R', 'R'], 1.0]] #sex, alleles, fert_cost - females homozygous sterile
    
    recomb_d = [[50, 1], [50, 1]]

    ## note: here, an individual that has a V on one chromosome and an R on the second is NOT sterile

    f_c = []

    file_name = "mutation_data/rd_1_suppression_onePartner_femSterile_worseFertile2"

    for maternal_carryover in [0, 0.3]:
        for clvr_cost in [0, 0.05, 0.1, 0.15]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            hf_c = [[0, ['C', 'W', 'X'], 1.0, [['R']]],
                    [1, ['C', 'W', 'X'], 1.0, []],
                    [0, ['C', 'V', 'X'], 1.0, [['R']]],
                    [1, ['C', 'V', 'X'], 1.0, []], # maternal resuce happens via 'R'
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]] # haploid fitness cost

            run_label = f'mc_prob_{maternal_carryover}_FC_{clvr_cost}'
            run_stochastic_sim(alleles, 2*NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            r_d= recomb_d,
                            mut_flag= "recomb_dist_1", run_label= run_label,
                            file_name= file_name, k=POP_MAX, mc_prob=maternal_carryover)

    return None

def rd1_clvr_suppression_fs_fivePartner_2():
    """runs simulations for multiple maternal carryovers and various haploid 
    fitness costs, for mating 1 female to 1 male"""
    num_partners = 5
    alleles = [['C', 'A'], ['V', 'W'], ['R', 'X']] # cleavable allele, clvr, grna/cargo 
    intro = [[1, 0, 0.1]] # release of homozygous males at 10% population frequency
    s_c = [[0, ['R', 'R'], 1.0]] #sex, alleles, fert_cost - females homozygous sterile
    
    recomb_d = [[50, 1], [50, 1]]

    ## note: here, an individual that has a V on one chromosome and an R on the second is NOT sterile

    f_c = []

    file_name = "mutation_data/rd_1_suppression_fivePartner_femSterile_worseFertile2"

    for maternal_carryover in [0, 0.3]:
        for clvr_cost in [0, 0.05, 0.1, 0.15]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            hf_c = [[0, ['C', 'W', 'X'], 1.0, [['R']]],
                    [1, ['C', 'W', 'X'], 1.0, []],
                    [0, ['C', 'V', 'X'], 1.0, [['R']]],
                    [1, ['C', 'V', 'X'], 1.0, []], # maternal resuce happens via 'R'
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]] # haploid fitness cost

            run_label = f'mc_prob_{maternal_carryover}_FC_{clvr_cost}'
            run_stochastic_sim(alleles, 2*NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            r_d= recomb_d,
                            mut_flag= "recomb_dist_1", run_label= run_label,
                            file_name= file_name, k=POP_MAX, mc_prob=maternal_carryover)

    return None

def rd1_clvr_suppression_fs_twentyPartner_2():
    """runs simulations for multiple maternal carryovers and various haploid 
    fitness costs, for mating 1 female to 1 male"""
    num_partners = 20
    alleles = [['C', 'A'], ['V', 'W'], ['R', 'X']] # cleavable allele, clvr, grna/cargo 
    intro = [[1, 0, 0.1]] # release of homozygous males at 10% population frequency
    s_c = [[0, ['R', 'R'], 1.0]] #sex, alleles, fert_cost - females homozygous sterile
    
    recomb_d = [[50, 1], [50, 1]]

    ## note: here, an individual that has a V on one chromosome and an R on the second is NOT sterile

    f_c = []

    file_name = "mutation_data/rd_1_suppression_twentyPartner_femSterile_worseFertile2"

    for maternal_carryover in [0, 0.3]:
        for clvr_cost in [0, 0.05, 0.1, 0.15]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            hf_c = [[0, ['C', 'W', 'X'], 1.0, [['R']]],
                    [1, ['C', 'W', 'X'], 1.0, []],
                    [0, ['C', 'V', 'X'], 1.0, [['R']]],
                    [1, ['C', 'V', 'X'], 1.0, []], # maternal resuce happens via 'R'
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]] # haploid fitness cost

            run_label = f'mc_prob_{maternal_carryover}_FC_{clvr_cost}'
            run_stochastic_sim(alleles, 2*NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            r_d= recomb_d,
                            mut_flag= "recomb_dist_1", run_label= run_label,
                            file_name= file_name, k=POP_MAX, mc_prob=maternal_carryover)

    return None


#endregion recombination distance

# endregion female sterility

##########################################
############# male sterility #############
##########################################

def resistant_suppression_malesterile_onePartner_MC():
    """runs simulations for multiple maternal carryovers and various haploid 
    fitness costs, for mating 1 female to 1 male"""
    num_partners = 1
    alleles = [['C', 'R', 'A'], ['V', 'W']] # a resistance allele is uncleavable
    intro = [[0, 0, 0.1], [0, 23, 0.025], [1, 23, 0.025]] # sex, genotype, frequency
    # genotype 0 = cc vv, genotype 23 = ra, ww (wt resistant)
    s_c = [[1, ['V', 'V'], 1.0]] #sex, alleles, fert_cost - males homozygous sterile
    f_c = []

    file_name = "mutation_data/resistant_suppression_onePartner_maleSterile"

    for maternal_carryover in [0, 0.1, 0.2, 0.3]:
        for clvr_cost in [0, 0.05, 0.1, 0.15, 0.2]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            hf_c = [[0, ['C', 'W'], 1.0, [['V']]],
                    [1, ['C', 'W'], 1.0, []],
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]] # haploid fitness cost

            run_label = f'mc_prob_{maternal_carryover}_FC_{clvr_cost}'
            run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=POP_MAX, mc_prob=maternal_carryover)

    return None

def lof_clvr_suppression_malesterile_onePartner_MC():
    """runs simulations for multiple maternal carryovers and various haploid 
    fitness costs, for mating 1 female to 1 male"""
    num_partners = 1
    alleles = [['C', 'A'], ['V', 'R', 'W']] # a resistance allele is uncleavable
    intro = [[0, 0, 0.08], [0, 13, 0.02]] # sex, genotype, frequency
    # NEED TO BE INTRODUCING FEMALES!! DUMMY
    # genotype 0 = cc vv, 13 = [['C', 'A'], ['R', 'R']]
    s_c = [[1, ['V', 'V'], 1.0],
           [1, ['V', 'R'], 1.0],
           [1, ['R', 'V'], 1.0],
           [1, ['R', 'R'], 1.0]] #sex, alleles, fert_cost - males homozygous sterile
    f_c = []

    file_name = "mutation_data/lof_clvr_suppression_onePartner_maleSterile"

    for maternal_carryover in [0, 0.1, 0.2, 0.3]:
        for clvr_cost in [0, 0.05, 0.1, 0.15, 0.2]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            hf_c = [[0, ['C', 'W'], 1.0, [['V']]],
                    [1, ['C', 'W'], 1.0, []],
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]] # haploid fitness cost

            run_label = f'mc_prob_{maternal_carryover}_FC_{clvr_cost}'
            run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=POP_MAX, mc_prob=maternal_carryover)

    return None

def lof_clvr_suppression_malesterile_fivePartner_MC():
    """runs simulations for multiple maternal carryovers and various haploid 
    fitness costs, for mating 1 female to 5 males"""
    num_partners = 5
    alleles = [['C', 'A'], ['V', 'R', 'W']] # a resistance allele is uncleavable
    intro = [[0, 0, 0.08], [0, 13, 0.02]] # sex, genotype, frequency
    # genotype 0 = cc vv, 13 = [['C', 'A'], ['R', 'R']]
    s_c = [[1, ['V', 'V'], 1.0],
           [1, ['V', 'R'], 1.0],
           [1, ['R', 'V'], 1.0],
           [1, ['R', 'R'], 1.0]] #sex, alleles, fert_cost - males homozygous sterile
    f_c = []

    file_name = "mutation_data/lof_clvr_suppression_fivePartner_maleSterile"

    for maternal_carryover in [0, 0.1, 0.2, 0.3]:
        for clvr_cost in [0, 0.05, 0.1, 0.15, 0.2]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            hf_c = [[0, ['C', 'W'], 1.0, [['V']]],
                    [1, ['C', 'W'], 1.0, []],
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]] # haploid fitness cost

            run_label = f'mc_prob_{maternal_carryover}_FC_{clvr_cost}'
            run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=POP_MAX, mc_prob=maternal_carryover)

    return None

def lof_clvr_suppression_malesterile_twentyPartner_MC():
    """runs simulations for multiple maternal carryovers and various haploid 
    fitness costs, for mating 1 female to 20 males"""
    num_partners = 20
    alleles = [['C', 'A'], ['V', 'R', 'W']] # a resistance allele is uncleavable
    intro = [[0, 0, 0.08], [0, 13, 0.02]] # sex, genotype, frequency
    # genotype 0 = cc vv, 13 = [['C', 'A'], ['R', 'R']]
    s_c = [[1, ['V', 'V'], 1.0],
           [1, ['V', 'R'], 1.0],
           [1, ['R', 'V'], 1.0],
           [1, ['R', 'R'], 1.0]] #sex, alleles, fert_cost - males homozygous sterile
    f_c = []

    file_name = "mutation_data/lof_clvr_suppression_twentyPartner_maleSterile"

    for maternal_carryover in [0, 0.1, 0.2, 0.3]:
        for clvr_cost in [0, 0.05, 0.1, 0.15, 0.2]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            hf_c = [[0, ['C', 'W'], 1.0, [['V']]],
                    [1, ['C', 'W'], 1.0, []],
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]] # haploid fitness cost

            run_label = f'mc_prob_{maternal_carryover}_FC_{clvr_cost}'
            run_stochastic_sim(alleles, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=POP_MAX, mc_prob=maternal_carryover)

    return None

##########################################
############ fitness benefit #############
############ heatmap stuff!  #############
##########################################

# region fb heatmap
def fitness_benefit_onePartner_heatmap_drive():
    "runs simulation of how a drive helps the spread of a beneficial allele"
    num_partners = 1
    alleles = [['C', 'A'], ['V', 'W']] # a resistance allele is uncleavable
    # genotype 0 = cc vv, 13 = [['C', 'A'], ['R', 'R']]
    s_c = [] #sex, alleles, fert_cost - no sterility!
    f_c = []

    file_name = "mutation_data/onePartner_heatmap_fb_drive"

    for intro_freq in [0.05, 0.1, 0.15, 0.2, 0.25, 0.3]:
        intro = [[1, 0, intro_freq]] # sex, genotype, frequency
        for edit_cost in [-0.1, -0.2, -0.3, -0.4, -0.5, -0.6]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            hf_c = [[0, ['C', 'W'], 1.0, [['V']]],
                    [1, ['C', 'W'], 1.0, []],
                    [0, ['C'], edit_cost, []], # haploid fitness cost
                    [1, ['C'], edit_cost, []]] # haploid fitness cost

            run_label = f'intro_freq_{intro_freq}_FC_{edit_cost}'
            run_stochastic_sim(ALLELES, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=POP_MAX)

    return None

def fitness_benefit_onePartner_heatmap_NoDrive():
    "runs simulation of how a drive helps the spread of a beneficial allele"
    num_partners = 1
    alleles = [['C', 'A'], ['V', 'W']] # a resistance allele is uncleavable
    s_c = [] #sex, alleles, fert_cost - no sterility!
    f_c = []

    file_name = "mutation_data/onePartner_heatmap_fb_nodrive"

    for intro_freq in [0.05, 0.1, 0.15, 0.2, 0.25, 0.3]:
        intro = [[1, 3, intro_freq]] # sex, genotype, frequency
        # genotype 3 should be cc ww
        for edit_cost in [-0.1, -0.2, -0.3, -0.4, -0.5, -0.6]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            hf_c = [[0, ['C'], edit_cost, []], # haploid fitness cost
                    [1, ['C'], edit_cost, []]] # haploid fitness cost

            run_label = f'intro_freq_{intro_freq}_FC_{edit_cost}'
            run_stochastic_sim(ALLELES, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=POP_MAX)

    return None

def fitness_benefit_onePartner_heatmap_diploidNoDrive():
    "runs simulation of how a drive helps the spread of a beneficial allele"
    num_partners = 1
    alleles = [['C', 'A'], ['V', 'W']]
    s_c = [] #sex, alleles, fert_cost - no sterility!
    hf_c = []

    file_name = "mutation_data/onePartner_heatmap_fb_diploidnodrive"

    for intro_freq in [0.05, 0.1, 0.15, 0.2, 0.25, 0.3]:
        intro = [[1, 3, intro_freq]] # sex, genotype, frequency
        # genotype 3 should be cc ww
        for edit_cost in [-0.1, -0.2, -0.3, -0.4, -0.5, -0.6]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            f_c = [[0, ['C'], edit_cost, []],
                   [1, ['C'], edit_cost, []]]

            run_label = f'intro_freq_{intro_freq}_FC_{edit_cost}'
            run_stochastic_sim(ALLELES, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=POP_MAX)

    return None

# region extras
def fb_diploidNoDrive():
    "runs simulation of how a drive helps the spread of a beneficial allele"
    num_partners = 1
    alleles = [['C', 'A'], ['V', 'W']]
    s_c = [] #sex, alleles, fert_cost - no sterility!
    hf_c = []

    file_name = "mutation_data/onePartner_heatmap_fb_diploidnodrive"

    for intro_freq in [0.2]:
        intro = [[1, 3, intro_freq]] # sex, genotype, frequency
        # genotype 3 should be cc ww
        for edit_cost in [-0.2, -0.3, -0.4, -0.5, -0.6]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            f_c = [[0, ['C'], edit_cost, []],
                   [1, ['C'], edit_cost, []]]

            run_label = f'intro_freq_{intro_freq}_FC_{edit_cost}'
            run_stochastic_sim(ALLELES, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=POP_MAX)

    return None

def fb_diploidNoDrive2():
    "runs simulation of how a drive helps the spread of a beneficial allele"
    num_partners = 1
    alleles = [['C', 'A'], ['V', 'W']]
    s_c = [] #sex, alleles, fert_cost - no sterility!
    hf_c = []

    file_name = "mutation_data/onePartner_heatmap_fb_diploidnodrive"

    for intro_freq in [0.25]:
        intro = [[1, 3, intro_freq]] # sex, genotype, frequency
        # genotype 3 should be cc ww
        for edit_cost in [-0.1, -0.2, -0.3, -0.4, -0.5, -0.6]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            f_c = [[0, ['C'], edit_cost, []],
                   [1, ['C'], edit_cost, []]]

            run_label = f'intro_freq_{intro_freq}_FC_{edit_cost}'
            run_stochastic_sim(ALLELES, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=POP_MAX)

    return None

def fb_diploidNoDrive3():
    "runs simulation of how a drive helps the spread of a beneficial allele"
    num_partners = 1
    alleles = [['C', 'A'], ['V', 'W']]
    s_c = [] #sex, alleles, fert_cost - no sterility!
    hf_c = []

    file_name = "mutation_data/onePartner_heatmap_fb_diploidnodrive"

    for intro_freq in [0.3]:
        intro = [[1, 3, intro_freq]] # sex, genotype, frequency
        # genotype 3 should be cc ww
        for edit_cost in [-0.1, -0.2, -0.3, -0.4, -0.5, -0.6]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            f_c = [[0, ['C'], edit_cost, []],
                   [1, ['C'], edit_cost, []]]

            run_label = f'intro_freq_{intro_freq}_FC_{edit_cost}'
            run_stochastic_sim(ALLELES, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=POP_MAX)

    return None
# endregion extras

# region explore FC types

def fb_1p_heatmap_diploid_dominant1():
    "runs simulation of how a drive helps the spread of a beneficial allele"
    num_partners = 1
    alleles = [['C', 'A'], ['V', 'W']]
    s_c = [] #sex, alleles, fert_cost - no sterility!
    hf_c = []

    file_name = "mutation_data2/onePartner_heatmap_fb_diploidnodrive_dom"

    for intro_freq in [0.05, 0.1]:
        intro = [[1, 3, intro_freq]] # sex, genotype, frequency
        # genotype 3 should be cc ww
        for edit_cost in [-0.1, -0.2, -0.3, -0.4, -0.5, -0.6]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            f_c = [[0, ['C'], edit_cost, []],
                   [1, ['C'], edit_cost, []]]

            run_label = f'intro_freq_{intro_freq}_FC_{edit_cost}'
            run_stochastic_sim(ALLELES, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=POP_MAX)

    return None

def fb_1p_heatmap_diploid_dominant2():
    "runs simulation of how a drive helps the spread of a beneficial allele"
    num_partners = 1
    alleles = [['C', 'A'], ['V', 'W']]
    s_c = [] #sex, alleles, fert_cost - no sterility!
    hf_c = []

    file_name = "mutation_data2/onePartner_heatmap_fb_diploidnodrive_dom"

    for intro_freq in [0.15, 0.2]:
        intro = [[1, 3, intro_freq]] # sex, genotype, frequency
        # genotype 3 should be cc ww
        for edit_cost in [-0.1, -0.2, -0.3, -0.4, -0.5, -0.6]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            f_c = [[0, ['C'], edit_cost, []],
                   [1, ['C'], edit_cost, []]]

            run_label = f'intro_freq_{intro_freq}_FC_{edit_cost}'
            run_stochastic_sim(ALLELES, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=POP_MAX)

    return None

def fb_1p_heatmap_diploid_dominant3():
    "runs simulation of how a drive helps the spread of a beneficial allele"
    num_partners = 1
    alleles = [['C', 'A'], ['V', 'W']]
    s_c = [] #sex, alleles, fert_cost - no sterility!
    hf_c = []

    file_name = "mutation_data2/onePartner_heatmap_fb_diploidnodrive_dom"

    for intro_freq in [0.25, 0.3]:
        intro = [[1, 3, intro_freq]] # sex, genotype, frequency
        # genotype 3 should be cc ww
        for edit_cost in [-0.1, -0.2, -0.3, -0.4, -0.5, -0.6]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            f_c = [[0, ['C'], edit_cost, []],
                   [1, ['C'], edit_cost, []]]

            run_label = f'intro_freq_{intro_freq}_FC_{edit_cost}'
            run_stochastic_sim(ALLELES, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=POP_MAX)

    return None


def fb_1p_heatmap_diploid_recessive1():
    "runs simulation of how a drive helps the spread of a beneficial allele"
    num_partners = 1
    alleles = [['C', 'A'], ['V', 'W']]
    s_c = [] #sex, alleles, fert_cost - no sterility!
    hf_c = []

    file_name = "mutation_data2/onePartner_heatmap_fb_diploidnodrive_rec"

    for intro_freq in [0.05, 0.1]:
        intro = [[1, 3, intro_freq]] # sex, genotype, frequency
        # genotype 3 should be cc ww
        for edit_cost in [-0.1, -0.2, -0.3, -0.4, -0.5, -0.6]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            f_c = [[0, ['C', 'C'], edit_cost, []],
                   [1, ['C', 'C'], edit_cost, []]]

            run_label = f'intro_freq_{intro_freq}_FC_{edit_cost}'
            run_stochastic_sim(ALLELES, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=POP_MAX)

    return None

def fb_1p_heatmap_diploid_recessive2():
    "runs simulation of how a drive helps the spread of a beneficial allele"
    num_partners = 1
    alleles = [['C', 'A'], ['V', 'W']]
    s_c = [] #sex, alleles, fert_cost - no sterility!
    hf_c = []

    file_name = "mutation_data2/onePartner_heatmap_fb_diploidnodrive_rec"

    for intro_freq in [0.15, 0.2]:
        intro = [[1, 3, intro_freq]] # sex, genotype, frequency
        # genotype 3 should be cc ww
        for edit_cost in [-0.1, -0.2, -0.3, -0.4, -0.5, -0.6]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            f_c = [[0, ['C', 'C'], edit_cost, []],
                   [1, ['C', 'C'], edit_cost, []]]

            run_label = f'intro_freq_{intro_freq}_FC_{edit_cost}'
            run_stochastic_sim(ALLELES, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=POP_MAX)

    return None

def fb_1p_heatmap_diploid_recessive3():
    "runs simulation of how a drive helps the spread of a beneficial allele"
    num_partners = 1
    alleles = [['C', 'A'], ['V', 'W']]
    s_c = [] #sex, alleles, fert_cost - no sterility!
    hf_c = []

    file_name = "mutation_data2/onePartner_heatmap_fb_diploidnodrive_rec"

    for intro_freq in [0.25, 0.3]:
        intro = [[1, 3, intro_freq]] # sex, genotype, frequency
        # genotype 3 should be cc ww
        for edit_cost in [-0.1, -0.2, -0.3, -0.4, -0.5, -0.6]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            f_c = [[0, ['C', 'C'], edit_cost, []],
                   [1, ['C', 'C'], edit_cost, []]]

            run_label = f'intro_freq_{intro_freq}_FC_{edit_cost}'
            run_stochastic_sim(ALLELES, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=POP_MAX)

    return None


def fb_1p_heatmap_diploid_additive1():
    "runs simulation of how a drive helps the spread of a beneficial allele"
    num_partners = 1
    alleles = [['C', 'A'], ['V', 'W']]
    s_c = [] #sex, alleles, fert_cost - no sterility!
    hf_c = []

    file_name = "mutation_data2/onePartner_heatmap_fb_diploidnodrive_add"

    for intro_freq in [0.05, 0.1]:
        intro = [[1, 3, intro_freq]] # sex, genotype, frequency
        # genotype 3 should be cc ww
        for edit_cost in [-0.1, -0.2, -0.3, -0.4, -0.5, -0.6]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            f_c = [[0, ['C', 'C'], edit_cost, []],
                   [1, ['C', 'C'], edit_cost, []]]

            run_label = f'intro_freq_{intro_freq}_FC_{edit_cost}'
            run_stochastic_sim(ALLELES, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=POP_MAX, FC_flag="additive")

    return None

def fb_1p_heatmap_diploid_additive2():
    "runs simulation of how a drive helps the spread of a beneficial allele"
    num_partners = 1
    alleles = [['C', 'A'], ['V', 'W']]
    s_c = [] #sex, alleles, fert_cost - no sterility!
    hf_c = []

    file_name = "mutation_data2/onePartner_heatmap_fb_diploidnodrive_add"

    for intro_freq in [0.15, 0.2]:
        intro = [[1, 3, intro_freq]] # sex, genotype, frequency
        # genotype 3 should be cc ww
        for edit_cost in [-0.1, -0.2, -0.3, -0.4, -0.5, -0.6]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            f_c = [[0, ['C', 'C'], edit_cost, []],
                   [1, ['C', 'C'], edit_cost, []]]

            run_label = f'intro_freq_{intro_freq}_FC_{edit_cost}'
            run_stochastic_sim(ALLELES, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=POP_MAX, FC_flag="additive")

    return None

def fb_1p_heatmap_diploid_additive3():
    "runs simulation of how a drive helps the spread of a beneficial allele"
    num_partners = 1
    alleles = [['C', 'A'], ['V', 'W']]
    s_c = [] #sex, alleles, fert_cost - no sterility!
    hf_c = []

    file_name = "mutation_data2/onePartner_heatmap_fb_diploidnodrive_add"

    for intro_freq in [0.25, 0.3]:
        intro = [[1, 3, intro_freq]] # sex, genotype, frequency
        # genotype 3 should be cc ww
        for edit_cost in [-0.1, -0.2, -0.3, -0.4, -0.5, -0.6]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            f_c = [[0, ['C', 'C'], edit_cost, []],
                   [1, ['C', 'C'], edit_cost, []]]

            run_label = f'intro_freq_{intro_freq}_FC_{edit_cost}'
            run_stochastic_sim(ALLELES, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=POP_MAX, FC_flag="additive")

    return None

def fb_1p_heatmap_diploid_additive1_real():
    "runs simulation of how a drive helps the spread of a beneficial allele"
    num_partners = 1
    alleles = [['C', 'A'], ['V', 'W']]
    s_c = [] #sex, alleles, fert_cost - no sterility!
    hf_c = []

    file_name = "mutation_data2/onePartner_heatmap_fb_diploidnodrive_add2"

    for intro_freq in [0.05, 0.1]:
        intro = [[1, 3, intro_freq]] # sex, genotype, frequency
        # genotype 3 should be cc ww
        for edit_cost in [-0.1, -0.2, -0.3, -0.4, -0.5, -0.6]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            f_c = [[0, ['C'], edit_cost, []],
                   [1, ['C'], edit_cost, []]]

            run_label = f'intro_freq_{intro_freq}_FC_{edit_cost}'
            run_stochastic_sim(ALLELES, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=POP_MAX, FC_flag="additive")

    return None

def fb_1p_heatmap_diploid_additive2_real():
    "runs simulation of how a drive helps the spread of a beneficial allele"
    num_partners = 1
    alleles = [['C', 'A'], ['V', 'W']]
    s_c = [] #sex, alleles, fert_cost - no sterility!
    hf_c = []

    file_name = "mutation_data2/onePartner_heatmap_fb_diploidnodrive_add2"

    for intro_freq in [0.15, 0.2]:
        intro = [[1, 3, intro_freq]] # sex, genotype, frequency
        # genotype 3 should be cc ww
        for edit_cost in [-0.1, -0.2, -0.3, -0.4, -0.5, -0.6]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            f_c = [[0, ['C'], edit_cost, []],
                   [1, ['C'], edit_cost, []]]

            run_label = f'intro_freq_{intro_freq}_FC_{edit_cost}'
            run_stochastic_sim(ALLELES, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=POP_MAX, FC_flag="additive")

    return None

def fb_1p_heatmap_diploid_additive3_real():
    "runs simulation of how a drive helps the spread of a beneficial allele"
    num_partners = 1
    alleles = [['C', 'A'], ['V', 'W']]
    s_c = [] #sex, alleles, fert_cost - no sterility!
    hf_c = []

    file_name = "mutation_data2/onePartner_heatmap_fb_diploidnodrive_add2"

    for intro_freq in [0.25, 0.3]:
        intro = [[1, 3, intro_freq]] # sex, genotype, frequency
        # genotype 3 should be cc ww
        for edit_cost in [-0.1, -0.2, -0.3, -0.4, -0.5, -0.6]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            f_c = [[0, ['C'], edit_cost, []],
                   [1, ['C'], edit_cost, []]]

            run_label = f'intro_freq_{intro_freq}_FC_{edit_cost}'
            run_stochastic_sim(ALLELES, NUM_REPS, NUM_GENS, intro,
                            f_c, hf_c, s_c, num_partners, 
                            mut_flag= "NA", run_label= run_label,
                            file_name= file_name, k=POP_MAX, FC_flag="additive")

    return None


# endregion explore FC types

# endregion fb heatmap

def main():
    function = sys.argv[1]
    print("executing function " + function)
    exec(function + "()")

if __name__ == "__main__":
    main()
