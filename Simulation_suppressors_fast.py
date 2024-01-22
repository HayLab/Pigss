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
from SuppressorMutationClasses_fast import *

import scipy.stats as stats
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

@profile
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

    index_haploid_dict = {}
    haploid_index_dict = {}
    for index, haploid in enumerate(Simulation.haplotypes):
        index_haploid_dict[index] = haploid
        haploid_index_dict[haploid] = index

    num_haplos = len(Simulation.haplotypes)


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

            pairs = np.zeros((num_haplos, num_haplos))
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
                    pollens = np.zeros(num_haplos)

                    # loop through 'fathers' to pool all pollens
                    for index in father_indexes:
                        father = grouped_pollens[index]
                        bucket = father.pop(0)
                        for pollen in bucket:
                            pollens[pollen] += 1

                    ###### apply fitness costs ##############
                    """# get random integers
                    rand_floats = rng.random(len(pollens))
                    # apply fitness costs to each pollen
                    surviving_pollens = [pollen for index, pollen in enumerate(pollens)
                                         if rand_floats[index] < Simulation.haplo_fitness[(pollen, 1)]]"""

                    surviving_pollens = np.zeros(num_haplos) # a list of len [number of types of pollens],
                                          # where each value is how many of that pollen exist in the population
                    #lower = 0.0
                    #upper = 1.0
                    options = np.array(pollens).nonzero()[0]
                    for index in options:
                        amount = pollens[index]
                        pollen = index_haploid_dict[index]
                        mu = Simulation.haplo_fitness[(pollen, 1)]
                        #sigma = (-19.99*(mu-0.5)**2+5)/100
                        #distribution = stats.truncnorm((lower - mu) / sigma,
                        #                            (upper - mu) / sigma,
                        #                            loc = mu, scale = sigma)
                        
                        # we use a normal distribution to approximate how many will survive
                        # and must truncate it between 1 and 0
                        #survival = min(rng.normal(mu, sigma), 1.0)
                        #survival = max(survival, 0)
                        
                        # multiply randomized fitness value by number of pollens
                        #surviving_pollens.append(round(distribution.rvs(1)[0] * amount))
                        surviving_pollens[index] = round(mu * amount)


                    # apply fitness costs to each ovule
                    #rand_floats = rng.random(len(ovules))
                    ovules_short = np.zeros(num_haplos)
                    for ovule in ovules:
                        ovules_short[ovule] += 1
                    # apply fitness costs to each ovule
                    #surviving_ovules = [ovule for index, ovule in enumerate(ovules)
                    #                    if rand_floats[index] < Simulation.haplo_fitness[(ovule, 0)]]
                    
                    surviving_ovules = np.zeros(num_haplos) # a list of len [number of types of pollens],
                                          # where each value is how many of that pollen exist in the population
                    
                    options = np.array(ovules_short).nonzero()[0]
                    for index in options:
                        amount = ovules_short[index]
                        ovule = index_haploid_dict[index]
                        mu = Simulation.haplo_fitness[(ovule, 1)]
                        surviving_ovules[index] = round(mu * amount)

                    # select mating pollens
                    sperm = np.zeros(num_haplos)
                    # note: there should ALWAYS be an excess of pollen
                    if sum(surviving_pollens) >= sum(surviving_ovules):
                        #sperm = list(rng.choice(surviving_pollens, len(surviving_ovules),
                        #                        replace = False))
                        options = surviving_pollens.nonzero()[0]
                        ratio = sum(surviving_ovules) / sum(surviving_pollens)
                        for index in options:
                            sperm[index] = round(surviving_pollens[index] * ratio)
                        while sum(sperm) > sum(surviving_ovules):
                            options = np.array(surviving_pollens).nonzero()[0]
                            sperm[rng.choice(options)] += -1

                        while sum(sperm) < sum(surviving_ovules):
                            options = np.array(surviving_pollens).nonzero()[0]
                            sperm[rng.choice(options)] += 1
                    else:
                        # the only case where this should not be true, is if males infertile.
                        sperm = surviving_pollens

                    """for index, single_sperm in enumerate(sperm):
                        mother_egg = surviving_ovules[index]
                        father_sperm = single_sperm
                        individual_index = Simulation.cross_dict[(mother_egg, father_sperm)]
                        temp_adults[individual_index] += 1"""

                    sperm_long = [y for index, count in enumerate(sperm) for y in [index] * int(count)]
                    ovules_long = [y for index, count in enumerate(surviving_ovules) for y in [index] * int(count)]
                    random.shuffle(ovules_long)
                    if len(sperm_long) != sum(sperm):
                        print("BIGTIME ERROR")
                    if len(ovules_long) != sum(surviving_ovules):
                        print("BIGTIME ERROR")
                    for index, single_sperm in enumerate(sperm_long):
                        # number of sperm should be equal to or less than number of ovules, so
                        # as long as ovules/eggs are shuffled, it shouldn't matter if the
                        # sperm list is shuffled or not !
                        egg_index = ovules_long[index]
                        pairs[egg_index][single_sperm] += 1

            mothers, fathers = pairs.nonzero()
            for i in range(len(mothers)):
                pair_index = mothers[i], fathers[i]
                pair_count = pairs[pair_index]
                pair = (index_haploid_dict[mothers[i]], index_haploid_dict[fathers[i]])
                individual_index = Simulation.cross_dict[pair]
                temp_adults[individual_index] += pair_count
                    
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

def RS_o01percent_femalesterile_onePartner_MC():
    """runs simulations for multiple maternal carryovers and various haploid 
    fitness costs, for mating 1 female to 1 male"""
    num_partners = 1
    alleles = [['C', 'R', 'A'], ['V', 'W']] # a resistance allele is uncleavable
    intro = [[1, 0, 0.1], [0, 23, 0.00005], [1, 23, 0.00005]] # sex, genotype, frequency
    # genotype 0 = cc vv, genotype 23 = ra, ww (wt resistant)
    s_c = [[0, ['V', 'V'], 1.0]] #sex, alleles, fert_cost - females homozygous sterile
    f_c = []

    file_name = "large_pop/small_pop_profile"
    num_reps_test = 1
    num_gens_test = 50
    #"mutation_data/RS_o01percent_onePartner_femSterile"

    K = 100000 #000 # 1 000 000

    for maternal_carryover in [0]: #, 0.3]:
        for clvr_cost in [0]: #, 0.05, 0.1, 0.15]:
            # fitness costs take the form [sex, required alleles, fitness cost, rescueAlleles]
            hf_c = [[0, ['C', 'W'], 1.0, [['V']]],
                    [1, ['C', 'W'], 1.0, []],
                    [0, ['V'], clvr_cost, []], # haploid fitness cost
                    [1, ['V'], clvr_cost, []]] # haploid fitness cost

            run_label = f'mc_prob_{maternal_carryover}_FC_{clvr_cost}'
            run_stochastic_sim(alleles, num_reps_test, num_gens_test, intro,
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

def main():
    RS_o01percent_femalesterile_onePartner_MC()
    #function = sys.argv[1]
    #print("executing function " + function)
    #exec(function + "()")

if __name__ == "__main__":
    main()
