"""
Author(s)     : Michelle Johnson and Tobin Ivy
Lab           : Hay Lab
Description   : Performing stochastic simulations on gene drive, where individual gametes
are tracked. Additionally, this module uses classes to keep track of .. certain things

NOTE: does not currently handle XY sex determination, only autosomal
"""
import copy
import numpy as np

rng = np.random.default_rng()

def product(*args, repeat=1):
    """modified code from Python Software Foundations description of itertools' product function,
    product produces the "Cartesian product of input iterables".
    
    params:
        *args   -- list of loci, each loci is an arg. example: [['C', 'A'], ['V', 'W']]
        repeat  -- unknown
    returns:
        NA      -- this function is a generator, it 'yields' instead"""

    pools = [pool for pool in args] * repeat
    result = [[]]
    for pool in pools:
        result = [x+[y] for x in result for y in pool]

    for prod in result:
        yield prod

def product_index(*args, repeat=1):
    """modified code from Python Software Foundations description of itertools' product function,
    product_index tracks the indices of a given product
    
    params:
        *args   -- list of loci, each loci is an arg. example: [['C', 'A'], ['V', 'W']]
        repeat  -- unknown
    returns:
        NA      -- this function is a generator, it 'yields' instead"""
    
    pools = [pool for pool in args] * repeat
    result = [[]]
    for pool in pools:
        result = [x+[j] for i, x in enumerate(result) for j, y in enumerate(pool)]

    for prod in result:
        yield prod

def all_option(subset, overset):
    """Determines whether all elements of the subset are in the overset.
    NOTE: if subset = [], returns true

    params:
        subset  -- list of values
        overset -- list of values
    returns:
        check   -- int, 1 if true, 0 if false"""
    
    # make copies of subset and overset, as to not edit them
    subsetcopy = copy.deepcopy(subset)
    oversetcopy = copy.deepcopy(overset)

    # preset check to true
    check = 1
    for item in subsetcopy:
        if item in oversetcopy:
            oversetcopy.remove(item)
            
        # if an item in subset is not in overset, return false
        else:
            check = 0
            break

    return(check)

class Diploid:
    """Immutable diploid object, containing genotype and alleles
    param:
        genotype    -- list of loci lists of alleles, example: [['C', 'C'], ['V', 'W']]
        alleles     -- flattened version of genotype, example: ['C', 'C', 'V', 'W']
    """
    # TODO: docstring for each function
    def __init__(self, genotype, alleles) -> None:
        """initialize instance of a diploid"""
        self._genotype: list = copy.deepcopy(genotype)
        self._alleles: list = copy.deepcopy(alleles)
    
    def __eq__(self, __value: object) -> bool:
        """diploids equal if genotype and alleles are equal"""
        return ((self._genotype == __value.genotype) &
                (self._alleles == __value.alleles))
    
    def __ne__(self, __value: object) -> bool:
        """diploids not equal if genotype or alleles not equal"""
        return ((self._genotype != __value.genotype) |
                (self._alleles != __value.alleles))
    
    def __hash__(self) -> bool:
        """diploids hashable by being condensed into strng"""
        return hash(str(self._genotype) + str(self._alleles))

    def __str__(self) -> str:
        """see string version of the haploid"""
        return str(self._genotype)

    @property
    def genotype(self):
        return self._genotype
    @property
    def alleles(self):
        return self._alleles

class Haploid:
    """ Immutable haploid object
    param:
        alleles -- list of alleles, of form ['1', '2']
        parent  -- list of alleles of the parent, of the form [['locus1-1', 'locus1-2'], ['locus2-1', 'locus2-2']]
    """

    def __init__(self, alleles, parent):
        """initialize instance of a haploid class given alleles and parent"""
        self._alleles: list = copy.deepcopy(alleles)
        self._parent: Diploid = copy.deepcopy(parent)

    def __eq__(self, other: object) -> bool:
        """haploids equal if alleles and parent are equal"""
        return ((self._alleles == other.alleles) &
                (self._parent == other.parent))
    
    def __ne__(self, other:object) -> bool:
        """haploids unequal if alleles or parent unequal"""
        return ((self._alleles != other.alleles) |
                (self._parent != other.parent))
    
    def __hash__(self) -> bool:
        """haploids hashable by being represented as a string"""
        return hash(str(self._alleles) + str(self._parent))
    
    def __str__(self) -> str:
        """see string version of the haploid"""
        return str(self._alleles) + ": " + str(self._parent)

    @property
    def alleles(self):
        return self._alleles
    
    @property
    def parent(self):
        return self._parent
    
    def cross(self, other: object) -> Diploid:
        """crosses two haploids, self and other, to generate single diploid individual.
        Used later in cross_dict"""
        mom_a, mom_b = self._alleles
        dad_a, dad_b = other.alleles
        cross_genotype = [[mom_a, dad_a], [mom_b, dad_b]]
        cross_alleles = [mom_a, dad_a, mom_b, dad_b]
        return Diploid(cross_genotype, cross_alleles)


class StochasticSim:
    """Class used to represent a single stochastic simulation.
    initial params:
        num_gens        -- int number of generations for each
        alleles         -- list of loci lists of alleles in the simulation, of the form
                            [[locus1-1, locus1-2], [locus2-1, locus2-2]] or [['C', 'W'], ['V', 'A']]
        intro           -- introduction information, of the form [[int:sex, list:genotype, float:frequency]]
        fitness_costs   -- list of fitness costs associated with diploids. Each takes the form
                            [sex, list of alleles required for fitness cost, fitness cost, list of alleles that rescue]
        haplo_fitness_costs -- list of fitness costs associated with haploids. Each takes the form
                                [sex, list of alleles required for fitness cost, fitness cost, list of PARENT alleles that rescue],
                                parent alleles rescuing via maternal carryover
        mc_prob         -- float probability that a haploid is rescued by maternal carryover
        sterility_costs -- list of sterility costs, each sterility cost takes the form
                            [sex, list of alleles required, sterility cost]
        cross_dict      -- dict storing index in genotypes of the offspring produced for each mother/father haploid pair
        gametogenesis_dict  -- dict storing probabilities of each haploid gamete being produced given a specific parent
        recomb_distances    -- list of int recombination distances, values range from 0 (co-inherited) to 50 (separate chromosomes)
        add_intro       -- list of additional intros, same form as intro + a release generation
        cleave_efficiency   -- float probability that cleavage occurs
        k               -- int carrying capacity
        n_ovules        -- int number of ovules per plant, similar to num_offspring
        growth_factor   -- int growth factor, represents low-density multiplier for population growth
        num_partners    -- number of males each female plant mates with
        mutation_flag   -- mutation type flag, tells what kind of mutation you want to simulate

    created params / attributes:
        num_pollen      -- number of pollen each male produces
        all_alleles     -- flattened list version of alleles
        genotypes       -- list of all diploids possible, given the possible alleles
        haplotypes      -- list of all haploids possible, given the possible alleles & genotypes
        fitness         -- list of lists (0 = female, 1 = male) mapping genotype[index] to a fitness
        haplo_fitness   -- dictionary mapping haploid & sex to their fitness (0 = always dies, 1 = always survives)
        fertility       -- dictionary mapping genotypes to fertility (modifier on number of offspring produced)
        n_recomb_distances 
        recomb_d_index_list
        adults          -- list of lists (0 = female, 1 = male)[generation] of how many of each genotype of individual exist in a generation

    Can test using form:
    object = StochasticSim(50, [['C', 'A'], ['V', 'W']], [[1, 0, 0.2]], [], [], 0, [], {}, {},
                             [50], [[]], 0.95, 100, 30, 6, 1, "no") 
    object = StochasticSim(50, [['C', 'A'], ['V', 'W'], ['R', 'X']], [[1, 0, 0.2]], [], [], 0, [], {}, {},
                           [[50, 1], [50, 1]], [[]], 0.95, 100, 30, 6, 1, "recomb_dist_1", 'NA')
    """

    def __init__(self, num_gens, alleles, intro, fitness_costs, haplo_fitness_costs, mc_prob, sterility_costs, cross_dict, gametogenesis_dict,
                 recomb_distances, add_intro, cleave_efficiency, k, n_ovules, growth_factor, num_partners, mutation_flag, FC_flag):
        """given inputs, initialize a StochasticSim class object, that contains everything needed to perform a simulation."""
        self.num_ovules = n_ovules # approx. number of seeds in a seed pod
        self.num_pollen = 100 # approx. number of pollen/anther
        # 100 chosen as > 2 * ovules, because we always want excess pollen

        self.num_partners = num_partners
        self.add_intro = add_intro

        self.all_alleles: list = [allele for locus in alleles for allele in locus]
        self.num_gens: int = num_gens # number of Generations
        self.genotypes = []
        self.haplotypes = []
        # females take position [0] and males take position [1] in fitness and fertility
        self.fitness: list = [[], []]
        self.haplo_fitness: dict = {}
        self.fertility: dict = {}
        # must generate fitness&fertility first, because generate genotype fills them with 1's
        self.__generate_genotype(alleles)
        if FC_flag == "additive":
            self.__edit_fitnesses_additive(fitness_costs)
        else:
            self.__edit_fitnesses(fitness_costs)
        self.haplo_fitness_costs = haplo_fitness_costs
        self.mc_prob = mc_prob
        self.__edit_haplo_fitnesses()
        self.__edit_fertilities(sterility_costs)
        self.cross_dict = cross_dict #TODO: rename
        self.__update_cross_dict()
        self.gametogenesis_dict = gametogenesis_dict
        self.n_recomb_distances = len(alleles) - 1
        self.recomb_d_index_list = []
        self.__edit_recomb_distances(recomb_distances)
        self.cleave_efficiency = cleave_efficiency
        self.n_ovules: int = n_ovules
        self.growth_factor = growth_factor

        if growth_factor > (1/2) * self.num_ovules:
            raise Exception(f'Growth factor can not be greater than half the number of ovules \
                            \ninput growth factor: {growth_factor} \
                            \nnumber of ovules: {self.num_ovules}')
        self.k: int = k
        self.mutation_flag = mutation_flag

        self.test1 = 0

        if self.mutation_flag == "recomb_dist_1":
            self.produce_ovules = self.produce_ovules_rd
            self.produce_pollen = self.produce_pollen_rd
            self.test1 = 1

        self.additonal_release_list = []
        if add_intro[0] != []:
            for add_release in add_intro:
                self.additonal_release_list.append([add_release[3] + add_release[4]*g for g in range(add_release[5])])
        
        self.__initialize_adults(intro)

    def __generate_genotype(self, alleles):
        """Generates genotypes and haplotypes list from input alleles"""
        # get diploid versions of each loci
        diploid_loci = [list(product(allele, allele)) for allele in alleles]
        # turn list of diploid loci into all possible offspring genotypes
        genotypes = list(product(*diploid_loci))

        # initialize SET to store haploids in
        temp_haplotypes = set()
        # iterate through all possible individual genotypes
        for genes in genotypes:
            # get alleles, and two diploid objects
            alleles = [a for l in genes for a in l]
            individual = Diploid(genes, alleles)
            self.genotypes.append(individual)

            # fill fertility with all 1's
            self.fertility[(individual, 0)] = 1.0
            self.fertility[(individual, 1)] = 1.0

            # generate haplotypes
            # add modification if necessary
            copy_genes = copy.deepcopy(genes)
            if ('V' in alleles) & ('C' not in alleles):
                copy_genes[0].append('C')
            # get all possible haplotypes
            haplotypes = list(product(*copy_genes))
            # initialize list to store them in 
            temp_gene_haplotypes = []
            # create a haploid for each haplotype
            for haplotype in haplotypes:
                child = Haploid(haplotype, alleles)
                temp_gene_haplotypes.append(child)
            # add haploids to set: this will remove duplicates
            temp_haplotypes.update(temp_gene_haplotypes)
        
        # add all haplotypes to official storage
        self.haplotypes = list(temp_haplotypes)

    def __edit_fitnesses(self, f_c):
        """creates dictionary mapping genotypes to fitness costs IN DIPLOIDS
            f_c --list of fitness costs [sex, required alleles, fitness cost]"""
        
        # make lists appropriate lengths
        self.fitness = [[1.0]*len(self.genotypes), [1.0]*len(self.genotypes)]
        # iterate through all haploid fitness costs
        for sex, alleles, fit_cost, rescue_alleles in f_c:
            # iterate through all genotypes
            for index, diploid in enumerate(self.genotypes):
                # if ahplotype alleles match
                if all_option(alleles, diploid.alleles):
                    # no rescue option:
                    self.fitness[sex][index] *= (1 - fit_cost)

    def __edit_fitnesses_additive(self, f_c):
        """edit_fitnesses, but for recessive fitness costs"""
        # make lists appropriate lengths
        self.fitness = [[1.0]*len(self.genotypes), [1.0]*len(self.genotypes)]
        # iterate through all haploid fitness costs
        for sex, alleles, fit_cost, rescue_alleles in f_c:
            # iterate through all genotypes
            for index, diploid in enumerate(self.genotypes):
                # if ahplotype alleles match
                if all_option(alleles, diploid.alleles):
                    # no rescue option:
                    self.fitness[sex][index] *= (1 - fit_cost*diploid.alleles.count(alleles[0]))



    def __edit_haplo_fitnesses(self):
        """creates list mapping genotypes to fitness costs IN HAPLOIDS
            hf_c --list of fitness costs [sex, required alleles, fitness cost]
            
            haplo_fitness[(haploid, sex)] = fitness cost"""
        # initialize dictionary
        for haploid in self.haplotypes:
            self.haplo_fitness[(haploid, 0)] = 1.0
            self.haplo_fitness[(haploid, 1)] = 1.0

        # iterate through all haploid fitness costs
        for sex, alleles, fit_cost, rescue_alleles in self.haplo_fitness_costs:
            # iterate through all haplotypes
            for haploid in self.haplotypes:
                # if ahplotype alleles match
                if all_option(alleles, haploid.alleles):
                    # check rescue
                    if (rescue_alleles != []):
                        not_rescued = True
                        # check all possible 
                        for required_set in rescue_alleles:
                            if all_option(required_set, haploid.parent):
                                # rescue present! dependent on maternal carryover
                                self.haplo_fitness[(haploid, sex)] *= self.mc_prob
                                not_rescued = False
                        if not_rescued:
                            self.haplo_fitness[(haploid, sex)] *= (1 - fit_cost)
                                
                    else:
                        self.haplo_fitness[(haploid, sex)] *= (1 - fit_cost)

    def __edit_fertilities(self, s_c):
        """creates dictionary mapping genotypes to fertility
            s_c --list of fertility costs [sex, required alleles, fertility cost]"""
        # fertility is list containing [female fert costs], [male fert costs]
        for sex, alleles, fert_cost in s_c:
            for diploid in self.genotypes:
                if all_option(alleles, diploid.alleles):
                    self.fertility[(diploid, sex)] *= (1 - fert_cost)

    def __update_cross_dict(self):
        """make a cross dictionary, for each haploid mother/father cross, what is the 
        diploid child?"""
        if (self.cross_dict == {}):
            for mother in self.haplotypes:
                for father in self.haplotypes:
                    mom_alleles = mother.alleles
                    dad_alleles = father.alleles
                    cross_genotype = [[mom_alleles[i], dad_alleles[i]] for i in range(len(mom_alleles))]
                    cross_alleles = [allele for locus in cross_genotype for allele in locus]
                    diploid = Diploid(cross_genotype, cross_alleles)
                    self.cross_dict[(mother,father)] = self.genotypes.index(diploid)

    def __edit_recomb_distances(self, r_d):
        """ Generate recombination distances list, given recombination distances list.
        This is for HAPLOID recombinations"""
        # Set recombination rates based on recombination distances, for HAPLOID !
        # r_d should take the form [recombination distance between L1 & L2,
        #                           recombination distance between L2 & L3, etc.]

        # begin setting up recombination distances list
        r_d_full = []

        if self.n_recomb_distances > 0: #if we have multiple loci
            for rec in r_d: # for each possible recombination
                r_d_temp = []
                [r_d_temp.append([0.5 + 0.5*(1-r/50), 0.5*(r/50)]) for r in rec]
                r_d_full.append(r_d_temp)

            for rec in r_d_full[0:2]:
                r_d_full.append(list(product(*rec)))

            for rec in r_d_full[2:4]:
                r_d_full.append([np.prod(rd) for rd in rec])

            for rec in r_d_full[4:6]:
                r_d_full.append([rd / 2 for rd in rec])
            # sum of all possibilities must be 1/2, not 1
            # because there are two ways to do each
                # you can co-inherit by giving both from your mom
                # OR by giving both from your dad
            
        else:
            r_d_full = [[1]]*8

        self.recomb_d_index_list = r_d_full[7]

    def __initialize_adults(self, intro):
        """list of lists. Each list is a sex, first female second is male.
         Each of those lists has (x number of generations) of lists.
         This smallest list is an array of values, where index = index of associated
         genotype"""
        self.adults = [[[]], [[]]]
        # fill adults with empty values, 0 of each type of individuals
        self.adults[0][0] = [0] * len(self.genotypes)
        self.adults[1][0] = [0] * len(self.genotypes)
        # add in wildtype individuals, assume last one is wildtype
        self.adults[0][0][-1] += int(self.k/2)
        self.adults[1][0][-1] += int(self.k/2)

        # add introduction
        for sub_intro in intro:
            # unpack introduction
            sex = sub_intro[0]
            gene_index = sub_intro[1]
            percent = sub_intro[2]
            # add that amount
            self.adults[sex][0][gene_index] += int(percent * self.k)

    def get_gametes(self, mother_counts, father_counts) -> list:  
        """given a list of counts for mothers and fathers, with indicies corresponding
        to the genotypes in self.genotypes, generates ovules and pollens
        
        returns
            ovules  -- list of lists, where each list represents a mother and contains
                        the set of ovules that single mother has
            pollens -- list of lists of lists, where each list represents a father and contains
                        the set of pollens, broken up into Simulation.num_partners buckets"""

        # initialize empty ovule and pollen vars
        ovules = []
        pollens = []

        # get ovules
        for index, mother_count in enumerate(mother_counts):
            mother = self.genotypes[index]
            new_ovules = self.produce_ovules(mother, mother_count)
            ovules.extend(new_ovules)

        for index, father_count in enumerate(father_counts):
            father = self.genotypes[index]
            new_pollens = self.produce_pollen(father, father_count)
            pollens.extend(new_pollens)

        return ovules, pollens
        
    def produce_ovules(self, mother: Diploid, mother_count: int) -> list:
        """produces list of mother_count number of sublists, where each sublist is a set
        of ovules
        param:
            mother  -- a single diploid individual, acting as a mother
            mother_count    -- how many of those mothers exist in the population"""
        # if mother has not been seen before, we must solve for gamete chances
        if (mother, 0) not in self.gametogenesis_dict.keys():
            # generate list of possible gametes
            possible_gametes = list(product(*mother.genotype))
            # each gamete should have equal chance, by mendelian genetics
            chance = 1 / len(possible_gametes)
            # track all allele possibilities and their probability in a list
            possible_gametes_list = [[alleles, chance] for alleles in possible_gametes]

            # note: these modifications are for the specified arabidopsis
            # loop through possible gametes & unpack
            for index, [alleles, prob] in enumerate(possible_gametes_list):
                # perform cleavage in mother
                if ('V' in mother.alleles) & ('A' in alleles):
                    # get modified offspring
                    new_alleles = ['C', alleles[1]]
                    possible_gametes_list.append([new_alleles, prob * self.cleave_efficiency])
                    # edit original offspring
                    possible_gametes_list[index] = [alleles, prob * (1-self.cleave_efficiency)]
            
            # loop through every possible haplotype, giving chances for production
            all_chances = np.zeros(len(self.haplotypes))
            for alleles, prob in possible_gametes_list:
                haploid = Haploid(alleles, mother.alleles)
                h_index = self.haplotypes.index(haploid)
                all_chances[h_index] += prob

            # add to dictionary
            self.gametogenesis_dict[(mother, 0)] = all_chances
        
        # pull probability vector from dictionary
        poss_gamete_chances = self.gametogenesis_dict[(mother, 0)]
        # get number of offspring PER MOTHER
        num_offspring = int(self.num_ovules*self.fertility[(mother,0)])
        # randomly choose offspring for ALL mothers
        ovules = rng.choice(self.haplotypes, size=num_offspring*mother_count, p = poss_gamete_chances)
        # split ovules into subsets, one grouped list for each individal
        ovules_sliced = [list(ovules[num_offspring*x:num_offspring*(x+1)]) for x in range(mother_count)]
            
        return ovules_sliced
    
    def produce_pollen(self, father: Diploid, father_count: int) -> list:
        """produces list of father_count number of sublists, where each sublist has self.num_partners
        sub-sublists, and each sub-sublist is a set of pollen
        param:
            father  -- single diploid invididual, acting as father
            father_count    -- how many of these fathers exist in the population"""
        # if father has not been seen before, we must solve for gamete chances
        if (father, 1) not in self.gametogenesis_dict.keys():
            # generate list of possible gametes
            possible_gametes = list(product(*father.genotype))
            # each gamete should have equal chance, by mendelian genetics
            chance = 1 / len(possible_gametes)
            # track all allele possibilities and their probability in a list
            possible_gametes_list = [[alleles, chance] for alleles in possible_gametes]

            # note: these modifications are for the specified arabidopsis
            # loop through possible gametes & unpack
            for index, [alleles, prob] in enumerate(possible_gametes_list):
                # perform cleavage in mother
                if ('V' in father.alleles) & ('A' in alleles):
                    # get modified offspring
                    new_alleles = ['C', alleles[1]]
                    possible_gametes_list.append([new_alleles, prob * self.cleave_efficiency])
                    # edit original offspring
                    possible_gametes_list[index] = [alleles, prob * (1-self.cleave_efficiency)]
            
            # loop through every possible haplotype, giving chances for production
            all_chances = np.zeros(len(self.haplotypes))
            for alleles, prob in possible_gametes_list:
                haploid = Haploid(alleles, father.alleles)
                # we may have generated new haploids by modification above: add those here
                h_index = self.haplotypes.index(haploid)
                all_chances[h_index] += prob

            # add to dictionary
            self.gametogenesis_dict[(father, 1)] = all_chances
        
        # pull probability vector from dictionary
        poss_gamete_chances = self.gametogenesis_dict[(father, 1)]
        # randomly choose offspring PER FATHER
        num_offspring = int(self.num_pollen*self.fertility[(father,1)])
        # randomly choose offspring for ALL FATHERS
        pollen = rng.choice(self.haplotypes, size=num_offspring*father_count, p = poss_gamete_chances)

        # remove empty sets - I think we can remove this ?
        #while [] in pollen:
        #    pollen.remove([])

        # split pollen into subsets
        pollen_sliced = [list(pollen[num_offspring*x:num_offspring*(x+1)]) for x in range(father_count)]

        # split subsets (broods) into buckets
        for index, brood in enumerate(pollen_sliced):
            bucket_len = len(brood)/self.num_partners
            split_brood = [brood[int(x*bucket_len):int((x+1)*bucket_len)]
                           for x in range(self.num_partners)]
            pollen_sliced[index] = split_brood
            
        return pollen_sliced
    
    def produce_ovules_rd(self, mother: Diploid, mother_count: int) -> list:
        """produces list of mother_count number of sublists, where each sublist is a set
        of ovules
        param:
            mother  -- a single diploid individual, acting as a mother
            mother_count    -- how many of those mothers exist in the population"""
        # if mother has not been seen before, we must solve for gamete chances

        # TODO: remove this variable after proven to work
        if self.test1 == 1:
            print("CALLED RD CORRECTLY")
            self.test1 = 0
        if (mother, 0) not in self.gametogenesis_dict.keys():
            # generate list of possible gametes
            possible_gametes = list(product(*mother.genotype))
            possible_gametes_inds = list(product_index(*mother.genotype))

            possible_gametes_list = []

            # track all allele possibilities and their probability in a list
            # probabilities based on recombination distance
            for index, gamete in enumerate(possible_gametes):
                gam_inds = possible_gametes_inds[index]
                r_d_ind = ''
                for r_event in range(self.n_recomb_distances):
                    if gam_inds[r_event] == gam_inds[r_event+1]:
                        r_d_ind += '0'
                    else:
                        r_d_ind += '1'

                chance = self.recomb_d_index_list[int(r_d_ind, 2)]
                possible_gametes_list.append([gamete, chance])


            # note: these modifications are for the specified arabidopsis
            # loop through possible gametes & unpack
            for index, [alleles, prob] in enumerate(possible_gametes_list):
                # perform cleavage in mother
                if ('V' in mother.alleles) & ('R' in mother.alleles) & ('A' in alleles):
                    # get modified offspring
                    new_alleles = ['C', alleles[1], alleles[2]]
                    possible_gametes_list.append([new_alleles, prob * self.cleave_efficiency])
                    # edit original offspring
                    possible_gametes_list[index] = [alleles, prob * (1-self.cleave_efficiency)]
            
            # loop through every possible haplotype, giving chances for production
            all_chances = np.zeros(len(self.haplotypes))
            for alleles, prob in possible_gametes_list:
                haploid = Haploid(alleles, mother.alleles)
                # we may have generated new haploids by modification above: add those here
                h_index = self.haplotypes.index(haploid)
                all_chances[h_index] += prob

            # add to dictionary
            self.gametogenesis_dict[(mother, 0)] = all_chances
        
        # pull probability vector from dictionary
        poss_gamete_chances = self.gametogenesis_dict[(mother, 0)]
        # get number of offspring PER MOTHER
        num_offspring = int(self.num_ovules*self.fertility[(mother,0)])
        # randomly choose offspring for ALL mothers
        ovules = rng.choice(self.haplotypes, size=num_offspring*mother_count, p = poss_gamete_chances)
        # split ovules into subsets, one grouped list for each individal
        ovules_sliced = [list(ovules[num_offspring*x:num_offspring*(x+1)]) for x in range(mother_count)]
            
        return ovules_sliced
    
    def produce_pollen_rd(self, father: Diploid, father_count: int) -> list:
        """produces list of father_count number of sublists, where each sublist has self.num_partners
        sub-sublists, and each sub-sublist is a set of pollen
        param:
            father  -- single diploid invididual, acting as father
            father_count    -- how many of these fathers exist in the population"""
        # if father has not been seen before, we must solve for gamete chances
        if (father, 1) not in self.gametogenesis_dict.keys():
            # generate list of possible gametes
            possible_gametes = list(product(*father.genotype))
            possible_gametes_inds = list(product_index(*father.genotype))

            possible_gametes_list = []

            # track all allele possibilities and their probability in a list
            # probabilities based on recombination distance
            for index, gamete in enumerate(possible_gametes):
                gam_inds = possible_gametes_inds[index]
                r_d_ind = ''
                for r_event in range(self.n_recomb_distances):
                    if gam_inds[r_event] == gam_inds[r_event+1]:
                        r_d_ind += '0'
                    else:
                        r_d_ind += '1'

                chance = self.recomb_d_index_list[int(r_d_ind, 2)]
                possible_gametes_list.append([gamete, chance])

            # note: these modifications are for the specified arabidopsis
            # loop through possible gametes & unpack
            for index, [alleles, prob] in enumerate(possible_gametes_list):
                # perform cleavage in mother
                if ('V' in father.alleles) & ('R' in father.alleles) & ('A' in alleles):
                    # get modified offspring
                    new_alleles = ['C', alleles[1], alleles[2]]
                    possible_gametes_list.append([new_alleles, prob * self.cleave_efficiency])
                    # edit original offspring
                    possible_gametes_list[index] = [alleles, prob * (1-self.cleave_efficiency)]
            
            # loop through every possible haplotype, giving chances for production
            all_chances = np.zeros(len(self.haplotypes))
            for alleles, prob in possible_gametes_list:
                haploid = Haploid(alleles, father.alleles)
                # we may have generated new haploids by modification above: add those here
                h_index = self.haplotypes.index(haploid)
                all_chances[h_index] += prob

            # add to dictionary
            self.gametogenesis_dict[(father, 1)] = all_chances
        
        # pull probability vector from dictionary
        poss_gamete_chances = self.gametogenesis_dict[(father, 1)]
        # randomly choose offspring PER FATHER
        num_offspring = int(self.num_pollen*self.fertility[(father,1)])
        # randomly choose offspring for ALL FATHERS
        pollen = rng.choice(self.haplotypes, size=num_offspring*father_count, p = poss_gamete_chances)

        # remove empty sets
        # while [] in pollen:
        #     pollen.remove([])

        # split pollen into subsets
        pollen_sliced = [list(pollen[num_offspring*x:num_offspring*(x+1)]) for x in range(father_count)]

        # split subsets (broods) into buckets
        for index, brood in enumerate(pollen_sliced):
            bucket_len = len(brood)/self.num_partners
            split_brood = [brood[int(x*bucket_len):int((x+1)*bucket_len)]
                           for x in range(self.num_partners)]
            pollen_sliced[index] = split_brood
            
        return pollen_sliced
    
class Fast_StochasticSim:
    """Class used to represent a single stochastic simulation.
    initial params:
        num_gens        -- int number of generations for each
        alleles         -- list of loci lists of alleles in the simulation, of the form
                            [[locus1-1, locus1-2], [locus2-1, locus2-2]] or [['C', 'W'], ['V', 'A']]
        intro           -- introduction information, of the form [[int:sex, list:genotype, float:frequency]]
        fitness_costs   -- list of fitness costs associated with diploids. Each takes the form
                            [sex, list of alleles required for fitness cost, fitness cost, list of alleles that rescue]
        haplo_fitness_costs -- list of fitness costs associated with haploids. Each takes the form
                                [sex, list of alleles required for fitness cost, fitness cost, list of PARENT alleles that rescue],
                                parent alleles rescuing via maternal carryover
        mc_prob         -- float probability that a haploid is rescued by maternal carryover
        sterility_costs -- list of sterility costs, each sterility cost takes the form
                            [sex, list of alleles required, sterility cost]
        cross_dict      -- dict storing index in genotypes of the offspring produced for each mother/father haploid pair
        gametogenesis_dict  -- dict storing probabilities of each haploid gamete being produced given a specific parent
        recomb_distances    -- list of int recombination distances, values range from 0 (co-inherited) to 50 (separate chromosomes)
        add_intro       -- list of additional intros, same form as intro + a release generation
        cleave_efficiency   -- float probability that cleavage occurs
        k               -- int carrying capacity
        n_ovules        -- int number of ovules per plant, similar to num_offspring
        growth_factor   -- int growth factor, represents low-density multiplier for population growth
        num_partners    -- number of males each female plant mates with
        mutation_flag   -- mutation type flag, tells what kind of mutation you want to simulate

    created params / attributes:
        num_pollen      -- number of pollen each male produces
        all_alleles     -- flattened list version of alleles
        genotypes       -- list of all diploids possible, given the possible alleles
        haplotypes      -- list of all haploids possible, given the possible alleles & genotypes
        fitness         -- list of lists (0 = female, 1 = male) mapping genotype[index] to a fitness
        haplo_fitness   -- dictionary mapping haploid & sex to their fitness (0 = always dies, 1 = always survives)
        fertility       -- dictionary mapping genotypes to fertility (modifier on number of offspring produced)
        n_recomb_distances 
        recomb_d_index_list
        adults          -- list of lists (0 = female, 1 = male)[generation] of how many of each genotype of individual exist in a generation

    Can test using form:
    object = StochasticSim(50, [['C', 'A'], ['V', 'W']], [[1, 0, 0.2]], [], [], 0, [], {}, {},
                             [50], [[]], 0.95, 100, 30, 6, 1, "no") 
    object = StochasticSim(50, [['C', 'A'], ['V', 'W'], ['R', 'X']], [[1, 0, 0.2]], [], [], 0, [], {}, {},
                             [50, 1], [[]], 0.95, 100, 30, 6, 1, "recomb_dist_1"))
    """

    def __init__(self, num_gens, alleles, intro, fitness_costs, haplo_fitness_costs, mc_prob, sterility_costs, cross_dict, gametogenesis_dict,
                 recomb_distances, add_intro, cleave_efficiency, k, n_ovules, growth_factor, num_partners, mutation_flag, FC_flag):
        """given inputs, initialize a StochasticSim class object, that contains everything needed to perform a simulation."""
        self.num_ovules = n_ovules # approx. number of seeds in a seed pod
        self.num_pollen = 100 # approx. number of pollen/anther
        # 100 chosen as > 2 * ovules, because we always want excess pollen

        self.num_partners = num_partners
        self.add_intro = add_intro

        self.all_alleles: list = [allele for locus in alleles for allele in locus]
        self.num_gens: int = num_gens # number of Generations
        self.genotypes = []
        self.haplotypes = []
        # females take position [0] and males take position [1] in fitness and fertility
        self.fitness:list = [np.ones(0), np.ones(0)]
        self.haplo_fitness: list = [np.ones(0), np.ones(0)]
        self.fertility: dict = {}
        # must generate fitness&fertility first, because generate genotype fills them with 1's
        self.__generate_genotype(alleles)
        if FC_flag == "additive":
            self.__edit_fitnesses_additive(fitness_costs)
        else:
            self.__edit_fitnesses(fitness_costs)
        self.haplo_fitness_costs = haplo_fitness_costs
        self.mc_prob = mc_prob
        self.__edit_haplo_fitnesses()
        self.__edit_fertilities(sterility_costs)
        self.cross_array = cross_dict #TODO: rename
        self.__update_cross_dict()
        self.gametogenesis_pdfs = []
        self.__initialize_gametogenesis_pdfs()
        self.father_pdfs = {}
        self.gametogenesis_dict = gametogenesis_dict
        self.n_recomb_distances = len(alleles) - 1
        self.recomb_d_index_list = []
        self.__edit_recomb_distances(recomb_distances)
        self.cleave_efficiency = cleave_efficiency
        self.n_ovules: int = n_ovules
        self.growth_factor = growth_factor

        if growth_factor > (1/2) * self.num_ovules:
            raise Exception(f'Growth factor can not be greater than half the number of ovules \
                            \ninput growth factor: {growth_factor} \
                            \nnumber of ovules: {self.num_ovules}')
        self.k: int = k
        self.mutation_flag = mutation_flag

        self.test1 = 0

        if self.mutation_flag == "recomb_dist_1":
            self.produce_ovules = self.produce_ovules_rd
            self.produce_pollen = self.produce_pollen_rd
            self.test1 = 1

        self.additonal_release_list = []
        if add_intro[0] != []:
            for add_release in add_intro:
                self.additonal_release_list.append([add_release[3] + add_release[4]*g for g in range(add_release[5])])
        
        self.__initialize_adults(intro)

    def __generate_genotype(self, alleles):
        """Generates genotypes and haplotypes list from input alleles"""
        # get diploid versions of each loci
        diploid_loci = [list(product(allele, allele)) for allele in alleles]
        # turn list of diploid loci into all possible offspring genotypes
        genotypes = list(product(*diploid_loci))

        # initialize SET to store haploids in
        temp_haplotypes = set()
        # iterate through all possible individual genotypes
        for index, genes in enumerate(genotypes):
            # get alleles, and two diploid objects
            alleles = [a for l in genes for a in l]
            individual = Diploid(genes, alleles)
            self.genotypes.append(individual)

            # fill fertility with all 1's
            self.fertility[(index, 0)] = 1.0
            self.fertility[(index, 1)] = 1.0

            # generate haplotypes
            # add modification if necessary
            copy_genes = copy.deepcopy(genes)
            if (('V' in alleles) | ('V1' in alleles)) & ('C' not in alleles):
                copy_genes[0].append('C')
            # get all possible haplotypes
            haplotypes = list(product(*copy_genes))
            # initialize list to store them in 
            temp_gene_haplotypes = []
            # create a haploid for each haplotype
            for haplotype in haplotypes:
                child = Haploid(haplotype, alleles)
                temp_gene_haplotypes.append(child)
            # add haploids to set: this will remove duplicates
            temp_haplotypes.update(temp_gene_haplotypes)
        
        # add all haplotypes to official storage
        self.haplotypes = list(temp_haplotypes)

    def __edit_fitnesses(self, f_c):
        """creates dictionary mapping genotypes to fitness costs IN DIPLOIDS
            f_c --list of fitness costs [sex, required alleles, fitness cost]"""
        
        # make lists appropriate lengths
        self.fitness = np.ones((len(self.genotypes), len(self.genotypes)))
        # iterate through all haploid fitness costs
        for sex, alleles, fit_cost, rescue_alleles in f_c:
            # iterate through all genotypes
            for index, diploid in enumerate(self.genotypes):
                # if ahplotype alleles match
                if all_option(alleles, diploid.alleles):
                    # no rescue option:
                    self.fitness[sex][index] *= (1 - fit_cost)

    def __edit_fitnesses_additive(self, f_c):
        """edit_fitnesses, but for recessive fitness costs"""
        # make lists appropriate lengths
        self.fitness = [np.ones(len(self.genotypes)), np.ones(len(self.genotypes))]
        # iterate through all haploid fitness costs
        for sex, alleles, fit_cost, rescue_alleles in f_c:
            # iterate through all genotypes
            for index, diploid in enumerate(self.genotypes):
                # if ahplotype alleles match
                if all_option(alleles, diploid.alleles):
                    # no rescue option:
                    self.fitness[sex][index] *= (1 - fit_cost*diploid.alleles.count(alleles[0]))



    def __edit_haplo_fitnesses(self):
        """creates list mapping genotypes to fitness costs IN HAPLOIDS
            hf_c --list of fitness costs [sex, required alleles, fitness cost]
            
            haplo_fitness[(haploid, sex)] = fitness cost"""
        # initialize arrays
        self.haplo_fitness = [np.ones(len(self.haplotypes)), 
                              np.ones(len(self.haplotypes))]
        
        self.haplo_fitness[0][-1] = 0
        self.haplo_fitness[1][-1] = 0

        # iterate through all haploid fitness costs
        for sex, alleles, fit_cost, rescue_alleles in self.haplo_fitness_costs:
            # iterate through all haplotypes
            for index, haploid in enumerate(self.haplotypes):
                # if ahplotype alleles match
                if all_option(alleles, haploid.alleles):
                    # check rescue
                    if (rescue_alleles != []):
                        not_rescued = True
                        # check all possible 
                        for required_set in rescue_alleles:
                            if all_option(required_set, haploid.parent):
                                # rescue present! dependent on maternal carryover
                                self.haplo_fitness[sex][index] *= self.mc_prob
                                not_rescued = False
                        if not_rescued:
                            self.haplo_fitness[sex][index] *= (1 - fit_cost)
                                
                    else:
                        self.haplo_fitness[sex][index] *= (1 - fit_cost)

    def __edit_fertilities(self, s_c):
        """creates dictionary mapping genotypes to fertility
            s_c --list of fertility costs [sex, required alleles, fertility cost]"""
        # fertility is list containing [female fert costs], [male fert costs]
        for sex, alleles, fert_cost in s_c:
            for index, diploid in enumerate(self.genotypes):
                if all_option(alleles, diploid.alleles):
                    self.fertility[(index, sex)] *= (1 - fert_cost)

    def __update_cross_dict(self):
        """make a cross dictionary, for each haploid mother/father cross, what is the 
        diploid child?"""
        if (self.cross_array == {}):
            self.cross_array = np.zeros((len(self.haplotypes), len(self.haplotypes)))
            for m_index, mother in enumerate(self.haplotypes):
                for f_index, father in enumerate(self.haplotypes):
                    mom_alleles = mother.alleles
                    dad_alleles = father.alleles
                    cross_genotype = [[mom_alleles[i], dad_alleles[i]] for i in range(len(mom_alleles))]
                    cross_alleles = [allele for locus in cross_genotype for allele in locus]
                    diploid = Diploid(cross_genotype, cross_alleles)
                    self.cross_array[m_index][f_index] = self.genotypes.index(diploid)

    def __edit_recomb_distances(self, r_d):
        """ Generate recombination distances list, given recombination distances list.
        This is for HAPLOID recombinations"""
        # Set recombination rates based on recombination distances, for HAPLOID !
        # r_d should take the form [recombination distance between L1 & L2,
        #                           recombination distance between L2 & L3, etc.]

        # begin setting up recombination distances list
        r_d_full = []

        if self.n_recomb_distances > 0: #if we have multiple loci
            for rec in r_d: # for each possible recombination
                r_d_temp = []
                [r_d_temp.append([0.5 + 0.5*(1-r/50), 0.5*(r/50)]) for r in rec]
                r_d_full.append(r_d_temp)

            for rec in r_d_full[0:2]:
                r_d_full.append(list(product(*rec)))

            for rec in r_d_full[2:4]:
                r_d_full.append([np.prod(rd) for rd in rec])

            for rec in r_d_full[4:6]:
                r_d_full.append([rd / 2 for rd in rec])
            # sum of all possibilities must be 1/2, not 1
            # because there are two ways to do each
                # you can co-inherit by giving both from your mom
                # OR by giving both from your dad
            
        else:
            r_d_full = [[1]]*8

        self.recomb_d_index_list = r_d_full[7]

    def __initialize_adults(self, intro):
        """list of lists. Each list is a sex, first female second is male.
         Each of those lists has (x number of generations) of lists.
         This smallest list is an array of values, where index = index of associated
         genotype"""
        self.adults = [[[]], [[]]]
        # fill adults with empty values, 0 of each type of individuals
        self.adults[0][0] = np.zeros(len(self.genotypes))
        self.adults[1][0] = np.zeros(len(self.genotypes))
        # add in wildtype individuals, assume last one is wildtype
        self.adults[0][0][-1] += int(self.k/2)
        self.adults[1][0][-1] += int(self.k/2)

        # add introduction
        for sub_intro in intro:
            # unpack introduction
            sex = sub_intro[0]
            gene_index = sub_intro[1]
            percent = sub_intro[2]
            # add that amount
            self.adults[sex][0][gene_index] += int(percent * self.k)

    def __initialize_gametogenesis_pdfs(self):
        # need to be searchable by sex, genotype
        # contains pdfs (), each of length len(haplotypes)
        self.gametogenesis_pdfs = np.ones((2, len(self.genotypes), len(self.haplotypes)))

    def get_gametes(self, sex, diploid_index):
        """ diploid is an int INDEX values"""
        
        pdf = self.gametogenesis_pdfs[sex][diploid_index]
        if int(sum(pdf)) == pdf.size:

            if (diploid_index, sex) not in self.gametogenesis_dict.keys():
                self.update_gametogenesis_dict(diploid_index, sex)

            gamete_chances = self.gametogenesis_dict[(diploid_index, sex)] 
            # should combine to 1
            haplo_fitness = self.haplo_fitness[sex]
            # maximum value is 1

            combined_fitness = np.multiply(haplo_fitness, gamete_chances)
            self.gametogenesis_pdfs[sex][diploid_index] = combined_fitness
            pdf = combined_fitness

        return pdf

            
    def update_gametogenesis_dict(self, diploid_index, sex):
        # generate list of possible gametes
        diploid = self.genotypes[diploid_index]
        possible_gametes = list(product(*diploid.genotype))
        # each gamete should have equal chance, by mendelian genetics
        chance = 1 / len(possible_gametes)
        # track all allele possibilities and their probability in a list
        possible_gametes_list = [[alleles, chance] for alleles in possible_gametes]

        # note: these modifications are for the specified arabidopsis
        # loop through possible gametes & unpack
        for index, [alleles, prob] in enumerate(possible_gametes_list):
            # perform cleavage in mother
            if (('V' in diploid.alleles) | ('V1' in diploid.alleles)) & ('A' in alleles):
                # get modified offspring
                new_alleles = ['C'] + alleles[1:]
                possible_gametes_list.append([new_alleles, prob * self.cleave_efficiency])
                # edit original offspring
                possible_gametes_list[index] = [alleles, prob * (1-self.cleave_efficiency)]
        
        # loop through every possible haplotype, giving chances for production
        all_chances = np.zeros(len(self.haplotypes))
        for alleles, prob in possible_gametes_list:
            haploid = Haploid(alleles, diploid.alleles)
            h_index = self.haplotypes.index(haploid)
            all_chances[h_index] += prob

        # add to dictionary
        self.gametogenesis_dict[(diploid_index, sex)] = all_chances
