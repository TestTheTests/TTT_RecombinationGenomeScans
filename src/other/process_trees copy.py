import msprime, pyslim

import numpy as np

import matplotlib.pyplot as plt



# Load the .trees file

ts = pyslim.load("results/2388682558411.trees")


# Calculate tree heights, giving uncoalesced sites the maximum time
def tree_heights(ts):
    heights = np.zeros(ts.num_trees + 1)
    for tree in ts.trees():
        if tree.num_roots > 1: # not fully coalesced
            heights[tree.index] = ts.slim_generation
        else:
            children = tree.children(tree.root)
            real_root = tree.root if len(children) > 1 else children[0]
            heights[tree.index] = tree.time(real_root)
    heights[-1] = heights[-2] # repeat the last entry for plotting
    return heights
 





# Recapitate!

recap = ts.recapitate(recombination_rate=1e-05, Ne=1e3, random_seed=1)

# Mutation! 

mutated = msprime.mutate(recap, rate=1e-6, random_seed=1, keep=True)
mutated.dump("results/2388682558411_recap_mutate.trees")




# Plot tree hieghts after recapitation
breakpoints = list(recap.breakpoints())
heights = tree_heights(recap)
plt.step(breakpoints, heights, where='post')

# Plot tree heights before recapitation
breakpoints1 = list(ts.breakpoints())
heights1 = tree_heights(ts)
plt.step(breakpoints1, heights1, where='post')

plt.show()

#
mutated.num_sites
mutated.num_samples


af = np.array([])
for variant in mutated.variants():
    result = np.sum(variant.genotypes)/(mutated.num_samples)
    af = np.append(af, result)

keep = np.logical_and(af>0.01, af < 0.99)


# I tried writing this code to get mutations that are not fixed,
# But I'm not sure if the trees format can handle it
def strip_MAF(ts):
    tables = ts.dump_tables()
    tables.sites.clear()
    tables.mutations.clear()
    for tree in ts.trees():
        for site in tree.sites():
            assert len(site.mutations) == 1  # Only supports infinite sites muts.
            mut = site.mutations[0]
            if np.logical_and(tree.num_samples(mut.node) > 0.01*ts.num_samples, tree.num_samples(mut.node) < 0.01*ts.num_samples): #
                site_id = tables.sites.add_row(position=site.position, ancestral_state=site.ancestral_state)
                tables.mutations.add_row(site=site_id, node=mut.node, derived_state=mut.derived_state)
    return tables.tree_sequence()

mutated_noMAF = strip_MAF(mutated)

with open("results/2388682558411.vcf", "w") as vcf_file:
    mutated.write_vcf(vcf_file, 2)
