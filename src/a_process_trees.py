import msprime, pyslim, argparse

import numpy as np

import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--seed", required = True)

args = vars(parser.parse_args())
seed = args['seed']

# Load the .trees file
ts = pyslim.load(seed + "_Invers.trees")

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
# mutated.dump("results/2388682558411_recap_mutate.trees")

# Plot tree hieghts after recapitation
# breakpoints = list(recap.breakpoints())
# heights = tree_heights(recap)
# plt.step(breakpoints, heights, where='post')

# Plot tree heights before recapitation
# breakpoints1 = list(ts.breakpoints())
# heights1 = tree_heights(ts)
# plt.step(breakpoints1, heights1, where='post')
# plt.show()

with open(seed + "_Invers.vcf", "w") as vcf_file:
    mutated.write_vcf(vcf_file, 2)
