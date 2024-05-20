
# ## Section B rebuild the tree
# - Neighbor joining tree based on distance matrix, maily refer to https://www.tenderisthebyte.com/blog/2022/08/31/neighbor-joining-trees/
# - This jupyter notebook requires a distance matrix generated from Section A
# - Other reference: https://github.com/Alirezafathian/phylogenetic_tree/blob/master/notebooks/phylogenetic_tree.ipynb
# 

import numpy as np
import pandas as pd
import itertools


# TODO: you need to convert the initial D_r numpy array into dictionary form {(node i, node j): link value, .....} 
# --> this will then affects your "net_divergence" function


with open('D_r.npy','rb') as f:
        D_r = np.load(f)
with open('D_x.npy','rb') as f:
        D_x = np.load(f)

K_leaf = [e for e in range(D_x.shape[0])]





def main():
    r = net_divergence(D_r, K_leaf)
    D_adjusted = adjusted_distance(D_r, K_leaf, r)
    pair_min = min(D_adjusted, key=D_adjusted.get)
    
    

# Below are functions used ========================================================


def net_divergence(D, K_leaf):
    # Net divergence r for a node i with the rest of all active nodes in K_leaf
    # D is the latest distance matrix (written in the numpy array matrix form)
    # K_leaf : the set of active leaf nodes
    
    # output: dict contains node name and its corresponding net divergence value.
    
    L = float(len(K_leaf)) 
    r = {}
    for i in K_leaf:
        # compute i's net divergence w.r.t. all other nodes
        r[i] = 0
        for j in K_leaf:
            if  i == j: # do not compute divergence of i with itself
                break
            else:
                r[i] += D[i][j]
            r[i] = (1/(L-2)) *r[i]
    
    return r
    


def adjusted_distance(D, K_leaf, r):
    # This function is to identify which two nodes are about to be linked towards a parent node
    # output: a dict
    pairs = list(itertools.combinations(K_leaf, 2)) # get all possible 2 nodes pairs from the active leaf list, method refers to: https://stackoverflow.com/questions/20762574/combinations-with-two-elements
    D_adjusted = {}
    for p in pairs:
        D_adjusted[p] = D[p[0]][p[1]] - (r[p[0]] + r[p[1]])
    
    return D_adjusted


def distance_from_child_to_parent(picked_pair, D, r, K_leaf):
    # this is to build 2 new branches from two identified child nodes linking to a same parent node

    # if len(middle_nodes) ==0:
    #     new_inter_node_lb = max(K_leaf) + 1
    # else:
    #     new_inter_node_lb = max(middle_nodes) + 1
    new_parent_lb = max(K_leaf)+1 # lb of new parent node
    K_leaf.append(new_parent_lb)
    K_leaf.remove(picked_pair[0])
    K_leaf.remove(picked_pair[1])
    
        
    # update the distance matrix by removing the picked pair (2 leaf nodes) and adding one parent node (treated as a new leaf node, although it is the intermediate node when we look at the whole tree graph)
    D_updated = {}
    '''
    - Compute distance from child (removed pairs) to parent (newly created intermediate node)
    - For removed pair <i,j> --> use d(ik) = [d(ij) + r(i) - r(j)] / 2    &   d(jk) = [d(ij) + r(j) - r(i)] / 2
    '''    
    D_updated[(p[0], new_parent_lb)] = (D[p[0]][p[1]] + r[p[0]] - r[p[1]])/2
    D_updated[(p[1], new_parent_lb)] = (D[p[0]][p[1]] + r[p[1]] - r[p[0]])/2
    
    
    '''
    - Compute the rest of leaf (non-child) to new node (aka, the parent of node removed pair <i,j>)
    - Use d(mk) = [d(im) + d(jm) - d(ij)] / 2y
    '''
    # K_leaf is now contains newly added parent node labelled with the largest label value
    pairs_new = list(itertools.combinations(K_leaf, 2)) # get all possible 2 nodes pairs from the active leaf list, method refers to: https://stackoverflow.com/questions/20762574/combinations-with-two-elements
    for p in pairs_new:
        D_adjusted[p] = D[p[0]][p[1]] - (r[p[0]] + r[p[1]])
    
    return D_adjusted

  
    
    






if __name__ == '__main__':
    main()
