import numpy as np
import copy
import random
import networkx
import matplotlib.pyplot as plt
from scipy.sparse.csgraph import dijkstra,connected_components
from scipy import sparse
#cyptoscape / networkx

FILE = "./ppi.txt"
LABELS = "./protein_labels.txt"
GAMMA = 0.3
DRAWTIMES = 1000

#TODO: add labels to each node, as stationary frequencies
#implement shortest paths/average shortest paths between annotated proteins
#background score: choose random set of nodes and
#compare average shortest path between them
#create 1000 random sets of proteins and measure shortest path metric
#plot distribution / calculate p-val
# 
#now do RWR starting from one annotated protein to 
#calc stationary freqs for all other cells (take average)
#create background rate by creating 1000 random sets of proteins and calc each set's prox
#using RWR. Describe background distribution and p-value for RWR approach
#degree preserving edge-swap: 

#All of the above should be done, all that's left is running it and writing the report

#float truncation for use in truncating labels
def truncate(f, n):
    '''Truncates/pads a float f to n decimal places without rounding'''
    s = '{}'.format(f)
    if 'e' in s or 'E' in s:
        return '{0:.{1}f}'.format(f, n)
    i, p, d = s.partition('.')
    return '.'.join([i, (d+'0'*n)[:n]])

#calc change in scores to know when to stop
def calc_deltas(scoring, bef):
    for i in range(len(scoring)):
        delt = abs(scoring[i] - bef[i])
        if delt > 0.001: #how big should our delta be?
            return False
    return True

# Random walks with restarts, takes a graph, starting node, and gamma constant
# not for use on disjoint graphs (?)
def rand_walk_restarts(graph, node, gamma):
    #have some gamma value between 0 and 1 for restart chance,
    # want some value for where we go next
    curr = node
    # each 'step' in our case is simply an iteration in our while loop
    # we go until we break bc of deltas of nodes stabilizing
    count = [0] * len(graph)
    score = [0] * len(graph)
    before = [0] * len(graph)
    iterv = 0

    total_travel = 0
    while True:
        #print(curr)
        prob_sum = 0
        iterv += 1 #interval for when we checking scores
        neighbors = [0] * len(graph)
        total_travel += 1
        count[curr] += 1
        before[curr] = score[curr]
        score[curr] = count[curr] / total_travel
        #print(score)

        #function here to handle how we stop (calculate delta values in all of score)
        # and stop when they stabilize past a threshold
        
        if iterv % 10:
            if calc_deltas(score, before) or iterv >= 100000:
                break
        
        if random.random() <= gamma:
            curr = node
            #print("restart")    

        # check what our neighbors are and their scores
        for i in range(len(graph)):
            prev = prob_sum
            if graph[curr, i]:
                prob_sum += graph[curr, i]
                neighbors[i] = ( i, prev, prob_sum)

        choice_taken = random.randint(0 , prob_sum)
        #print("choice" , choice_taken , curr)

        for i in range(len(neighbors)):
            if (neighbors[i] and choice_taken <= neighbors[i][2] 
                and choice_taken > neighbors[i][1]):
                curr = neighbors[i][0]
                break
                #print(curr , neighbors[i])
    return score

# creates an adjacency matrix from a file representing an edgelist
# returns the generated adjacency matrix and a dictionary mapping protein names to
# indices in the matrix
def generate_matrix(file):
    dic = {}
    count = 0
    labels = []
    with open(file, 'r') as data:
        edges = data.readlines()
        for i in edges:
            s = i.split(' ')
            s[1] = s[1].strip()
            print(s)
            if s[0] not in dic:
                dic.update({s[0] : count})
                count += 1
            if s[1] not in dic:
                dic.update({s[1] : count})
                count += 1
        
    with open(LABELS, 'r') as labe:
        lab = labe.readlines()
        for i in lab:
            s = i.split('\n')
            labels.append(s[0])

    matrix = np.zeros((count, count),dtype=int)
    for i in edges:
        l = i.split(' ')
        l[1] = l[1].strip()
        print(l)
        first = dic.get(l[0])
        sec = dic.get(l[1])
        print((first, sec))
        matrix[first][sec] = 1
        matrix[sec][first] = 1
    
    return (matrix , dic, labels)

# compare the ranking of nodes based on degree / score 
# takes an adj matrix and a score list as returned by rwr       
def ranking(mat, scr):
    lst = [(0,0)] * len(mat)
    res1 = []
    res2 = []
    for i in range(len(mat)):
        for j in range(len(mat)):
            if mat[i,j]:
                lst[i] = (lst[i][0] + 1, i)
    sort_lst = sorted(lst, key=lambda x: x[0]) 
    for i in range(len(scr)):
        lst[i] = (scr[i], i)
    scr_lst = sorted(lst, key=lambda x: x[0]) 
    for i in sort_lst:
        res1.append(i[1])
    for j in scr_lst:
        res2.append(j[1])
    # list will be sorted by (degree, nodenum)
    print("degree: ", res1, "\nscr: ", res2 )

#create background sets and scorings given original matrix and size of sample set
def random_set(matr, sample_size):
    shortest_score = []
    rwr_score = []
    dist = dijkstra(matr)

    for times in range(DRAWTIMES + 1):
        sample = [0] * (sample_size + 1)
        for i in range(sample_size + 1):
            sample[i] = random.randint(0, len(matr[0]) - 1)
        #print(sample)
        short = average_shortest(dist, sample)
        shortest_score.append(short)
        scr = rand_walk_restarts(matr, sample[0], GAMMA)
        rwr_score.append(rand_walk_prox(scr, sample))
        
    return (rwr_score, shortest_score)

def rand_walk_prox(scr, sample_set):
    scoring  = 0
    for i in sample_set:
        scoring += scr[i]

    return (scoring/float(len(sample_set)))


def average_shortest(dist, sample_set):
    paths = 0
    total = 0
    for i in range(len(sample_set)):
        for j in range( i+1, len(sample_set)):
            first = sample_set[i]
            second = sample_set[j]
            paths += 1
            total += dist[first,second]
    return total/paths
    

# create the graph as shown in class as an adjacency matrix
def base_matrix():
    matrix = np.matrix([[0,1,0,0,0,0,0,0,0,0,0], 
                       [1,0,1,1,0,1,1,0,0,0,0], 
                       [0,1,0,0,1,0,0,0,0,0,0], 
                       [0,1,0,0,1,0,0,0,0,0,0], 
                       [0,0,1,1,0,0,0,0,0,0,0], 
                       [0,1,0,0,0,0,1,1,0,0,0], 
                       [0,1,0,0,0,1,0,0,0,1,0], 
                       [0,0,0,0,0,1,0,0,1,0,0], 
                       [0,0,0,0,0,0,0,1,0,0,1], 
                       [0,0,0,0,0,0,1,0,0,0,1], 
                       [0,0,0,0,0,0,0,0,1,1,0]])
    return matrix

# generate graph based on example from class and perform rwr on it
# basically part 1
def example_graphing():
    mat = base_matrix()
    scr = rand_walk_restarts(mat, 1, GAMMA)
    g = networkx.from_numpy_array(mat)
    #format universal position for example arr
    pos = [np.array([-0.4, 0.5]), np.array([-0.1, 0.3]), 
           np.array([-0.3, -0.3]), np.array([-0.5, -0.3]),
           np.array([-0.4, -0.6]), np.array([0.1, 0.3]),
           np.array([0.1, -0.1]), np.array([0.3, -0.1]),
           np.array([0.3, -0.3]), np.array([0.1, -0.3]),
           np.array([0.2, -0.7])]
    labl = {}
    for i in range(len(scr)):
        labl.update({i : truncate(scr[i], 3)})

    networkx.draw(g, pos, labels=labl, node_color=scr, with_labels=True, linewidths=15, cmap=plt.cm.YlOrRd)
    plt.show()
    ranking(mat, scr)

def main():
    res = generate_matrix(FILE)
    matr = res[0]
    cpy = copy.deepcopy(matr)
    g = networkx.from_numpy_array(cpy)

    shuffled = networkx.double_edge_swap(g, nswap=10, max_tries=100)
    name = res[1]
    labels = res[2] # labels as strings to index into names
    #example_graphing()
    dist = dijkstra(matr)

    label_indx = []
    print(labels)
    print(networkx.is_connected(g))
    for i in labels:
        label_indx.append(name[i])
    
    scr = rand_walk_restarts( matr, label_indx[0], GAMMA)
    
    actual_rwr = rand_walk_prox(scr, label_indx)
    print("completed randwalk")
    actual_shortest = average_shortest(dist, label_indx)
    print(actual_shortest, actual_rwr)
    rand_set = random_set(matr, len(label_indx))
    rwr_set = rand_set[1]
    short_set = rand_set[0]
    rwr_extreme_count = 0 #how many results of the rand rwr set are more extreme than our actual
    short_count = 0 #likewise for shortest path 

    for i in range(DRAWTIMES):
        if rwr_set[i] > actual_rwr:
            rwr_extreme_count += 1
        if short_set[i] > actual_shortest:
            short_count += 1

    pval_rwr = rwr_extreme_count / DRAWTIMES
    pval_short = short_count / DRAWTIMES
    print()
    print(f"PVALUES for rwr: {pval_rwr} and for shortest: {pval_short}")
    
    #transform into np arrays to save as csv file to turn into graph
    pr_rwr = np.asarray(rwr_set)
    pr_shortest = np.asarray(short_set)
    np.savetxt("rwr.csv", pr_rwr, delimiter=",")
    np.savetxt("short.csv", pr_shortest, delimiter=",")


    shuffled = networkx.to_numpy_array(shuffled)
    shuffled = np.asmatrix(shuffled)
    scr = rand_walk_restarts(shuffled, label_indx[0], GAMMA)
    actual_rwr = rand_walk_prox(scr, label_indx)
    dis = dijkstra(shuffled)
    actual_shortest = average_shortest(dis, label_indx)
    rand_set = random_set(shuffled, len(label_indx))
    rwr_set = rand_set[1]
    short_set = rand_set[0]

    rwr_extreme_count = 0 #how many results of the rand rwr set are more extreme than our actual
    short_count = 0 #likewise for shortest path 
    for i in range(DRAWTIMES + 1):
        if rwr_set[i] > actual_rwr:
            rwr_extreme_count += 1
        if short_set[i] > actual_shortest:
            short_count += 1
    pval_rwr = rwr_extreme_count / DRAWTIMES
    pval_short = short_count / DRAWTIMES
    print()
    print(f"UPDATED PVALUES for rwr: {pval_rwr} and for shortest: {pval_short}") 
    #intuition says that pvals should be higher

    print("==============done=============")


if __name__ == '__main__': 
    main()