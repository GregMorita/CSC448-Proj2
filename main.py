import numpy as np
import random
import networkx
import matplotlib.pyplot as plt
#cyptoscape / networkx

FILE = "./ppi.txt"
LABELS = "./protein_labels.txt"
GAMMA = 0.3
ITLIM = 10000000000



def calc_deltas(scoring, bef):
    for i in range(len(scoring)):
        if scoring[i] == 0 and GAMMA != 1:
            return False
        delt = abs(scoring[i] - bef[i])
        if delt > 0.0001: #how big should our delta be?
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
    iterm = 0

    total_travel = 0
    while True:
        prob_sum = 0
        iterm += 1
        neighbors = [0] * len(graph)
        total_travel += 1
        count[curr] += 1
        before[curr] = score[curr]
        score[curr] = count[curr] / total_travel
        print(score)

        #function here to handle how we stop (calculate delta values in all of score)
        # and stop when they stabilize past a threshold

        
        if iterm % 100:
            if calc_deltas(score, before):
                break
        
        if random.random() <= gamma:
            curr = node
            print("restart")    

        # check what our neighbors are and their scores
        for i in range(len(graph)):
            prev = prob_sum
            if graph[curr, i]:
                prob_sum += graph[curr, i]
                neighbors[i] = ( i, prev, prob_sum)

        choice_taken = random.randint(0 , prob_sum)
        print("choice" , choice_taken , curr)

        for i in range(len(neighbors)):
            if (neighbors[i] and choice_taken <= neighbors[i][2] 
                and choice_taken > neighbors[i][1]):
                curr = neighbors[i][0]
                print(curr , neighbors[i])
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
            print(s)
            labels.append(dic.get(s[0]))

    matrix = np.zeros((count, count))
    for i in edges:
        l = i.split(' ')
        l[1] = l[1].strip('\n')
        first = dic.get(l[0])
        sec = dic.get(l[1])
        print(first, sec, l)
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
    print("degree: ", res1, "\nscr: ", res2 )
    # list will be sorted by (degree, nodenum)
    
# shortest path alg, takes an adjacency matrix and a start and end node
# and returns the shortest path between two nodes
def shortest_path():
    pass

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


def main():
    res = generate_matrix(FILE)
    matr = res[0]
    name = res[1]
    labels = res[2] 
    mat = base_matrix()
    scr = rand_walk_restarts(mat, 1, GAMMA)
    g = networkx.from_numpy_array(mat)
    pos = [np.array([-0.4, 0.5]), np.array([-0.1, 0.3]), 
           np.array([-0.3, -0.3]), np.array([-0.5, -0.3]),
           np.array([-0.4, -0.6]), np.array([0.1, 0.3]),
           np.array([0.1, -0.1]), np.array([0.3, -0.1]),
           np.array([0.3, -0.3]), np.array([0.1, -0.3]),
           np.array([0.2, -0.7])]
    networkx.draw(g,  pos, node_color=scr, linewidths=15, cmap=plt.cm.YlOrRd)
    plt.show()
    print(scr)
    ranking(mat, scr)
    print("done")


if __name__ == '__main__': 
    main()