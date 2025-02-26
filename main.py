import numpy as np
import random
#cyptoscape / networkx

FILE = "./ppi.txt"
GAMMA = 0.5


def calc_deltas(scoring, bef):
    for i in range(len(scoring)):
        if scoring[i] == 0:
            return False
        delt = abs(scoring[i] - bef[i])
        if delt > 0.01:
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
    count = [0] * len(graph[0])
    score = [0] * len(graph[0])
    before = [0] * len(graph[0])

    total_travel = 0
    while True:
        prob_sum = 0
        neighbors = []
        total_travel += 1
        count[curr] += 1
        before[curr] = score[curr]
        score[curr] = count[curr] / total_travel

        #function here to handle how we stop (calculate delta values in all of score)
        # and stop when they stabilize past a threshold
        if calc_deltas(score, before):
            break
        if random.random() <= gamma:
            curr = node

        # check what our neighbors are and their scores
        for i in range(graph[curr]):
            before = prob_sum
            if graph[curr][i]:
                prob_sum += graph[curr][i]
                neighbors[i] = (before , prob_sum)

        choice_taken = random.randint(0 , prob_sum)
        for i in range(neighbors):
            #FLAWED LOGIC here FIX LATER PLEASE!!!
            if choice_taken <= neighbors[i][2] and choice_taken > neighbors[i][1]:
                curr = i
    return score

# creates an adjacency matrix from a file representing an edgelist
# returns the generated adjacency matrix and a dictionary mapping protein names to
# indices in the matrix
def generate_matrix(file):
    dic = {}
    count = 0
    with open(file, 'r') as data:
        edges = data.readlines()
        for i in edges:
            s = i.split('\t')
            if s[0] not in dic:
                dic.update({s[0] : count})
                count += 1
            if s[1] not in dic:
                dic.update({s[1] : count})
                count += 1
        
    matrix = np.zeros((count, count))
    for i in edges:
        l = i.split('\t')
        first = dic.get(l[0])
        sec = dic.get(l[1])
        matrix[first][sec] = 1
        matrix[sec][first] = 1
    
    return (matrix , dic)
            

# create the graph as shown in class as an adjacency matrix
def base_matrix():
    matrix = np.matrix([0,1,0,0,0,0,0,0,0,0,0], 
                       [0,0,1,1,0,1,1,0,0,0,0], 
                       [0,1,0,0,1,0,0,0,0,0,0], 
                       [0,1,0,0,1,0,0,0,0,0,0], 
                       [0,0,1,1,0,0,0,0,0,0,0], 
                       [0,1,0,0,0,0,1,1,0,0,0], 
                       [0,1,0,0,0,1,0,0,0,1,0], 
                       [0,0,0,0,0,1,0,0,1,0,0], 
                       [0,0,0,0,0,0,0,1,0,0,1], 
                       [0,0,0,0,0,0,1,0,0,0,1], 
                       [0,0,0,0,0,0,0,0,1,1,0])
    return matrix


def main():
    #res = generate_matrix(FILE)
    #matr = res[0]
    #name = res[1]
    # 
    mat = base_matrix()
    scr = rand_walk_restarts(mat, 1, GAMMA)


if __name__ == '__main__': 
    main()