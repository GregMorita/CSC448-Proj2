import numpy as np
import random
#cyptoscape / networkx

def calc_deltas(scoring, travel_time):
    for i in scoring:
        if i == 0:
            return False
    return True

# Random walks with restarts, takes a graph, starting node, and gamma constant
def rand_walk_restarts(graph, node, gamma):
    #have some gamma value between 0 and 1 for restart chance,
    # want some value for where we go next
    curr = node
    # each 'step' in our case is simply an iteration in our while loop
    # we go until we break bc of deltas of nodes stabilizing
    score = [0] * len(graph[0])
    total_travel = 0
    while True:
        prob_sum = 0
        choices = []
        total_travel += 1
        score[curr] += 1
        #function here to handle how we stop (calculate delta values in all of score)
        if calc_deltas(score, total_travel):
            break

        if random.random() <= gamma:
            curr = node

        for i in range(graph[curr]):
            before = prob_sum
            prob_sum += graph[curr][i]
            choices[i] = (before , prob_sum)

        choice_taken = random.randint(0 , prob_sum)
        for i in range(choices):
            if choice_taken < choices[i][2] and choice_taken > choices[i][1]:
                curr = i
        

        #do something until we break
        pass
    pass

def generate_matrix():
    pass

def main():
    #generate_matrix()
    pass

if __name__ == '__main__': 
    main()