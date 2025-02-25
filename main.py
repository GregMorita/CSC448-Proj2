import numpy as np
#cyptoscape / networkx

# Random walks with restarts 
def rand_walk_restarts(graph, node, gamma):
    #have some gamma value between 0 and 1 for restart chance,
    curr = node
    while True:
        sum = 0
        choices = {}
        for i in range(graph[curr]):
            sum += graph[curr][i]
            choices.update({i : graph[curr][i]})
        

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