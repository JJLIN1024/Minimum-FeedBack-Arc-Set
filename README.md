# Minimum-FeedBack-Arc-Set

## Overview

This program is used to break cycle, that is, to find the minimum feedback arc set in the following 4 types of graphs: 
        1. undirected unweighted graph 
        2. undirected weighted graph
        3. directed unweighted graph
        4. directed weighted graph

Proposed solution:
    1. For undirected unweighted graph, we can use DFS to find all cycles, and break it arbitrary. Time complexity is O(V + E).
    2. For undirected weighted graph, we can reverse all edges' weight and then use Kruskal's Algorithm to find the MST, then 
    the remaining edges are the minimum feed back arc set. Time complexity is $O(E log E + E \alpha(V))$ = O(E log E), where $\alpha(V)$ is the Ackermann function used in disjoint set data structure.
    3. For directed unweighted graph and directed weighted graph, the problem becomes NP-hard, so we adopt greedy approximation algorithm. 
    
Source of abovementioned greedy approximation algorithm: Simpson, Michael, Venkatesh Srinivasan, and Alex Thomo. "Efficient computation of feedback arc set at web-scale." Proceedings of the VLDB Endowment 10.3 (2016): 133-144.


## Usage

1. `make`
2. `./cb <path_to_your_input_file> <path_to_your_output_file> <br>`

