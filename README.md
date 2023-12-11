# Minimum-FeedBack-Arc-Set

## Overview

This program is used to find the minimum feedback arc set, that is, to break cycles in the following 4 types of graphs with minimum cost(number of edges removed) if possible: 
        1. undirected unweighted graph 
        2. undirected weighted graph
        3. directed unweighted graph
        4. directed weighted graph

Proposed solution:
    1. For an undirected unweighted graph, we can use DFS to find and break cycles. Time complexity is O(V + E).
    2. For an undirected weighted graph, we can reverse all edges' weight and then use Kruskal's Algorithm to find the MST, then 
    the remaining edges are the minimum feedback arc set. Time complexity is $O(E log E + E \alpha(V))$ = O(E log E), where $\alpha(V)$ is the Ackermann function used in the disjoint set data structure.
    3. For a directed unweighted graph and a directed weighted graph, the problem becomes NP-hard, so we adopt a greedy approximation algorithm. 
    
Source of the abovementioned greedy approximation algorithm: Simpson, Michael, Venkatesh Srinivasan, and Alex Thomo. "Efficient computation of feedback arc set at web-scale." Proceedings of the VLDB Endowment 10.3 (2016): 133-144.


## Usage

1. `make`
2. `./cb <path_to_your_input_file> <path_to_your_output_file> <br>`
