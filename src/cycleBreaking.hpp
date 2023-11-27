#pragma once

#include <vector>       // vector
#include <utility>      // pair 
#include <tuple>        // tuple

using namespace std;

struct Edge
{
    int src, dest, weight;
};


class Graph
{
private:
    int V; // number of vertices
    vector<Edge> edges;  // all edges
    vector<vector<pair<int, int>>> adj;  // adjacency list
    vector<int> inDegree; // for greedy_FAS
    vector<int> outDegree; // for greedy_FAS
    vector<int> color; // dfs traversal, detect cycle
public:
    Graph(int V);
    void add_undirected_edge(int u, int v, int w);
    void add_directed_edge(int u, int v, int w);
    vector<int> find_cycle();
    void reverse_graph();
    bool remove_edge(int u, int v);
    void print_graph();
    void dfs_pruning(vector<tuple<int, int, int>>& fas, int& cost);
    void dfs_pruning_util(vector<tuple<int, int, int>>& fas, int& cost, int node, int parent);
    void kruskal_MST(vector<tuple<int, int, int>>& fas, int& cost);
    void greedy_FAS(vector<tuple<int, int, int>>& fas, int& cost);
};


class DisjointSet
{
private:
    int V;
    vector<int> parent;
    vector<int> rank;
public:
    DisjointSet(int V);
    bool Union(int x, int y);
    int find(int x);
};


