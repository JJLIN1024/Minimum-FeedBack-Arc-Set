//  Created by Jimmy Lin on 2019/12/14.
//  Last modified by Jimmy Lin on 2023/12/1
//  Copyright © 2023 Jimmy Lin. All rights reserved.

/*

Problem statement: 
    find minimum feedback arc set in the following 4 types of graphs: 
        1. undirected unweighted graph 
        2. undirected weighted graph
        3. directed unweighted graph
        4. directed weighted graph
Proposed solution:
    1. For undirected unweighted graph, we can use DFS to find all cycles, and break it arbitrary. Time complexity is O(V + E).
    2. For undirected weighted graph, we can reverse all edges' weight and then use Kruskal's Algorithm to find the MST, then 
    the remaining edges are the minimum feed back arc set. Time complexity is O(E log E + E alpha(V)) = O(E log E), where alpha(V) is the Ackermann function used 
    in disjoint set data structure.
    3. For directed unweighted graph and directed weighted graph, the problem becomes NP-hard, so we adopt approximation algorithm.

    The equivalent of a minimum spanning tree in a directed graph is called an optimum branching or a minimum-cost arborescence. The classical algorithm for solving this problem is the Chu-Liu/Edmonds algorithm
*/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <string>
#include <unordered_set>

#include "cycleBreaking.hpp"

using namespace std;

DisjointSet::DisjointSet(int V) {
    this->V = V;
    this->parent.resize(V);
    for(int i = 0; i < V; i++) {
        parent[i] = i;
    }

    this->rank.resize(V, 1);
}

bool DisjointSet::Union(int x, int y) {
    x = find(x), y = find(y);
    if(x != y) {
        if(rank[x] > rank[y]) {
            parent[y] = x;
            rank[x] += 1;
        } else {
            parent[x] = y;
            rank[y] += 1;
        }
        return false;
    }
    return true;
}

int DisjointSet::find(int x) {
    if(parent[x] != x) {
        return parent[x] = find(parent[x]);
    }
    return parent[x]; 
}

Graph::Graph(int V) {
    this->V = V;
    this->adj.resize(V);
    this->inDegree.resize(V, 0);
    this->outDegree.resize(V, 0);
    this->color.resize(V, 0);
}

void Graph::add_undirected_edge(int u, int v, int w) {
    adj[u].push_back(make_pair(v, w));
    adj[v].push_back(make_pair(u, w));
    
    // kruskal's algorithm only needs one edge 
    struct Edge e = {u, v, w};
    edges.push_back(e);

    inDegree[v] += w;
    outDegree[u] += w;
    inDegree[u] += w;
    outDegree[v] += w;
}

void Graph::add_directed_edge(int u, int v, int w) {
    adj[u].push_back(make_pair(v, w));
    
    struct Edge e1 = {u, v, w};
    edges.push_back(e1);

    inDegree[v] += w;
    outDegree[u] += w;
}

bool Graph::remove_edge(int u, int v) {

    vector<pair<int, int>>:: iterator i;
    for (i = adj[u].begin(); i != adj[u].end(); i++) {
        int vertex = i->first;
        if (v == vertex)
        {
            adj[u].erase(i);
            return true;
        }
    }
    return false;
}

void Graph::reverse_graph() {
    vector<vector<pair<int, int>>> adj_reverse;
    adj_reverse.resize(this->V);
    for(int i = 0; i < adj.size(); i++) {
        for(int j = 0; j < adj[i].size(); j++) {
            int v = adj[i][j].first;
            int w = adj[i][j].second;
            adj_reverse[v].push_back({i, w});
        }
    }
    this->adj = adj_reverse;
}

vector<int> Graph::find_cycle() {

    unordered_map<int, int> parents;
    unordered_map<int, int> visited;
    unordered_map<int, int> onStack;
    vector<int> stack;
    for(int i = 0; i < this->adj.size(); i++) {
        if(visited[i]) continue;
        stack.push_back(i);
        while(!stack.empty()) {
            int n = stack.back();
            if(!visited[n]) {
                visited[n] = 1;
                onStack[n] = 1;
            } else {
                onStack[n] = 0;
                stack.pop_back();
            }

            for (auto const &neighbor : adj[n]) {
                int v = neighbor.first;
                if (!visited[v]) {
                    stack.push_back(v);
                    parents[v] = n;
                } else if (onStack[v]) {
                    vector<int> cycle = {v};
                    while(n != v) {
                        cycle.push_back(n);
                        n = parents[n];
                        return cycle;
                    }
                }
            }
        }
    }
    return {};
}

void Graph::print_graph() {   
    // neighbors
    for(int i = 0; i < adj.size(); i++) {
        for(int j = 0; j < adj[i].size(); j++) {
            int dest = adj[i][j].first;
            int w = adj[i][j].second;
            cout << i << " " << dest << " " << w << endl;
        }
    }
}

void Graph::dfs_pruning(vector<tuple<int, int, int>>& fas, int& cost) {
    for(int i = 0; i < this->V; i++) {
        if(color[i] == 0) {
            dfs_pruning_util(fas, cost, i, i);
        }
    }
}

void Graph::dfs_pruning_util(vector<tuple<int, int, int>>& fas, int& cost, int node, int parent) {
    color[node] = 1;
    for(auto neighbor: adj[node]) {
        int v = neighbor.first;
        int w = neighbor.second;
        if(v == parent) continue;
        if(color[v] == 0) {
            dfs_pruning_util(fas, cost, v, node);
        } else if (color[v] == 1) {
            cost += w;
            fas.push_back({node, v, w});
        }
    }
    color[node] = 2;
}

void Graph::kruskal_MST(vector<tuple<int, int, int>>& fas, int& cost) {
    DisjointSet d(this->V);
    vector<Edge> edges_transform;
    // negate the weight to find the maximum spanning tree
    for(auto& edge: edges) {
        edges_transform.push_back({edge.dest, edge.src, -edge.weight});
    }
    sort(edges_transform.begin(), edges_transform.end(), [](const Edge& e1, const Edge& e2) {return e1.weight < e2.weight;});
    for(auto& edge: edges_transform) {
        bool formCycle = d.Union(edge.dest, edge.src);
        if(formCycle) {
            fas.push_back({edge.dest, edge.src, -edge.weight}); // negate the weight back to its original value
            cost += -edge.weight;
        }
    }
}

void Graph::greedy_FAS(vector<tuple<int, int, int>>& fas, int& cost) {

    vector<int> s1 = {};
    vector<int> s2 = {};

    vector<int> sinks = {};
    vector<int> sources = {};

    for(int i = 0; i < inDegree.size(); i++) {
        if(inDegree[i] == 0) {
            sources.push_back(i);
        }
    }

    for(int i = 0; i < outDegree.size(); i++) {
        if(outDegree[i] == 0) {
            sinks.push_back(i);
        }
    }
    vector<int> nodes;
    for(int i = 0; i < this->V; i++) {
        nodes.push_back(i);
    }
    while(!nodes.empty()) {
        while(!sinks.empty()) {
            int s = sinks.back();
            sinks.pop_back();
            s1.push_back(s);
            for(int i = 0; i < adj.size(); i++) {
                if(i == s) continue;
                for(int j = 0; j < adj[i].size(); j++) {
                    int v = adj[i][j].first;
                    if(v == s) {
                        remove_edge(i, v);
                        outDegree[i] -= 1;
                        if(outDegree[i] == 0) sinks.push_back(i);
                    }
                }
            }
            auto it = find(nodes.begin(), nodes.end(), s);
            if (it != nodes.end()) {
                nodes.erase(it);
            }

        }
        while(!sources.empty()) {
            int c = sources.back();
            sources.pop_back();
            s2.push_back(c);
            for(int i = 0; i < adj[c].size(); i++) {
                int v = adj[c][i].first;
                remove_edge(c, v);
                inDegree[v] -= 1;
                if(inDegree[v] == 0) sources.push_back(v);
            }
            auto it = find(nodes.begin(), nodes.end(), c);
            if (it != nodes.end()) {
                nodes.erase(it);
            }
        }

        if(!nodes.empty()) {
            int maxDiff = INT_MIN;
            int node_choosen = -1;
            for(int i = 0; i < nodes.size(); i++) {
                int diff = outDegree[nodes[i]] - inDegree[nodes[i]];
                if(diff > maxDiff) {
                    maxDiff = diff;
                    node_choosen = nodes[i];
                }
            }
            if(node_choosen != -1)
                s1.push_back(node_choosen);
            auto it = find(nodes.begin(), nodes.end(), node_choosen);
            if (it != nodes.end()) {
                nodes.erase(it);
            }
        }
    }
    vector<int> s12;
    s12.reserve( s1.size() + s2.size() ); // preallocate memory
    s12.insert( s12.end(), s1.begin(), s1.end() );
    s12.insert( s12.end(), s2.begin(), s2.end() );

    reverse_graph();

    unordered_set<string> st;
    for(int i = 0; i < adj.size(); i++) {
        for(int j = 0; j < adj[i].size(); j++) {
            int v = adj[i][j].first;
            st.insert(to_string(i) + "->" + to_string(v));
        }
    }

    for(int i = 0; i < s12.size(); i++) {
        for(int j = i + 1; j < s12.size(); j++) {
            string possible_inversion = to_string(s12[i]) + "->" + to_string(s12[j]);
            if(st.find(possible_inversion) != st.end()) {
                // there's a inversion! the inversion is one of our feedback arc
                for(int k = 0; k < adj[s12[i]].size(); k++) {
                    if(adj[s12[i]][k].first == s12[j]) {
                        int w = adj[s12[i]][k].second;
                        fas.push_back({s12[j], s12[i], w});
                        cost += w;
                    }
                }
            }
        }
    }
    reverse_graph();
}

int main(int argc, const char * argv[])  {

    if (argc < 3) {
        cerr << "Error: command format should be --> ./cb <input file name> <output file name>" << endl;
        return 1;
    }

    ifstream infile;
    infile.open(argv[1], ios::in);

    if (!infile.is_open()) {
        cerr << "Error opening file!" << endl;
        return 1;
    }

    char direction;
    int totalVertices;
    int totalEdges;

    infile >> direction;
    infile >> totalVertices;
    infile >> totalEdges;

    Graph g(totalVertices);

    bool weighted = false;
    for (int i = 0; i < totalEdges; i++) {
        int u, v, w;
        infile >> u;
        infile >> v;
        infile >> w;
        if(w != 1) {
            weighted = true;
        }
        if(direction == 'u')
            g.add_undirected_edge(u, v, w);
        else if(direction == 'd')
            g.add_directed_edge(u, v, w);
    }
    infile.close();

    vector<tuple<int, int, int>> feed_back_arc_set;
    int cost = 0;
    if (direction == 'u') {
        if(weighted) { // weighted undirected graph
            g.kruskal_MST(feed_back_arc_set, cost);
        } else {
            g.dfs_pruning(feed_back_arc_set, cost);
        }
    } else if (direction == 'd') {
        g.greedy_FAS(feed_back_arc_set, cost);
    }

    for(int i = 0; i < feed_back_arc_set.size(); i++) {
        int u = get<0>(feed_back_arc_set[i]);
        int v = get<1>(feed_back_arc_set[i]);
        g.remove_edge(u, v);
    }

    vector<int> cycle = g.find_cycle();
    if(!cycle.empty()) {
        cout << "there's still cycle to break!" << endl;
    }
    

    // write feed back arc set to output file
    ofstream outfile;
    outfile.open(argv[2], ios::out);

    outfile << cost << "\n";
    for(int i = 0; i < feed_back_arc_set.size(); i++) {
        int u = get<0>(feed_back_arc_set[i]);
        int v = get<1>(feed_back_arc_set[i]);
        int w = get<2>(feed_back_arc_set[i]);
        outfile << u << " " << v << " " << w << "\n";
    }
    return 0;
}



