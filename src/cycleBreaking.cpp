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
    For undirected unweighted graph, we can use DFS to find all cycles, and break it arbitrary. For undirected weighted graph, we
    can reverse all edges' weight and then use Kruskal's Algorithm to find the MST, then the remaining edges are the minimum feed back
    arc set. For directed unweighted graph and directed weighted graph, the problem becomes NP-hard, so we adopt approximation algorithm.
*/

#include <iostream>
#include <cstring>
#include <fstream>
#include <vector>
#include <string>
#include <list>
#include <queue>
#include <iomanip>
#include <utility>
#include <tuple>
#include <set>
#include <unordered_set>
#include <iterator>
#include <queue>
#include <algorithm>


using namespace std;
// A structure to represent a subset for union-find
class subset
{
    public:
    int parent;
    int rank;
};

// class Edge
// {
//     public:
//     int src, dest, weight;
// };

// A utility function to find set of an element i
// (uses path compression technique)
int find(subset subsets[], int i)
{
    // find root and make root as parent of i
    // (path compression)
    if (subsets[i].parent != i)
        subsets[i].parent = find(subsets, subsets[i].parent);

    return subsets[i].parent;
}

// A function that does union of two sets of x and y
// (uses union by rank)
void Union(subset subsets[], int x, int y)
{
    int xroot = find(subsets, x);
    int yroot = find(subsets, y);

    // Attach smaller rank tree under root of high
    // rank tree (Union by Rank)
    if (subsets[xroot].rank < subsets[yroot].rank)
        subsets[xroot].parent = yroot;
    else if (subsets[xroot].rank > subsets[yroot].rank)
        subsets[yroot].parent = xroot;

    // If ranks are same, then make one as root and
    // increment its rank by one
    else
    {
        subsets[yroot].parent = xroot;
        subsets[xroot].rank++;
    }
}


class Graph
{
    int V;
    vector<tuple<int, int, int>> *adj;
    vector<int> *inDegree;
    vector<int> *outDegree;
    public:
    Graph(int V);

    void undirectedAddEdge(int u, int v, int w);
    void directedAddEdge(int u, int v, int w);
    int findEdgeWeight(int u, int v);
    int findRealEdgeWeight(int u, int v);
    bool removeEdge(int u, int v, int w);
    bool decreaseEdge(int u, int v, int w);
    tuple<vector<tuple<int, int, int>>, int> kruskalMST();
    void FAS(vector<tuple<int, int, int>>& answerList);
    void feedbackEdgeSetUtil(int v, int parent, int color[], vector<tuple<int, int, int>>& circles, int parentList[], int& cyclenumber);
    void feedbackEdgeSet(vector<tuple<int, int, int>>& circles);
    bool directedDetectCycleUtil(int v,int color[]);
    bool directedDetectCycle();
    void printDegree();
};

Graph::Graph(int V){
    this->V = V;
    adj = new vector<tuple<int, int, int>>[V];
    inDegree = new vector<int>(V, 0);
    outDegree = new vector<int>(V, 0);
}

void Graph::undirectedAddEdge(int u, int v, int w) {
    adj[u].push_back(tuple<int, int, int>(-w, v, w));
    // adj[v].push_back(tuple<int, int, int>(-w, u, w));

}

void Graph::directedAddEdge(int u, int v, int w) {
    adj[u].push_back(tuple<int, int, int>(w, v, w));
    inDegree->at(v) = inDegree->at(v) + 1;
    outDegree->at(u) = outDegree->at(u) + 1;
}

void Graph::printDegree() {
    cout << "Indegree :\n";
    vector<int>::iterator i;
    for(i = inDegree->begin(); i != inDegree->end(); i++) {
        cout << *i << " ";
    }
    cout << "\n";

    cout << "outdegree :\n";
    vector<int>::iterator j;
    for(j = outDegree->begin(); j != outDegree->end(); j++) {
        cout << *j << " ";
    }
    cout << "\n";
}

int Graph::findEdgeWeight(int u, int v) {

    vector<tuple<int, int, int>> :: iterator i;
    for (i = adj[u].begin(); i != adj[u].end(); i++) {
        int vertex = get<1>(*i);
        int weight = get<0>(*i);
        if (v == vertex)
        {
            return weight;
        }
    }
    return -1;
}

int Graph::findRealEdgeWeight(int u, int v) {

    vector<tuple<int, int, int>> :: iterator i;
    for (i = adj[u].begin(); i != adj[u].end(); i++) {
        int vertex = get<1>(*i);
        int weight = get<2>(*i);
        if (v == vertex)
        {
            return weight;
        }
    }
    return -1;
}

bool Graph::removeEdge(int u, int v, int w){
    vector<tuple<int, int, int>> :: iterator i;
    for (i = adj[u].begin(); i != adj[u].end(); i++) {
        int vertex = get<1>(*i);
        if (v == vertex)
        {
            adj[u].erase(i);
            return true;
        }
    }
    return false;
}

bool Graph::decreaseEdge(int u, int v, int w){
    vector<tuple<int, int, int>> :: iterator i;
    for (i = adj[u].begin(); i != adj[u].end(); i++) {
        int vertex = get<1>(*i);
        int weight = get<0>(*i);
        if (v == vertex)
        {
            int weightToBe = weight - w;
            get<0>(*i) = weightToBe;
            return true;
        }
    }
    return false;
}

tuple<vector<tuple<int, int, int>>, int> Graph::kruskalMST(){

    vector<tuple<int, int, int>> Answer;

    //collect all the edges
    vector<tuple<int, int, int>> Edges;
    for(int j = 0; j < V; j++) {
        vector<tuple<int, int, int>>::iterator i;
        for (i = adj[j].begin(); i != adj[j].end(); ++i)
        {
            int v = get<1>(*i);
            int w = get<0>(*i);
            Edges.push_back(make_tuple(w, j, v));
        }
    }
    sort(Edges.begin(), Edges.end());
//    reverse(Edges.begin(), Edges.end());

    // Allocate memory for creating V subsets
    subset *subsets = new subset[( V * sizeof(subset) )];

    // Create V subsets with single elements
    for (int v = 0; v < V; ++v)
    {
        subsets[v].parent = v;
        subsets[v].rank = 0;
    }
    int e = 0; // An index variable, used for result[]
    int i = 0; // An index variable, used for sorted edges
    // Number of edges to be taken is equal to V-1
    int numberOfEdges = Edges.size();
    while (e < V - 1 && i < numberOfEdges)
    {
//        (w, j, v)
        tuple<int, int, int> nextEdge = Edges[i];
        i++;
//        int w = get<0>(nextEdge);
        int u = get<1>(nextEdge);
        int v = get<2>(nextEdge);

        int x = find(subsets, u);
        int y = find(subsets, v);

        if (x != y)
        {
            Answer.push_back(nextEdge);
            e++;
            Union(subsets, x, y);
        }

    }
    // cout<<"Following are the edges in the constructed MST\n";
    for(int i = 0; i < Answer.size(); i++) {
        int w = get<0>(Answer[i]);
        int u = get<1>(Answer[i]);
        int v = get<2>(Answer[i]);
//        cout << u << " " << v << " " << w << endl;
        removeEdge(u, v, w);
    }

    int totalCost = 0;
    vector<tuple<int, int, int>> finalAnswer;
    for(int j = 0; j < V; j++) {
        vector<tuple<int, int, int>>::iterator i;
        for (i = adj[j].begin(); i != adj[j].end(); ++i)
        {
            int weight = get<0>(*i);
            int v = get<1>(*i);
//            removeEdge(v, j, weight);
            totalCost += -weight;
            finalAnswer.push_back(tuple<int, int, int>(j, v, -weight));
        }
    }

    tuple<vector<tuple<int, int, int>>, int> A(finalAnswer, totalCost);
    return A;

}

void Graph::FAS(vector<tuple<int, int, int>> &answerList){

    sort(adj[0].begin(), adj[0].end());
    // turn edge weight to R+
    for(int m = 0; m < V; m++){
        vector<tuple<int, int, int>>::iterator n;
        for(n = adj[m].begin(); n != adj[m].end(); n++) {
            int w = get<0>(*n);
            int newWeight = w + 101;
            get<0>(*n) = newWeight;
        }
    }

    while(directedDetectCycle() == true) {
        vector<tuple<int, int, int>> circles;
        feedbackEdgeSet(circles);
        // can't find cycle anymore
        if (circles.size() == 0){
            cout << "siez is 0! " << endl;
            break;
        }
        vector<tuple<int, int, int>>::iterator j;
        for(j = circles.begin(); j != circles.end(); j++) {
            answerList.push_back(*j);
        }
        for(int i = 0; i < circles.size(); i++) {
            int u = get<0>(circles[i]);
            int v = get<1>(circles[i]);
            int w = get<2>(circles[i]);
            removeEdge(u, v, w);
        }
    }
}
void Graph::feedbackEdgeSet(vector<tuple<int, int, int>> &circles){

    // color == 0 means this vertex has not been visited yet.
    // color == 1 means partially visited.
    // color == 2 means complete vistited.

    int parentList[V];
    int cyclenumber = 0;
    int *color = new int[V];
    int *inCycle = new int[V];
    for(int i = 0; i < V; i++) {
        color[i] = 0;
        inCycle[i] = 0;
    }
    for(int i = 0; i < V; i++) {
        if (color[i] == 0) {
            feedbackEdgeSetUtil(i, i, color, circles, parentList, cyclenumber);
        }
    }

}


void Graph::feedbackEdgeSetUtil(int v, int parent, int color[], vector<tuple<int, int, int>>& circles, int parentList[], int& cyclenumber)
{
    if (color[v] == 2) {
//        cout << "finish: " << v << endl;
        return;
    }
    // we find a back edge
    if (color[v] == 1) {
        // cout << "find it: " << v << endl;
        cyclenumber ++;
        int cur = parent;

        // backtrack the vertex which are
        // in the current cycle thats found
        vector<int> temp;
        temp.push_back(v);
        while (cur != v) {
            // cout << "cur " << cur << endl;
            temp.push_back(cur);
            cur = parentList[cur];
        }

        temp.push_back(v);
        reverse(temp.begin(), temp.end());
        int minWeight = 999;
        int start = -1;
        int end = -1;
        for(int k = 0; k < temp.size() - 1; k++) {
            int w = findEdgeWeight(temp[k], temp[k+1]);
            if( w <= minWeight and w != 0) {
                minWeight = w;
                start = temp[k];
                end = temp[k+1];
            }
        }

        int realWeight = findRealEdgeWeight(start, end);
        circles.push_back(tuple<int, int, int> (start, end, realWeight));
        return;
    }

    parentList[v] = parent;
    color[v] = 1;

    // Recur for all the vertices adjacent to this vertex
    vector<tuple<int, int, int>>::iterator i;
    for (i = adj[v].begin(); i != adj[v].end(); ++i)
    {
        // If an adjacent is not visited, then recur for that adjacent

        int vertex = get<1>(*i);
        // if (vertex == parentList[v]) {
        //     continue;
        // }
        feedbackEdgeSetUtil(vertex, v, color, circles, parentList, cyclenumber);

    }
    color[v] = 2;
}


bool Graph::directedDetectCycle(){

    // color == 0 means this vertex has not been visited yet.
    // color == 1 means partially visited.
    // color == 2 means complete vistited.

    int *color = new int[V];
    for(int i = 0; i < V; i++) {
        color[i] = 0;
    }
    for(int i = 0; i < V; i++) {
        if (color[i] == 0) {
            if (directedDetectCycleUtil(i, color) == true){
                // cout << "graph contain cycle! " << endl;
                return true;
            }
        }
    }
    // cout << "graph doesn't contain cycle! " << endl;
    return false;

}

bool Graph::directedDetectCycleUtil(int v, int color[])
{


    color[v] = 1;

    // Recur for all the vertices adjacent to this vertex
    vector<tuple<int, int, int>>::iterator i;
    for (i = adj[v].begin(); i != adj[v].end(); ++i)
    {

        int vertex = get<1>(*i);
        if (color[vertex] == 1) {
            return true;
        }
        if (color[vertex] == 0 and directedDetectCycleUtil(vertex, color)) {
            return true;
        }
    }
    color[v] = 2;
    return false;
}

int main(int argc, const char * argv[]) {

    if (argc < 3) {
        cerr << "Error: command format should be --> ./mps <input file name> <output file name>" << endl;
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

    if (direction == 'u') {
        for (int i = 0; i < totalEdges; i++) {
            int u, v, w;
            infile >> u;
            infile >> v;
            infile >> w;
            g.undirectedAddEdge(u, v, w);
        }
    }

    else if (direction == 'd') {
        for (int i = 0; i < totalEdges; i++) {
            int u, v, w;
            infile >> u;
            infile >> v;
            infile >> w;
            g.directedAddEdge(u, v, w);
        }
    }

    infile.close();
    if (direction == 'u') {
        tuple<vector<tuple<int, int, int>>, int> A;
        // cout << "start of undirected." << endl;
        A = g.kruskalMST();

        ofstream outfile;
        outfile.open(argv[2], ios::out);
        vector<tuple<int, int, int>> finalAnswer;
        int totalCost;
        finalAnswer = get<0>(A);
        totalCost = get<1>(A);

        // cout << "totalCost: " << totalCost << endl;
        outfile << totalCost << "\n";

        // the input file is already a DAG
        if (totalCost == 0) {
            return 0;
        }

        for(int i = 0; i < finalAnswer.size(); i++) {
            int u = get<0>(finalAnswer[i]);
            int v = get<1>(finalAnswer[i]);
            int w = get<2>(finalAnswer[i]);
            outfile << u << " " << v << " " << w << "\n";
        }
    }

    else if (direction == 'd') {
        // cout << "start of directed." << endl;
//        g.directedDetectCycle();

        vector<tuple<int, int, int>> answerList;
        g.FAS(answerList);
        // erase duplicates
        set<tuple<int, int, int>> s;
        unsigned size = answerList.size();
        for( unsigned i = 0; i < size; ++i ) s.insert( answerList[i] );
        answerList.assign( s.begin(), s.end() );

        int totalCost2 = 0;
        for(int i = 0; i < answerList.size(); i++) {
            int w = get<2>(answerList[i]);
            totalCost2 += w;
        }

        ofstream outfile;
        outfile.open(argv[2], ios::out);

        // cout << "totalCost2: " << totalCost2 << endl;
        outfile << totalCost2 << "\n";

        if (totalCost2 == 0) {
            return 0;
        }

        for(int i = 0; i < answerList.size(); i++) {
            int u = get<0>(answerList[i]);
            int v = get<1>(answerList[i]);
            int w = get<2>(answerList[i]);
            outfile << u << " " << v << " " << w << "\n";
            // cout <<  u << " " << v << " " << w << "\n";
        }
    }
    return 0;
}



