// This was taken from https://www.geeksforgeeks.org/connected-components-in-an-undirected-graph/
#ifndef GRAPHS_H
#define GRAPHS_H
#endif
#include <list>
using namespace std;

// Graph class represents a undirected graph
// using adjacency list representation
class Graph {
    int V; // No. of vertices

    // Pointer to an array containing adjacency lists
    list<int>* adj;

    // A function used by DFS
    void DFSUtil(int v, bool visited[]);

public:
    Graph(int V); // Constructor
    ~Graph();
    void addEdge(int v, int w);
    int connectedComponents();
//   void connectedComponents();
};

