// Got this from https://www.geeksforgeeks.org/connected-components-in-an-undirected-graph/
#include"graphs.h"
#include<iostream>
#include<vector>
using namespace std;
extern double cellvol;
extern int maxclustersize;
extern double M;
int clustersize, numcomp;
extern vector<int> clustercelllabel, distinctclusternum;
extern int ncell[3];

// Method to print connected components in an undirected graph
int Graph::connectedComponents() {
    int v;
     numcomp = 0;

// Mark all the vertices as not visited
    bool* visited = new bool[V];
    for (v = 0; v < V; v++) visited[v] = false;
//    cout << "V = " << V << endl;
    

    for (v = 0; v < V; v++) if (visited[v] == false) {
            // print all reachable vertices
            // from v
            clustersize = 0;
            DFSUtil(v, visited);
//            cout << "\n";
            numcomp++;
            if (clustersize > maxclustersize) maxclustersize = clustersize;
  
           M += clustersize*cellvol /numcomp;
    }
      
//    cout << "Number of connected components = " << numcomp << endl;
    delete[] visited;
    
    return numcomp;
}


// depth-first search
// Need figure out how to get the connected-component sizes
void Graph::DFSUtil(int v, bool visited[]) {
    // Mark the current node as visited and print it
    visited[v] = true;
// print the next connected component
//    cout << v << " ";
    clustersize++;  // increment the size of the current void
   
    distinctclusternum[clustercelllabel[v]] = numcomp;
    // Recur for all the vertices adjacent to this vertex
    list<int>::iterator i;
    for (i = adj[v].begin(); i != adj[v].end(); ++i)
        if (!visited[*i]) {
//            voidsize++;
            DFSUtil(*i, visited);
        }
}


// Constructor
Graph::Graph(int V) {
    this->V = V;
    adj = new list<int>[V];
}


// Destructor
Graph::~Graph() { delete[] adj; }


// method to add an undirected edge
void Graph::addEdge(int v, int w)
{
    adj[v].push_back(w);
    adj[w].push_back(v);
}






