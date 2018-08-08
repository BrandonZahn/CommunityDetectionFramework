#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

//typedef pair<int, int> iPair;

// To represent Disjoint Sets
struct DisjointSets
{
  int *parent, *rnk;
  int n;

  // Constructor.
  DisjointSets(int n)
  {
    // Allocate memory
    this->n = n;
    parent = new int[n + 1];
    rnk = new int[n + 1];

    // Initially, all vertices are in
    // different sets and have rank 0.
    for (int i = 0; i <= n; i++)
    {
      rnk[i] = 0;

      //every element is parent of itself
      parent[i] = i;
    }
  }

  // Find the parent of a node 'u'
  // Path Compression
  int find(int u)
  {
    /* Make the parent of the nodes in the path
    from u--> parent[u] point to parent[u] */
    if (u != parent[u])
      parent[u] = find(parent[u]);
    return parent[u];
  }

  // Union by rank
  void merge(int x, int y)
  {
    x = find(x), y = find(y);

    /* Make tree with smaller height
    a subtree of the other tree  */
    if (rnk[x] > rnk[y])
      parent[y] = x;
    else // If rnk[x] <= rnk[y]
      parent[x] = y;

    if (rnk[x] == rnk[y])
      rnk[y]++;
  }
};

class Graph
{
public:
  int numVertices;
  int numEdges;
  vector<vector<int>>* adjLists;
  vector<vector<int>> kruskalAdjMatrix;
  vector<pair<float, pair<int, int> > > edges;
  vector<bool>* visited;

  Graph(int V);
  Graph(int V, int E);
  ~Graph();
  void addEdge(int src, int dest);
  void addWeightedEdge(int u, int v, float w);
  int kruskalMST();
  void DFS(int vertex);
  bool checkConnected();
  int GetVisitedCount();
};

Graph::Graph(int vertices)
{
  numVertices = vertices;
  adjLists = new vector<vector<int>>();
  for (int i = 0; i < vertices; i++)
    adjLists->push_back(vector<int>(0));
  visited = new vector<bool>(vertices, false);
}

Graph::Graph(int vertices, int edges)
{
  numVertices = vertices;
  numEdges = edges;
  adjLists = new vector<vector<int>>();
  for (int i = 0; i < vertices; i++)
    adjLists->push_back(vector<int>(0));
  for (int i = 0; i < vertices; i++)
    kruskalAdjMatrix.push_back(vector<int>(vertices, 0));
  visited = new vector<bool>(vertices, false);
}

Graph::~Graph()
{
  delete visited;
  delete adjLists;
  adjLists = NULL;
  visited = NULL;
}

void Graph::addEdge(int src, int dest) { adjLists->at(src).push_back(dest); }

void Graph::addWeightedEdge(int src, int dest, float weight) { edges.push_back({ weight,{ src, dest } }); }

int Graph::kruskalMST()
{
  int mst_wt = 0; // Initialize result
  int edgeCount = 0;
  // Sort edges in increasing order on basis of cost
  sort(edges.begin(), edges.end());

  // Create disjoint sets
  DisjointSets ds(numVertices);

  // Iterate through all sorted edges
  vector< pair<float, pair<int, int> > >::iterator it;
  for (it = edges.begin(); it != edges.end(); it++)
  {

    int u = it->second.first;
    int v = it->second.second;

    int set_u = ds.find(u);
    int set_v = ds.find(v);

    cout << "DFS with " << GetVisitedCount() << " visited. ";
    DFS(1850 % 50);

    // Check if the selected edge is creating
    // a cycle or not (Cycle is created if u
    // and v belong to same set)
    if (!checkConnected())
    {
      // Current edge will be in the MST
      cout << "Adding " << edgeCount << " edge to kruskal: " << u << " - " << v << endl;
      kruskalAdjMatrix.at(u).at(v) = 1;
      kruskalAdjMatrix.at(v).at(u) = 1;

      // Update MST weight
      mst_wt += it->first;

      // Merge two sets
      ds.merge(set_u, set_v);
      edgeCount++;
      //cout << edgeCount << endl;
    }
    else
    {
      kruskalAdjMatrix.at(u).at(v) = 0;
      kruskalAdjMatrix.at(v).at(u) = 0;
    }
  }

  return mst_wt;
}

void Graph::DFS(int vertex)
{
  visited->at(vertex) = true;
  if (kruskalAdjMatrix.at(vertex).size() > 0)
    for (size_t i = 0; i < kruskalAdjMatrix.at(vertex).size(); i++)
      if (!visited->at(kruskalAdjMatrix.at(vertex).at(i)))
        DFS(kruskalAdjMatrix.at(vertex).at(i));
}

bool Graph::checkConnected()
{
  bool connected = false;
  int visitedCount = 0;
  for (int i = 0; i < numVertices; i++)
  {
    if (visited->at(i) == true)
      visitedCount++;
  }
  if (visitedCount == numVertices)
    connected = true;
  return connected;
}

int Graph::GetVisitedCount()
{
  int visitedCount = 0;
  for (int i = 0; i < numVertices; i++)
  {
    if (visited->at(i))
      visitedCount++;
  }
  return visitedCount;
}