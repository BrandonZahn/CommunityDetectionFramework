#include <iostream>
#include <vector>
using namespace std;

class Graph
{
public:
	int numVertices;
	vector<vector<int> >* adjLists;
	vector<bool>* visited;

	Graph(int V);
	~Graph();
	void addEdge(int src, int dest);
	void DFS(int vertex);
	bool checkConnected();
	int GetVisitedCount();
};

Graph::Graph(int vertices)
{
	numVertices = vertices;
	adjLists = new vector<vector<int> >();
	for (int i = 0; i < vertices; i++)
	{
		adjLists->push_back(vector<int>(0));
	}
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

void Graph::DFS(int vertex)
{
	visited->at(vertex) = true;
	if (adjLists->at(vertex).size() > 0)
	{
		for (size_t i = 0; i < adjLists->at(vertex).size(); i++)
		{
			if (!visited->at(adjLists->at(vertex).at(i)))
				DFS(adjLists->at(vertex).at(i));
		}
	}
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