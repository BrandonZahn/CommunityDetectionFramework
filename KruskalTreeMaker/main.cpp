#include <fstream>
#include <vector>
#include <iostream>
#include <string>

#include "graph.h"

using namespace std;

int main(void)
{
  int numVert = 0;
  int numEdges = 0;
  vector<vector<int>> adjMat;
  vector<vector<float>> weightMat;

  cout << "Getting input file..." << endl;
  ifstream inputFile("input.txt");
  string inputData;

  getline(inputFile, inputData);
  numVert = atoi(inputData.c_str());

  for (int i = 0; i < numVert; i++)
  {
    cout << "AdjMat " << i << endl;
    getline(inputFile, inputData);
    string delimiter = " ";
    size_t position = 0;
    string token;
    vector<int> tempVec;

    while ((position = inputData.find(delimiter)) != inputData.npos)
    {
      token = inputData.substr(0, position);
      tempVec.push_back(atoi(token.c_str()));
      inputData.erase(0, position + delimiter.length());
    }
    tempVec.push_back(atoi(inputData.c_str()));

    adjMat.push_back(tempVec);
  }

  getline(inputFile, inputData);

  for (int i = 0; i < numVert; i++)
  {
    cout << "WeigMat " << i << endl;
    getline(inputFile, inputData);
    string delimiter = " ";
    size_t position = 0;
    string token;
    vector<float> tempVec;

    while ((position = inputData.find(delimiter)) != inputData.npos)
    {
      token = inputData.substr(0, position);
      tempVec.push_back((float) atof(token.c_str()));
      inputData.erase(0, position + delimiter.length());
    }
    tempVec.push_back((float) atof(inputData.c_str()));

    weightMat.push_back(tempVec);
  }

  for (int i = 0; i < numVert; i++)
  {
    for (int j = 0; j < numVert; j++)
    {
      if (i < j)
      {
        if (adjMat[i][j] = 1)
          numEdges++;
      }
    }
  }

  Graph g(numVert, numEdges);

  for (int i = 0; i < numVert; i++)
  {
    for (int j = 0; j < numVert; j++)
    {
      if (i < j)
      {
        if (adjMat[i][j] = 1)
          g.addWeightedEdge(i, j, weightMat[i][j]);
      }
    }
  }

  int mst_wt = g.kruskalMST();

  ofstream outputFile("kruskalTree.txt");
  for (int i = 0; i < numVert; i++)
  {
    for (int j = 0; j < (numVert - 1); j++)
    {
      outputFile << g.kruskalAdjMatrix[i][j] << " ";
    }
    outputFile << g.kruskalAdjMatrix[i][numVert - 1] << endl;
  }

  int end = 0;
  cout << "END?" << endl;
  cin >> end;
  return 0;
}