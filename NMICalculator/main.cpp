#include <vector>
#include <numeric>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <algorithm>

int main(void)
{
  std::string inputData;

	std::ifstream classesFile("classesCora.txt");
  
	std::getline(classesFile, inputData);
  int classCount = stoi(inputData);

	std::vector<std::string> classes;
	std::vector<int> classesRatio;
  for (int i = 0; i < classCount; i++)
  {
    getline(classesFile, inputData);
    classes.push_back(inputData);
		classesRatio.push_back(0);
  }

	std::vector<std::string> nodeClassClassification;
  while (!classesFile.eof())
  {
    getline(classesFile, inputData);
    nodeClassClassification.push_back(inputData);
		for (int i = 0; i < classCount; i++)
		{
			if (inputData == classes[i])
				classesRatio[i]++;
		}
  }
  classesFile.close();	

	std::ifstream clustersFile("clusters.txt");
	
	std::vector<std::string> nodeClusterClassification;
  while (!clustersFile.eof())
  {
    getline(clustersFile, inputData);
    nodeClusterClassification.push_back(inputData);
  }
	
	std::vector<std::pair<std::string, int>> clusters;
	clusters.push_back(make_pair(nodeClusterClassification[0], 1));
  for (int i = 1; i < nodeClusterClassification.size(); i++)
  {
		for (int j = 0; j < clusters.size(); j++)
		{
			if (clusters[j].first == nodeClusterClassification[i])
			{
				clusters[j].second++;
				break;
			}
			else if (j == clusters.size() - 1)
			{
				clusters.push_back(make_pair(nodeClusterClassification[i], 1));
				break;
			}
		}
  }
	
  // Begin calculation
	double entropyClass = 0.0;
	for (int i = 0; i < classCount; i++)
	{
		entropyClass -= ((double)classesRatio[i] / (double)nodeClassClassification.size()) *
			std::log((double)classesRatio[i] / (double)nodeClassClassification.size());
	}

	double MI = 0.0;
	double clusterIclassJ = 0.0;
	for (int i = 0; i < clusters.size(); i++) // clusters
	{
		for (int j = 0; j < classCount; j++) // classes
		{
			clusterIclassJ = 0.0;
			for (int k = 0; k < nodeClassClassification.size(); k++)
			{
				if (nodeClassClassification[k] == classes[j])
					if (nodeClusterClassification[k] == clusters[i].first)
						clusterIclassJ++;
			}

			if (clusterIclassJ != 0.0)
			{
				clusterIclassJ = clusterIclassJ / (double)nodeClassClassification.size();
				MI += clusterIclassJ * std::log(clusterIclassJ /
				  (((double)classesRatio[j] / (double)nodeClassClassification.size()) * 
					((double)clusters[i].second / (double)nodeClassClassification.size())));
			}
		}
	}
	
	double entropyClusters = 0.0;
	for (int i = 0; i < clusters.size(); i++)
	{
		entropyClusters -= ((double)clusters[i].second / (double)nodeClassClassification.size()) *
			std::log((double)clusters[i].second / (double)nodeClassClassification.size());
	}

	double NMI = 0.0;
	if (entropyClass > entropyClusters)
		NMI = MI / entropyClass;
	else
		NMI = MI / entropyClusters;
}