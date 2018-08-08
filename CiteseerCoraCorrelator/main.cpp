#include <vector>
#include <numeric>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <algorithm>

using namespace std;

int main(void)
{
	// Structures used for the papers
	struct CoraPaper
	{
		int paperId;
		int vocabularyArray[1433];
		string paperClass;
		vector<int>* cites;
		vector<int>* citedBy = new vector<int>();
	};

	struct CiteseerPaper
	{
		string paperId;
		int vocabularyArray[3703];
		string paperClass;
		vector<string>* cites;
		vector<string>* citedBy = new vector<string>();
	};

	vector<CoraPaper>* paperList = new vector<CoraPaper>();

	vector<vector<int>>* citesMatrix = new vector<vector<int>>();
	vector<vector<double>>* weightedMatrix = new vector<vector<double>>();
	string inputData;
	ifstream coraCites("cora.cites");

	string delimiter = "\t";
	string token;
	size_t position = 0;
	
	cout << "Building the list of papers...\n";

	vector<pair<int, int>>* citingsList = new vector<pair<int, int>>();

	while (!coraCites.eof())
	{
		getline(coraCites, inputData);
		position = inputData.find(delimiter);
		int nameOne = stoi(inputData.substr(0, position));
		inputData.erase(0, position + 1);
		int nameTwo = stoi(inputData);
		
		citingsList->push_back(make_pair(nameOne, nameTwo));
	}

	cout << "Paper list complete!\n";

	ifstream coraContent("cora.content");

	while (!coraContent.eof())
	{
		getline(coraContent, inputData);
		position = inputData.find(delimiter);
		int paperId = stoi(inputData.substr(0, position));
		CoraPaper tempPaper;
		tempPaper.paperId = paperId;
		inputData.erase(0, position + 1);

		for (int i = 0; i < 1433; i++)
		{
			string vocab = inputData.substr(0, 1);
			tempPaper.vocabularyArray[i] = stoi(vocab);
			inputData.erase(0, 2);
		}
		tempPaper.paperClass = inputData;

		paperList->push_back(tempPaper);
	}

	for (int i = 0; i < citingsList->size(); i++)
	{
		for (int j = 0; j < paperList->size(); j++)
		{
			if (paperList->at(j).paperId == citingsList->at(i).first)
			{
				paperList->at(j).citedBy->push_back(citingsList->at(i).second);
				break;
			}
		}
	}

	cout << "Building adjacency array....\n";
	for (int i = 0; i < 2708; i++)
	{
		vector<int> intVector;
		vector<double> doubleVector;
		citesMatrix->push_back(intVector);
		weightedMatrix->push_back(doubleVector);

		for (int j = 0; j < 2708; j++)
		{
			citesMatrix->at(i).push_back(0);
			weightedMatrix->at(i).push_back(0);
		}

		CoraPaper* paper = &paperList->at(i);
		size_t vectorSize = paper->citedBy->size();
		for (size_t j = 0; j < vectorSize; j++)
		{
			for (int k = 0; k < 2708; k++)
			{
				if (paper->citedBy->at(j) == paperList->at(k).paperId)
				{
					citesMatrix->at(i).at(k) = 1;
					break;
				}
			}
		}
	}
	
	// Create a text file containing all the names of the papers in
	// the Cora data set 
	ofstream names("labelsCora.txt");
	for (size_t i = 0; i < 2708 - 1; i++)
	{
		names << paperList->at(i).paperId << " ";
	}
	names << paperList->at(2707).paperId;

	names.close();

	// Create a text file containing all the classes of the papers in
	// the Cora data set 
	ofstream classes("classesCora.txt");
	for (size_t i = 0; i < 2708 - 1; i++)
	{
		classes << paperList->at(i).paperClass << endl;
	}
	classes << paperList->at(2707).paperClass;

	classes.close();

	cout << "Adjacency array built!\nBuilding weighted matrix...\n";
	for (int i = 0; i < 2708; i++)
	{
		for (int j = 0; j < 2708; j++)
		{
			for (int k = 0; k < 1433; k++)
			{
				if (paperList->at(i).vocabularyArray[k] == 1 && 
					paperList->at(j).vocabularyArray[k] == 1)
					weightedMatrix->at(i).at(j) += 1;
			}
		}
	}
	
	/*cout << "Refining weighted matrix...\n";
	for (int i = 0; i < 2708; i++)
	{
		for (int j = 0; j < 2708; j++)
		{
			if (weightedMatrix->at(i).at(j) > 5)
			{
				weightedMatrix->at(i).at(j) -= 5;
				weightedMatrix->at(i).at(j) *= 1.25;
			}
			else 
				weightedMatrix->at(i).at(j) = 0.0;
		}
	}*/

	cout << "Creating output file 'inputCora.dat'...\n";
	
	// Output the inputCora.dat file for the NSGAII algorithm
	ofstream outputFile("inputCora.dat");

	outputFile << "<INSERT RANDOM NUMER FOR SEED GENERATION>\n";
	outputFile << "2708" << endl;
	for (int i = 0; i < 2708; i++)
	{
		for (int j = 0; j < (2708 - 1); j++)
		{
			if (i != j)
				outputFile << citesMatrix->at(j).at(i) << " ";
			else
				outputFile << "0 ";
		}
		if (i != (2708 - 1))
			outputFile << citesMatrix->at(2708 - 1).at(i) << endl;
		else
			outputFile << "0\n";
	}

	outputFile << endl;

	for (int i = 0; i < 2708; i++)
	{
		for (int j = 0; j < (2708 - 1); j++)
		{
			if (i != j)
				outputFile << weightedMatrix->at(i).at(j) << " ";
			else
				outputFile << "0 ";
		}
		if (i != (2708 - 1))
			outputFile << weightedMatrix->at(i).at(2708 - 1) << endl;
		else
			outputFile << "0\n";
	}

	outputFile << endl;	
	
	/*vector<float> histogram;
	int maxValue = 0;
	for (int i = 0; i < 2708; i++)
	{
		for (int j = 0; j < 2708; j++)
		{
			if (weightedMatrix->at(i).at(j) > maxValue)
				maxValue = weightedMatrix->at(i).at(j);
		}
	}
	outputFile << "MAX IS " << maxValue << endl;
	
	for (int i = 0; i < maxValue; i ++)	{ histogram.push_back(0.0); }
	
	for (int i = 0; i < 2708; i++)
	{
		for (int j = 0; j < 2708; j++)
		{
			if (i != j && i < j)
			{
				histogram[(int)weightedMatrix->at(i).at(j)]++;
			}
		}
	}
	
	for (int i = 0; i < histogram.size(); i++)
	{
		outputFile << histogram[i] << "\t";
	}
	
	outputFile << endl;*/
	outputFile << "<INSERT POPULATION SIZE (MULTIPLE OF 4)>\n";
	outputFile << "<INSERT NUMBER OF GENERATIONS>\n";
	outputFile << "<INSERT NUMBER OF OBJECTIVES>\n";
	outputFile << "<INSERT NUMBER OF CONSTRAINTS>\n";
	outputFile << "<pmut_comm>\n";
	outputFile << "<pcross_comm>\n";
	outputFile << "<eta_c>\n";
	outputFile << "<eta_m>\n";

	coraCites.close();
	coraContent.close();
	outputFile.close();

	delete citesMatrix;
	delete weightedMatrix;
	delete paperList;
	delete citingsList;

	// ---------------------------------------------------
	// Start the Citeseer data set
	// ---------------------------------------------------
	ifstream citeseerContent("citeseer.content");
	
	vector<CiteseerPaper>* citeseerPaperList = new vector<CiteseerPaper>();

	while (!citeseerContent.eof())
	{
		getline(citeseerContent, inputData);
		position = inputData.find(delimiter);
		string paperId = inputData.substr(0, position);
		
		CiteseerPaper newPaper;
		newPaper.paperId = paperId;
		newPaper.cites = new vector<string>();
		newPaper.citedBy = new vector<string>();

		inputData.erase(0, position + 1);
		
		for (int i = 0; i < 3703; i++)
		{
			string vocab = inputData.substr(0, 1);
			newPaper.vocabularyArray[i] = stoi(vocab);
			inputData.erase(0, 2);
		}
		newPaper.paperClass = inputData;
		citeseerPaperList->push_back(newPaper);
	}

	vector<vector<int>>* citeseerCitesMatrix = new vector<vector<int>>();
	vector<vector<double>>* citeseerWeightedMatrix = new vector<vector<double>>();
	ifstream citeseerCites("citeseer.cites");

	delimiter = "\t";
	token;
	position = 0;

	cout << "Building the citings for Citeseer papers...\n";

	vector<pair<string, string>>* citeseerCitingsList = 
		new vector<pair<string, string>>();
	while (!citeseerCites.eof())
	{
		getline(citeseerCites, inputData);
		position = inputData.find(delimiter);
		string nameOne = inputData.substr(0, position);
		inputData.erase(0, position + 1);
		string nameTwo = inputData;
		
		citeseerCitingsList->push_back(make_pair(nameOne, nameTwo));
	}

	for (int i = 0; i < citeseerCitingsList->size(); i++)
	{
		for (int j = 0; j < citeseerPaperList->size(); j++)
		{
			if (citeseerPaperList->at(j).paperId == 
				citeseerCitingsList->at(i).first)
			{
				citeseerPaperList->at(j).citedBy->
					push_back(citeseerCitingsList->at(i).second);
				break;
			}
		}
	}

	cout << "Paper list complete!\nBuilding adjacency array...\n";
	for (int i = 0; i < 3312; i++)
	{
		vector<int> intVector;
		vector<double> doubleVector;
		citeseerCitesMatrix->push_back(intVector);
		citeseerWeightedMatrix->push_back(doubleVector);

		for (int j = 0; j < 3312; j++)
		{
			citeseerCitesMatrix->at(i).push_back(0);
			citeseerWeightedMatrix->at(i).push_back(0);
		}

		CiteseerPaper* paper = &citeseerPaperList->at(i);
		size_t vectorSize = paper->citedBy->size();
		for (size_t j = 0; j < vectorSize; j++)
		{
			for (int k = 0; k < 3312; k++)
			{
				if (paper->citedBy->at(j) == 
					citeseerPaperList->at(k).paperId)
				{
					citeseerCitesMatrix->at(i).at(k) = 1;
					break;
				}
			}
		}
	}

	cout << "Adjacency array built!\nBuilding weighted matrix...\n";
	for (int i = 0; i < 3312; i++)
	{
		for (int j = 0; j < 3312; j++)
		{
			for (int k = 0; k < 3703; k++)
			{
				if (citeseerPaperList->at(i).vocabularyArray[k] == 1 && 
					citeseerPaperList->at(j).vocabularyArray[k] == 1)
					citeseerWeightedMatrix->at(i).at(j) += 1;
			}
		}
	}

	cout << "Creating output file 'input.dat'...\n";
	
	// Output the inputCiteseer.dat file for the NSGAII calculations
	ofstream outputFile2("intputCiteseer.dat");

	outputFile2 << "<INSERT RANDOM NUMER FOR SEED GENERATION>\n";
	outputFile2 << "3312" << endl;
	for (int i = 0; i < 3312; i++)
	{
		for (int j = 0; j < (3312 - 1); j++)
		{
			if (i != j)
				outputFile2 << citeseerCitesMatrix->at(j).at(i) << " ";
			else
				outputFile2 << "0 ";
		}
		if (i != (3312 - 1))
			outputFile2 << citeseerCitesMatrix->at(3312 - 1).at(i) << endl;
		else
			outputFile2 << "0\n";
	}

	outputFile2 << endl;

	for (int i = 0; i < 3312; i++)
	{
		for (int j = 0; j < (3312 - 1); j++)
		{
			if (i != j)
				outputFile2 << citeseerWeightedMatrix->at(i).at(j) << " ";
			else
				outputFile2 << "0 ";
		}
		if (i != (3312 - 1))
			outputFile2 << citeseerWeightedMatrix->at(i).at(3312 - 1) << endl;
		else
			outputFile2 << "0\n";
	}

	outputFile2 << endl;

	outputFile2 << "<INSERT POPULATION SIZE (MULTIPLE OF 4)>\n";
	outputFile2 << "<INSERT NUMBER OF GENERATIONS>\n";
	outputFile2 << "<INSERT NUMBER OF OBJECTIVES>\n";
	outputFile2 << "<INSERT NUMBER OF CONSTRAINTS>\n";
	outputFile2 << "<pmut_comm>\n";
	outputFile2 << "<pcross_comm>\n";
	outputFile2 << "<eta_c>\n";
	outputFile2 << "<eta_m>\n";

	// Create a text file containing all the names of the papers in
	// the Citeseer data set 
	ofstream namesCite("labelsCite.txt");
	for (size_t i = 0; i < 3312 - 1; i++)
	{
		namesCite << citeseerPaperList->at(i).paperId << " ";
	}
	namesCite << citeseerPaperList->at(3311).paperId;

	namesCite.close();

	// Create a text file containing all the classes of the papers in
	// the Citeseer data set 
	ofstream classescite("classesCite.txt");
	for (size_t i = 0; i < 3312 - 1; i++)
	{
		classescite << citeseerPaperList->at(i).paperClass << endl;
	}
	classescite << citeseerPaperList->at(3311).paperClass;

	classescite.close();	
	outputFile2.close();
}