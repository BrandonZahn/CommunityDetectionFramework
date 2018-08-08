#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <iterator>
#include <algorithm>
#include <vector>
#include <functional>
#include <cmath>
#include <numeric>
#include <thread>

#include "graph.h"

using namespace std;

// Computes the distance between two std::vectors
template <typename T>
double vectors_distance(const vector<T>& a, const vector<T>& b)
{
	vector<double> auxiliary;

	transform(a.begin(), a.end(), b.begin(), back_inserter(auxiliary),
		[](T element1, T element2) { return pow((element1 - element2), 2); });

	return sqrt(accumulate(auxiliary.begin(), auxiliary.end(), 0.0));
}
// end template vectors_distance

template <typename T>
int sgn(T val) { return (T(0) < val) - (val < T(0)); }

double findDWCosine(vector<float> &A, vector<float> &B)
{
	float WHD = 0.0, WHDNumA = 0.0, WHDNumB = 0.0, WHDDenA = 0.0, 
	WHDDenB = 0.0, WDC = 0.0, WDCNum = 0.0, WDCDenA = 0.0, WDCDenB = 0.0;
	for (size_t it = 0; it < A.size(); it++)
	{
		WHDNumA += A[it] * (1 - sgn(A[it] * B[it]));
		WHDNumB += B[it] * (1 - sgn(B[it] * A[it]));
		WHDDenA += A[it];
		WHDDenB += B[it];
		WDCNum += A[it] * B[it];
		WDCDenA += A[it] * A[it];
		WDCDenB += B[it] * B[it];
	}
	WDCDenA = sqrt(WDCDenA);
	WDCDenB = sqrt(WDCDenB);
	WHD = (WHDNumA / WHDDenA) + (WHDNumB / WHDDenB);
	WDC = WDCNum / (WDCDenA * WDCDenB);
	return WDC / ((WHD * WHD) + WDC);
}

int main(void)
{
	struct Hotel
	{
		int id;
		string location = "unknown";
		int price = -1;
		vector<double>* ratings;
		vector<vector<float>>* aspects;
		int reviewCount;
	};

	// Get the names of the hotel files using a config file
	cout << "Getting configuration file config.txt..." << endl;
	
	// Use this hotelCount variable to change the number of nodes
	int hotelCount;
	double threshold;
	double featuredAspectWeighting;
	double differenceThreshold;

	ifstream inputFile1("config.txt");
	string inputData;

	cout << "Configuration file received!" << endl;
	cout << "Run with settings:" << endl;

	// Get the number of hotels to analyze (number of vertices)
	getline(inputFile1, inputData);
	hotelCount = atoi(inputData.c_str());
	cout << "Hotel count: " << hotelCount << endl;

	// Get the threshold value from the config file. It is used to 
	// set a starting point for the DFS later on
	getline(inputFile1, inputData);
	threshold = atof(inputData.c_str());
	cout << "Threshold: " << threshold << endl;

	// Get the weighting for featured words which will be used in
	// the weighted matrix calculations
	getline(inputFile1, inputData);
	featuredAspectWeighting = atof(inputData.c_str());
	cout << "Featured aspect weighting: " << featuredAspectWeighting << endl;

	getline(inputFile1, inputData);
	differenceThreshold = atof(inputData.c_str());
	cout << "Difference threshold: " << differenceThreshold << endl;

	// Get the name of the featured words list file
	getline(inputFile1, inputData);
	ifstream featureWordsFile(inputData);
	cout << "Getting featured aspects list file from " << inputData << endl;

	vector<vector<string>>* featureWordsArray = new vector<vector<string>>();
	cout << "Building the array of feature words" << endl;

	for (int i = 0; i < 7; i++)
	{
		getline(featureWordsFile, inputData);
		string delimiter = "\t";
		size_t position = 0;
		string token;

		vector<string> tempVector;

		while ((position = inputData.find(delimiter)) != inputData.npos)
		{
			token = inputData.substr(0, position);
			tempVector.push_back(token);
			inputData.erase(0, position + delimiter.length());
		}
		tempVector.push_back(inputData);

		featureWordsArray->push_back(tempVector);
	}

	cout << "Feature words array built!" << endl;

	// Initialize all the arrays used for calculations
	// An i x 8 array of doubles to store all the rating data. ith row is the ith hotel
	vector<Hotel>* hotels = new vector<Hotel>(hotelCount);

	// Get the filename that stores the names of the hotel files
	getline(inputFile1, inputData);
	cout << "Getting hotel file names from " << inputData << endl;
	cout << "Processed hotel " << inputData << endl;
	ifstream inputFile2(inputData); 
	// inputData should now store the name of the hotel listings file
	cout << "Hotel file names file received!" << endl;

	vector<string>* hotelFileNames = new vector<string>;

	cout << "Processing hotel file names..." << endl;

	// -------------------------------------------------------------
	// Process the list of hotel file names file
	// -------------------------------------------------------------
	for (int i = 0; i < 7; i++)
	{
		getline(inputFile2, inputData);
	}
	string delimiter = "hotel";
	size_t namePosition = 39;
	string token;

	for (int i = 0; i < hotelCount; i++)
	{
		getline(inputFile2, inputData);
		token = inputData.substr(namePosition, inputData.length());
		hotelFileNames->push_back(token);

		// Get the hotel id from the name of the file
		delimiter = "_";
		token.erase(0, 6);
		int idPosition = token.find(delimiter);
		token = token.substr(0, idPosition);

		hotels->at(i).id = atoi(token.c_str());
		cout << "Processed hotel " << i << endl;
	}

	cout << "Hotel file names processed!" << endl;

	// -------------------------------------------------------------
	// Start reading each hotel file
	// -------------------------------------------------------------
	for (int i = 0; i < hotelCount; i++)
	{
		vector<vector<float> >* aspectsArray =
		new vector<vector<float> >(8, vector<float>());

		vector<double>* ratingsArray = new vector<double>(8);
		double ratingsCount[8] = {0,0,0,0,0,0,0,0};

		string inputFile;
		inputFile = hotelFileNames->at(i);
		ifstream inputFileI(inputFile);

		cout << "File " << (i + 1) << " of " << hotelCount << ": ";
		cout << "Processing hotel file "<< inputFile << "... ";

		int reviewCounter = 0;
		while (!inputFileI.eof())
		{
			getline(inputFileI, inputData); 

			if (inputData.substr(0, 8) == "<Rating>")
			{
				inputData.erase(0, 8);
				delimiter = "\t";

				int k = 0;
				while ((namePosition = inputData.find(delimiter)) != inputData.npos)
				{
					double tempDouble = ratingsArray->at(k);
					token = inputData.substr(0, namePosition);

					if (token != "-1")
					{
						tempDouble += atof(token.c_str());
						ratingsCount[k]++;
					}

					ratingsArray->at(k) = tempDouble;
					inputData.erase(0, namePosition + delimiter.length());
					k++;
				}
			}
			else if (inputData == "<Aspects>")
			{
				delimiter = "\t";
				for (size_t j = 0; j < 8; j++)
				{
					getline(inputFileI, inputData);

					if (inputData.substr(0, 1) != "0\t")
					{
						namePosition = inputData.find(delimiter);
						inputData.erase(0, namePosition + delimiter.length());

						while ((namePosition = inputData.find(delimiter)) != inputData.npos)
						{
							int id;
							token = inputData.substr(0, namePosition);

							string delimiterTwo = "(";
							size_t wordPosition = token.find(delimiterTwo);
							string tempString = token.substr(0, wordPosition);
							id = atoi(tempString.c_str());
							token.erase(0, wordPosition + delimiterTwo.length());

							string delimiterThree = ")";
							size_t wordCountPosition = token.find(delimiterThree);

							string word = token.substr(0, wordCountPosition);
							token.erase(0, (wordCountPosition + delimiterThree.length() + 1));

							float count = stof(token);

							if (id + 1 > aspectsArray->at(j).size())
							{
								int difference = id + 1 - aspectsArray->at(j).size();
								for (int it = 0; it < difference; it++) { aspectsArray->at(j).push_back(0.0); }
							}

							if (j != 0)
							{
								size_t vectorSize = featureWordsArray->at(j-1).size();
								for (size_t n = 0; n < featureWordsArray->at(j-1).size(); n++)
								{
									if (featureWordsArray->at(j-1).at(n) == word)
									{
										count *= featuredAspectWeighting;
										break;
									}
								}
							}
							aspectsArray->at(j).at(id) += count;
							inputData.erase(0, namePosition + delimiter.length());
						}
					}
				}
			}
			reviewCounter++;
		}

		// Divide by the count values to get an average
		for (int j = 0; j < 8; j++)
		{
			double tempDivider = ratingsArray->at(j);
			tempDivider /= ratingsCount[j];
			ratingsArray->at(j) = tempDivider;
		}
		hotels->at(i).ratings = ratingsArray;
		hotels->at(i).reviewCount = reviewCounter;	
		hotels->at(i).aspects = aspectsArray;
		cout << "Review count: " << reviewCounter << endl;
	}

	cout << "All "<< hotelCount << " hotel files processed!" << endl;

	// Get the hotel's location and price
	getline(inputFile1, inputData);
	ifstream inputFile3(inputData);

	cout << "Getting price location listings file " << inputData << "..." << endl;
	getline(inputFile3, inputData);

	int hotelId;
	int price;
	string location;

	while (!inputFile3.eof())
	{
		getline(inputFile3, inputData);
		inputData.erase(0, 6);

		delimiter = "_";
		token = inputData.substr(0, inputData.find(delimiter));
		hotelId = atoi(token.c_str());
		delimiter = "\t";
		inputData.erase(0, inputData.find(delimiter) + delimiter.length());

		for (int i = 0; i < hotelCount; i++)
		{
			if (hotels->at(i).id == hotelId)
			{
				token = inputData.substr(0, inputData.find(delimiter));
				if (token == "unknown")
					hotels->at(i).price = -1;
				else
					hotels->at(i).price = atoi(token.c_str());

				inputData.erase(0, inputData.find(delimiter) + delimiter.length());
				hotels->at(i).location = inputData;
			}
		}
	}
	cout << "All " << hotelCount << " hotel files successfully updated with location and prices!" << endl;

	ofstream names("labels.txt");
	for (size_t i = 0; i < hotelCount - 1; i++)
	{
		names << hotels->at(i).id << " " << hotels->at(i).location << " " << hotels->at(i).price << " ";
	}
	names << hotels->at(hotelCount - 1).id << " " << hotels->at(hotelCount - 1).location << " " << hotels->at(hotelCount - 1).price;
	names.close();

	cout << "Building the adjacency and weighted matrices for the input.dat file..." << endl;

	vector<vector<double> >* adjacencyMatrix =
		new vector<vector<double> >(hotelCount, vector<double>(hotelCount));

	vector<vector<double> >* weightedMatrix =
		new vector<vector<double> >(hotelCount, vector<double>(hotelCount));

	double correlation;
	bool connected = false;

	ofstream euclideanDistanceFile("EUCDIS.txt");
	for (int i = 0; i < hotelCount; i++)
	{
		for (int j = 0; j < hotelCount; j++)
		{
			if (i != j)
			{
				vector<double>* hotelIRating = hotels->at(i).ratings;
				vector<double>* hotelJRating = hotels->at(j).ratings;

				correlation = vectors_distance(*hotelIRating, *hotelJRating);
				euclideanDistanceFile << correlation << " ";
			}
			else
				euclideanDistanceFile << "0 ";
		}
		euclideanDistanceFile << endl;
	}
	euclideanDistanceFile.close();

	while (!connected)
	{
		for (int i = 0; i < hotelCount; i++)
		{
			for (int j = 0; j < hotelCount; j++)
			{
				if (i != j)
				{
					vector<double>* hotelIRating = hotels->at(i).ratings;
					vector<double>* hotelJRating = hotels->at(j).ratings;

					correlation = vectors_distance(*hotelIRating, *hotelJRating);
					if (correlation <= threshold)
						adjacencyMatrix->at(i).at(j) = 1;
					else
						adjacencyMatrix->at(i).at(j) = 0;
				}
			}
		}
		int numEdges = 0;
		Graph graf(hotelCount);
		for (int i = 0; i < hotelCount; i++)
		{
			for (int j = 0; j < hotelCount; j++)
			{
				if (i < j && adjacencyMatrix->at(i).at(j) == 1)
				{
					graf.addEdge(i, j);
					numEdges++;
				}
			}
		}
		cout << "Edge count: " << numEdges << endl;
		graf.DFS(0);
		cout << "DFS visited " << graf.GetVisitedCount() << " nodes!" << endl;
		connected = graf.checkConnected();
		if (!connected)
		{
			threshold += 0.005;
			cout << "Threshold set to: " << threshold << endl;
		}
		else
			cout << "Connected! at threshold: " << threshold << endl;
	}

	for (int i = 0; i < hotelCount; i++)
	{
		for (int j = 0; j < hotelCount; j++)
		{
			if (i != j && i < j)
			{
				vector<vector<float> >* hotelIAspects = hotels->at(i).aspects;
				vector<vector<float> >* hotelJAspects = hotels->at(j).aspects;
				double dwCosine = 0.0;
				vector<float> hotelACounts, hotelBCounts;

				for (int k = 0; k < 8; k++)
				{
					if (hotelIAspects->at(k).size() != hotelJAspects->at(k).size())
					{
						if (hotelIAspects->at(k).size() < hotelJAspects->at(k).size())
						{
							int difference = hotelJAspects->at(k).size() - hotelIAspects->at(k).size();
							for (int it = 0; it < difference; it++) { hotelIAspects->at(k).push_back(0.0); }
						}
						else if (hotelJAspects->at(k).size() < hotelIAspects->at(k).size())
						{
							int difference = hotelIAspects->at(k).size() - hotelJAspects->at(k).size();
							for (int it = 0; it < difference; it++) { hotelJAspects->at(k).push_back(0.0); }
						}
					}
					for (int it = 0; it < hotelJAspects->at(k).size(); it++)
					{
						hotelACounts.push_back(hotelIAspects->at(k).at(it) / hotels->at(i).reviewCount);
						hotelBCounts.push_back(hotelJAspects->at(k).at(it) / hotels->at(j).reviewCount);
					}
					float result = findDWCosine(hotelACounts, hotelBCounts);
					if (result > dwCosine)
						dwCosine = result;
				}

				if (dwCosine <= 0.25)
				{
					weightedMatrix->at(i).at(j) = 0;
					weightedMatrix->at(j).at(i) = 0;
				}
				else
				{
					weightedMatrix->at(i).at(j) = dwCosine - 0.25;
					weightedMatrix->at(j).at(i) = dwCosine - 0.25;
				}
				cout << "Element [" << i << "][" << j << "] complete!..." << endl;
			}
		}
	}

	cout << "Matrices built!" << endl;

	cout << "Creating output file input.dat..." << endl;
	// Output the input.dat file for the NSGAII calculations
	ofstream outputFile("input.dat");

	outputFile << "<INSERT RANDOM NUMER FOR SEED GENERATION>" << endl;
	outputFile << hotelCount << endl;
	for (int i = 0; i < hotelCount; i++)
	{
		for (int j = 0; j < (hotelCount - 1); j++)
		{
			if (i != j)
				outputFile << adjacencyMatrix->at(i).at(j) << " ";
			else
				outputFile << "0 ";
		}
		
		if (i != (hotelCount - 1))
			outputFile << adjacencyMatrix->at(i).at(hotelCount - 1) << endl;
		else
			outputFile << "0" << endl;
	}

	outputFile << endl;

	for (int i = 0; i < hotelCount; i++)
	{
		for (int j = 0; j < (hotelCount - 1); j++)
		{
			if (i != j)
				outputFile << weightedMatrix->at(i).at(j) << " ";
			else
				outputFile << "0 ";
		}
		
		if (i != (hotelCount - 1))
			outputFile << weightedMatrix->at(i).at(hotelCount - 1) << endl;
		else
			outputFile << "0" << endl;
	}

	outputFile << endl;

	outputFile << "<INSERT POPULATION SIZE (MULTIPLE OF 4)>" << endl;
	outputFile << "<INSERT NUMBER OF GENERATIONS>" << endl;
	outputFile << "<INSERT NUMBER OF OBJECTIVES>" << endl;
	outputFile << "<INSERT NUMBER OF CONSTRAINTS>" << endl;
	outputFile << "<pmut_comm>" << endl;
	outputFile << "<pcross_comm>" << endl;
	outputFile << "<eta_c>" << endl;
	outputFile << "<eta_m>" << endl;
}