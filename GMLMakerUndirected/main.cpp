#include <iostream>
#include <fstream>
#include <vector>
#include <string>

int main()
{
  std::string line;
  std::ifstream fstream("inputTA.dat");

  int nodes;
  std::vector<std::vector<int>> edges;
  std::vector<std::vector<double>> values;
  std::vector<std::string> labels;
  std::vector<int> comClass;
  
  std::getline(fstream, line);
  nodes = std::stoi(line);
  std::getline(fstream, line);

  std::getline(fstream, line);
  std::string delimiter = " ";
  std::string token;
  size_t position = 0;

  while ((position = line.find(delimiter)) != line.npos)
  {
    token = line.substr(0, position);
    labels.push_back(token);
    line.erase(0, position + delimiter.length());
  }
  labels.push_back(line);

  std::getline(fstream, line);

  for (int i = 0; i < nodes; i++)
  {
    std::getline(fstream, line);
    std::vector<int> row;
    while ((position = line.find(delimiter)) != line.npos)
    {
      token = line.substr(0, position);
      row.push_back(stoi(token));
      line.erase(0, position + delimiter.length());
    }
    row.push_back(stoi(line));
    edges.push_back(row);
  }

  std::getline(fstream, line);

  for (int i = 0; i < nodes; i++)
  {
    std::getline(fstream, line);
    std::vector<double> row;
    while ((position = line.find(delimiter)) != line.npos)
    {
      token = line.substr(0, position);
      row.push_back(stod(token));
      line.erase(0, position + delimiter.length());
    }
    row.push_back(stod(line));
    values.push_back(row);
  }

  std::getline(fstream, line);

  std::getline(fstream, line);
  while ((position = line.find(delimiter)) != line.npos)
  {
    token = line.substr(0, position);
    comClass.push_back(stoi(token));
    line.erase(0, position + delimiter.length());
  }
  comClass.push_back(stoi(line));

  std::ofstream gephiFile("TA.gml");
  gephiFile << "graph\n[\nCreator \"Brandon Zahn\"\ndirected 0\n";

  for (int i = 0; i < nodes; i++)
  {
    gephiFile << "node\n[\n";
    gephiFile << "id " << i << std::endl;
    gephiFile << "label \"" << "Name: " << labels[3 * i] << ", Loc: " << labels[3 * i + 1] 
		<< ", Price: " << labels[3 * i + 2] << "\"" << std::endl;
    gephiFile << "ModularityClass \"" << comClass[i] << "\"" << std::endl;
    gephiFile << "]\n";
  }

  int edgeId = 0;
  for (int i = 0; i < nodes; i++)
  {
    for (int j = 0; j < nodes; j++)
    {
      if (edges[i][j] == 1 && i < j)
      {
        gephiFile << "edge\n[\n";
        gephiFile << "id " << edgeId << std::endl;
        gephiFile << "source " << i << std::endl;
        gephiFile << "target " << j << std::endl;
        gephiFile << "value " << values[i][j] << std::endl;
        gephiFile << "]\n";

        edgeId++;
      }
    }
  }

  gephiFile << "]";

  return 0;
}