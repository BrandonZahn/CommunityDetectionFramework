<h>CommunityDetectionFramework</h>

 - CiteseerCoraCorrelator, the correlation program and the data set of the Cora 
      and CiteSeer data sets.
      
 - CommunityDetectorLocalSearch, check lines 208 to 325 of the global.cpp file 
      for the modified code.
 
 - CommunityDetectorMemManage, check functions NSGA2::commcross and 
      NSGA2::advance for the modified code.
 
 - CommunityDetectorRefinedComms, check functions 
      individual::initialize, individual::comm_mutate and NSGA2::commcross for 
      the modified code.
 
 - GMLMakerDirected, a program to create a .gml file for analysis of directed 
      graphs using Gephi.
 
 - GMLMakerUndirected, a program to create a .gml file for analysis of 
      undirected graphs using Gephi.
 
 - KruskalTreeMaker, A program to find the modified kruskal tree from the 
      edges of a graph. Used to find a subset of edges with a minimised edge
      weighting.
      
 - NMICalculator, a program used to find the normalised mutual information of
      two input files. The classes file represents the predefined classes and
      the clusters input file represents the cluster tag ID found by the 
      community detector program.
      
 - TripAdvisorCorrelator, the correlation program to correlate the hotels with
      one another. It has a configuration file that includes some values and 
      file names of all the hotel file names.
 
