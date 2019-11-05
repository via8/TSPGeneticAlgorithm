#ifndef TSP_GENETIC_ALGORITHM_H
#define TSP_GENETIC_ALGORITHM_H

#include "TSPGenome.h"

class TSPGeneticAlgorithm final {
private:
	static const int MIN_GEN_SIZE = 4;	// minimum possible size of population to run the algorithm

	int const generations;		// number of generations will be generated until termination
	int const generationSize;	// number of genomes included in generation
	int const tournamentSize;	// number of genomes will be selected to compare during tournament
	int const crossoverRate;	// number of genomes will be selected to produce next generation
	double const mutationRate;	// value from [0, 1] interval - chance of random swap of 2 genes of single genome

	TSPGeneticAlgorithm(int generations, int generationSize, int crossoverRate, int tournamentSize, double mutationRate);
	TSPGenome** createNextGeneration(int citiesNumber, TSPGenome** currentGeneration = nullptr);
	TSPGenome** selection(TSPGenome** generation);
	void swapGenomes(TSPGenome** generation, int i, int j);
	void crossover(TSPGenome* parent1, TSPGenome* parent2, TSPGenome* child1, TSPGenome* child2, int genomeSize);
	int tournament(TSPGenome** generation, int generationSize);
	int findFittestGenome(TSPGenome const* const* generation);

public:
	// fabric method
	static TSPGeneticAlgorithm* create(int generations = 64, int generationSize = 128, int crossoverRate = 32, int tournamentSize = 3, double mutationRate = 0.1);
	int run(int citiesNumber, int const** pricesMatrix, int* sequence);
};

#endif /* TSP_GENETIC_ALGORITHM_H */