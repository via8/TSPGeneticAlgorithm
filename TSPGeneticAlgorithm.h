#ifndef TSP_GENETIC_ALGORITHM_H
#define TSP_GENETIC_ALGORITHM_H

#include "TSPGenome.h"

class TSPGeneticAlgorithm final {
private:
	static const int MIN_GEN_SIZE = 4;	// minimum possible size of population to run the algorithm

	int generationSize;			// number of genomes included in generation
	int tournamentSize;			// number of genomes will be selected to compare during tournament
	int crossoverRate;			// number of genomes will be selected to produce next generation
	double mutationRate;		// value from [0, 1] interval - chance of random swap of 2 genes of single genome

	TSPGeneticAlgorithm(int generationSize, int crossoverRate, int tournamentSize, double mutationRate);
	TSPGenome** Selection(TSPGenome** generation);
	int Tournament(TSPGenome** generation, int generationSize);
	void SwapGenomes(TSPGenome** generation, int i, int j);
	void Crossover(int const* keyParent1, gene_t const* parent2, gene_t* child1, gene_t* child2, int genomeSize);
	void Mutate(gene_t* genome, int genomeSize);

public:
	// fabric method
	static TSPGeneticAlgorithm* createTSPGeneticAlgorithm(int generationSize, int crossoverRate, int tournamentSize, double mutationRate);
	TSPGenome** createNextGeneration(int citiesNumber, TSPGenome** currentGeneration = nullptr);
};

#endif /* TSP_GENETIC_ALGORITHM_H */