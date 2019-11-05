#include "TSPGenome.h"
#include <random>

TSPGenome::TSPGenome(int citiesNumber) :
	genomeSize(citiesNumber - 1),
	fitness(-1) /* not calculated */ {

	this->genome  = new int[genomeSize];
	this->indexes = new int[citiesNumber];	// 0 element of this array is excess
}

void TSPGenome::randomizeGenome() {
	// to create random genome simply add all cities numbers to genome array except for 0
	for (int i = 0; i < genomeSize; ++i)
		genome[i] = i + 1;

	// and shuffle them using std::shuffle
	std::shuffle(genome, genome + genomeSize, std::default_random_engine(std::rand()));

	// thanks to this array of keys, index of k city in genome
	// will be found by the addres indexes[genome[k]]
	for (int i = 0; i < genomeSize; ++i)
		indexes[genome[i]] = i;
}

void TSPGenome::recalculateFitness(int const** pricesMatrix) {
	fitness = 0;
	// calculate whole path cost 
	for (int i = 0, currentCity = 0; i < genomeSize; ++i) {
		fitness += pricesMatrix[currentCity][genome[i]];
		currentCity = genome[i];
	}
	// including coming back from last city to the starting one
	fitness += pricesMatrix[genome[genomeSize - 1]][0];
}

void TSPGenome::mutate(double mutationRate) {
	// generate random double number from [0, 1] interval
	double random = (double)std::rand() / ((double)RAND_MAX + 1.0);

	// if this value exceeds mutation rate
	if (random > mutationRate) {
		int  firstIndex = std::rand() % genomeSize;
		int secondIndex = std::rand() % genomeSize;

		// if indexes matched then simply assign second to 0 gene if first isn't 0, 1 otherwise
		// (to avoid long search of another random for small genomes)
		if (firstIndex == secondIndex) {
			secondIndex = (firstIndex == 0) ? 1 : 0;
		}

		// swap genes of first and second random indexes
		// swap their keys in indexes array also
		int tempGene  = genome[firstIndex];
		int tempIndex = indexes[genome[firstIndex]];

		genome[firstIndex] = genome[secondIndex];
		indexes[genome[firstIndex]] = indexes[genome[secondIndex]];

		genome[secondIndex] = tempGene;
		indexes[genome[secondIndex]] = tempIndex;
	}
}