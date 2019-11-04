#include "TSPGenome.h"
#include <random>

#define TSP_GA_DEBUG

#ifdef TSP_GA_DEBUG
#include <iostream>
#define TSP_GA_DEBUG_REPORT(file, line) (std::cout << "ERROR: file " << file << ", line " << line << std::endl)
#endif /* TSP_DEBUG */

TSPGenome::TSPGenome(int citiesNumber, gene_t* genome, int* keyGenome) :
	genomeSize(citiesNumber - 1),
	genome(genome), 
	keyGenome(keyGenome),
	fitness(-1) /* not calculated */ {
	if (genome == nullptr) {
		this->genome = new gene_t[genomeSize];
		randomizeGenome();
	}
	if (keyGenome == nullptr) {
		this->keyGenome = new int[genomeSize];
		fillKeyGenome();
	}
}

int TSPGenome::getFitness() const {
#ifdef TSP_GA_DEBUG
	if (fitness == -1) {
		TSP_GA_DEBUG_REPORT(__FILE__, __LINE__);
		exit(-1);
	}
#endif /* TSP_GA_DEBUG */
	return fitness;
}

void TSPGenome::recalculateFitness(int const** pricesMatrix) {

#ifdef TSP_GA_DEBUG
	if (pricesMatrix == nullptr || genome == nullptr) {
		TSP_GA_DEBUG_REPORT(__FILE__, __LINE__);
		exit(-1);
	}
#endif /* TSP_GA_DEBUG */

	fitness = 0;
	// calculate whole path cost 
	for (int i = 0, currentCity = 0; i < genomeSize; ++i) {
		fitness += pricesMatrix[currentCity][genome[i]];
		currentCity = genome[i];
	}
	// including coming back from last city to the starting one
	fitness += pricesMatrix[genome[genomeSize - 1]][0];
}


void TSPGenome::randomizeGenome() {

#ifdef TSP_GA_DEBUG
	if (genome == nullptr) {
		TSP_GA_DEBUG_REPORT(__FILE__, __LINE__);
		exit(-1);
	}
#endif /* TSP_GA_DEBUG */

	//// to create random genome simply add all cities numbers to genome array
	//for (int i = 0; i < startingCity; ++i)
	//	genome[i] = i;
	//// except for starting city number
	//for (int i = startingCity + 1; i <= genomeSize; ++i)
	//	genome[i - 1] = i;

	// to create random genome simply add all cities numbers to genome array except for 0
	for (int i = 1; i <= genomeSize; ++i)
		genome[i - 1] = i;

	// and shuffle them using std::shuffle
	std::shuffle(genome, genome + genomeSize, std::default_random_engine(rand()));
}

void TSPGenome::fillKeyGenome() {

#ifdef TSP_GA_DEBUG
	if (genome == nullptr || keyGenome == nullptr) {
		TSP_GA_DEBUG_REPORT(__FILE__, __LINE__);
		exit(-1);
	}
#endif /* TSP_GA_DEBUG */

	// thank to this array of keys index of k city in genome
	// will be found by the addres keyGenome[genome[k] - 1]
	for (int i = 0; i < genomeSize; ++i)
		keyGenome[genome[i] - 1] = i;
}

#ifdef TSP_GA_DEBUG
#undef TSP_GA_DEBUG
#endif /* TSP_GA_DEBUG */
