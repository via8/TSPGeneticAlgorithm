#include "TSPGeneticAlgorithm.h"
#include <cstdlib>

#define TSP_GA_DEBUG

#ifdef TSP_GA_DEBUG
#include <iostream>
#include <algorithm>
#define TSP_GA_DEBUG_REPORT(file, line) (std::cout << "ERROR: file " << file << ", line " << line << std::endl)
#endif /* TSP_GA_DEBUG */

TSPGeneticAlgorithm::TSPGeneticAlgorithm(int generationSize, int crossoverRate, int tournamentSize, double mutationRate) :
	generationSize(generationSize),
	tournamentSize(tournamentSize),
	crossoverRate(crossoverRate),
	mutationRate(mutationRate) {
}

TSPGeneticAlgorithm* TSPGeneticAlgorithm::createTSPGeneticAlgorithm(int generationSize, int crossoverRate, int tournamentSize, double mutationRate) {
	// make sure that generation size is even as crossover operation always builds pairs of children
	if (generationSize % 2 != 0)
		generationSize++;
	
	if (generationSize - crossoverRate < tournamentSize - 1 ||	// otherwise it would be impossible to arrange tournament
		generationSize < MIN_GEN_SIZE			 ||	// minimum possible size of generation to run the algorithm
		generationSize < crossoverRate			 ||	// impossible to cross more genomes then we have in generation
		tournamentSize < 0 ||	// tournament size can't be negative
		mutationRate   < 0.0)	// mutation rate cant' be negative
		return nullptr;

	return new TSPGeneticAlgorithm(generationSize, crossoverRate, tournamentSize, mutationRate);
}

TSPGenome** TSPGeneticAlgorithm::createNextGeneration(int citiesNumber, TSPGenome** currentGeneration) {
	// create new array for next generation
	TSPGenome** nextGeneration = new TSPGenome* [generationSize];

	// creation of initial generation
	if (currentGeneration == nullptr) {
		for (int i = 0; i < generationSize; ++i)
			nextGeneration[i] = new TSPGenome(citiesNumber);
		return nextGeneration;
	}

	// get new array of selected genomes for crossing
	TSPGenome** selection = Selection(currentGeneration);

	// each step of cycle adds 2 new children to the next generation
	for (int genomeSize = citiesNumber - 1, i = 0; i + 1 < generationSize; i += 2) {

		// get 2 random parents for crossing
		int  firstParent = rand() % crossoverRate;
		int secondParent = rand() % crossoverRate;

		// if indexes matched then simply assign second to 0 if first isn't 0, 1 otherwise
		// (to avoid long search of another random)
		if (firstParent == secondParent)
			secondParent = (firstParent == 0) ? 1 : 0;

		// allocate memory for children genomes
		gene_t* child1 = new gene_t[genomeSize];
		gene_t* child2 = new gene_t[genomeSize];

		// fill chidren genomes array with new genes via crossover operation
		int const* keyParent1 = selection[firstParent]->getKeyGenome();
		gene_t const* parent2 = selection[secondParent]->getGenome();
		//Crossover(selection[firstParent]->getKeyGenome, selection[secondParent]->getGenome, child1, child2, genomeSize);
		Crossover(keyParent1, parent2, child1, child2, genomeSize);

		// apply mutation for minor changes (mutationRate supposed to be low)
		Mutate(child1, genomeSize);
		Mutate(child2, genomeSize);

		// add new children to the next generation
		nextGeneration[i]	  = new TSPGenome(citiesNumber, child1);
		nextGeneration[i + 1] = new TSPGenome(citiesNumber, child2);
	}

	// don't forget to free memory allocated for selected parents in Selection method
	delete[] selection;
	return nextGeneration;
}

TSPGenome** TSPGeneticAlgorithm::Selection(TSPGenome** generation) {
	// create array of pointers to TSPGenome of selected genomes to take part in further crossover
	TSPGenome** selection = new TSPGenome*[crossoverRate];

	// on each step selected one genome from the generation that still wasn't chosen
	for (int winner, i = 0; i < crossoverRate; ++i) {

#ifdef TSP_GA_DEBUG
		if (generationSize - i < tournamentSize) {
			TSP_GA_DEBUG_REPORT(__FILE__, __LINE__);
			exit(-1);
		}
#endif /* TSP_GA_DEBUG */

		// set index of new winner genome of tournament
		winner = Tournament(generation, generationSize - i);
		selection[i] = generation[winner];

		// swap winner and the last genome in generation array to avoid winner's multiple selection
		SwapGenomes(generation, winner, generationSize - (i + 1));
	}

	return selection;
}

int TSPGeneticAlgorithm::Tournament(TSPGenome** generation, int remainingPopulationSize) {
	// set initial minimum
	int minimum = generation[0]->getFitness();
	int winner  = 0;

#ifdef TSP_GA_DEBUG
	if (tournamentSize > remainingPopulationSize) {
		TSP_GA_DEBUG_REPORT(__FILE__, __LINE__);
		exit(-1);
	}
#endif /* TSP_GA_DEBUG */

	// select one genome with minimum fitness from tournamentSize random genomes
	for (int i = 0; i < tournamentSize; ++i) {
		// get random genomes excluding those in the beginning that have already been compared
		int random = std::rand() % (remainingPopulationSize - i) + i;

		// compare chosen genome's fitness with current minimum
		if (minimum > generation[random]->getFitness()) {
			minimum = generation[random]->getFitness();
			winner  = random;
		}
		// swap compared random genome with the first not compared to avoid multiple comparation
		SwapGenomes(generation, random, i);
	}

	return winner;
}

void TSPGeneticAlgorithm::SwapGenomes(TSPGenome** generation, int i, int j) {
	TSPGenome* temp = generation[i];
	generation[i] = generation[j];
	generation[j] = temp;
}

/* Improved crossover CX2 algorithm for TSP (source: https://www.hindawi.com/journals/cin/2017/7430125/):
Step 1. Choose two parents for mating.
Step 2. Select 1st bit from second parent as a 1st bit of first offspring.
Step 3. The selected bit from Step  2 would be found in first parent and pick the exact same position bit
		which is in second parent and that bit would be found again in the first parent and, finally,
		the exact same position bit which is in second parent will be selected for 1st bit of second offspring.
Step 4. The selected bit from Step  3 would be found in first parent and pick the exact same position bit
		which is in second parent as the next bit for first offspring.
		(Note: for the first offspring, we choose bits only with one move and two moves for second offspring’s bits.)
Step 5. Repeat Steps  3 and 4 till 1st bit of first parent will not come in second offspring (complete a cycle)
		and process may be terminated.
Step 6. If some bits are left, then the same bits in first parent and in second offspring till now and vice versa
		are left out from both parents. For remaining bits repeat Steps  2, 3, and 4 to complete the process.*/
void TSPGeneticAlgorithm::Crossover(int const* keyParent1, gene_t const* parent2, gene_t* child1, gene_t* child2, int genomeSize) {
	// create array of indexes of genes in parent genome, added to genome of first child
	bool* checked = new bool[genomeSize];

	// and initialize it
	for (int i = 0; i < genomeSize; ++i)
		checked[i] = false;

	// here algorithm starts with step 1 already accomplished with method call
	int firstGeneIndex = 0;
	int  currGeneIndex = 0;

	// step 2: select first bit of first offspring
	child1[0] = parent2[currGeneIndex];
	checked[0] = true;

	// step 3: select first bit of second offspring
	currGeneIndex = keyParent1[parent2[currGeneIndex] - 1];
	currGeneIndex = keyParent1[parent2[currGeneIndex] - 1];
	child2[0] = parent2[currGeneIndex];

	// repeat steps 3 - 4: select next bits for children
	for (int i = 1; i < genomeSize; ++i) {

		// step 5: start new cycle
		if (keyParent1[child2[i - 1] - 1] == firstGeneIndex) {
			int j;
			for (j = 0; j < genomeSize; ++j) {
				// step 6: repeat steps 2 and 3 - 4 for remaining bits
				if (checked[j] == false) {
					// set first gene index of new cycle
					firstGeneIndex = j;
					currGeneIndex  = j;

#ifdef TSP_GA_DEBUG
					if (checked[currGeneIndex] == true) {
						TSP_GA_DEBUG_REPORT(__FILE__, __LINE__);
						exit(-1);
					}
#endif /* TSP_GA_DEBUG */

					// step 2: select next bit of first offspring
					child1[i] = parent2[currGeneIndex];
					checked[currGeneIndex] = true;

					// step 3: select next bit of second offspring
					currGeneIndex = keyParent1[parent2[currGeneIndex] - 1];
					currGeneIndex = keyParent1[parent2[currGeneIndex] - 1];
					child2[i] = parent2[currGeneIndex];
					break;
				}
			}

#ifdef TSP_GA_DEBUG
			if (j == genomeSize) {
				TSP_GA_DEBUG_REPORT(__FILE__, __LINE__);
				exit(-1);
			}
#endif /* TSP_GA_DEBUG */

			continue;
		}
		
		// step 4
		currGeneIndex = keyParent1[parent2[currGeneIndex] - 1];

#ifdef TSP_GA_DEBUG
		if (checked[currGeneIndex] == true) {
			TSP_GA_DEBUG_REPORT(__FILE__, __LINE__);
			exit(-1);
		}
#endif /* TSP_GA_DEBUG */

		child1[i] = parent2[currGeneIndex];

		// mark current parent's bit as checked
		checked[currGeneIndex] = true;

		// step 3
		currGeneIndex = keyParent1[parent2[currGeneIndex] - 1];
		currGeneIndex = keyParent1[parent2[currGeneIndex] - 1];
		child2[i] = parent2[currGeneIndex];

	}

#ifdef TSP_GA_DEBUG
	gene_t* temp = new gene_t[genomeSize];

	std::memcpy(temp, child1, sizeof(gene_t) * genomeSize);
	std::sort(temp, temp + genomeSize);
	for (int i = 0; i < genomeSize - 1; ++i) {
		if (temp[i] + 1 != temp[i + 1]) {
			TSP_GA_DEBUG_REPORT(__FILE__, __LINE__);
			delete[] temp;
			exit(-1);
		}
	}

	std::memcpy(temp, child2, sizeof(gene_t) * genomeSize);
	std::sort(temp, temp + genomeSize);
	for (int i = 0; i < genomeSize - 1; ++i) {
		if (temp[i] + 1 != temp[i + 1]) {
			TSP_GA_DEBUG_REPORT(__FILE__, __LINE__);
			delete[] temp;
			exit(-1);
		}
	}

	delete[] temp;
#endif /* TSP_GA_DEBUG */
}

void TSPGeneticAlgorithm::Mutate(gene_t* genome, int genomeSize) {
#ifdef TSP_GA_DEBUG
	if (genome == nullptr) {
		TSP_GA_DEBUG_REPORT(__FILE__, __LINE__);
		exit(-1);
	}
#endif /* TSP_GA_DEBUG */

	// generate random double number from [0, 1] interval
	double random = (double)rand() / ((double)RAND_MAX + 1.0);

	// if this value exceeds mutation rate
	if (random > mutationRate) {
		int firstIndex  = rand() % genomeSize;
		int secondIndex = rand() % genomeSize;

		// if indexes matched then simply assign second to 0 gene if first isn't 0, 1 otherwise
		// (to avoid long search of another random for small genomes)
		if (firstIndex == secondIndex) {
			secondIndex = (firstIndex == 0) ? 1 : 0;
		}

		// swap genes with firstIndex and secondIndex
		gene_t temp = genome[firstIndex];
		genome[firstIndex] = genome[secondIndex];
		genome[secondIndex] = temp;
	}
}

#ifdef TSP_GA_DEBUG
#undef TSP_GA_DEBUG
#endif /* TSP_GA_DEBUG */


//TSPGenome** TSPGeneticAlgorithm::createNextGeneration(int citiesNumber, TSPGenome** currentGeneration) {
//	TSPGenome** nextGeneration = new TSPGenome * [generationSize];
//
//	if (currentGeneration == nullptr) {
//		for (int i = 0; i < generationSize; ++i)
//			nextGeneration[i] = new TSPGenome(citiesNumber);
//		return nextGeneration;
//	}
//
//	int genomeSize = citiesNumber - 1;
//	TSPGenome** selection = Selection(currentGeneration);
//	int(*crossedPairs)[2] = new int[crossoverRate][2];
//
//	for (int i = 0; i < crossoverRate; i += 2) {
//
//		int firstParent = rand() % crossoverRate;
//
//		int secondParent = rand() % crossoverRate;
//
//		for (int j = 0; j < i; ++j) {
//			// if such pair of parents has already been crossed 
//			// (it's enough no check their roles as first and second)
//			if (crossedPairs[j][0] == firstParent && crossedPairs[j][1] == secondParent) {
//				// if indexes matched then simply assign second to 0 if first isn't 0, 1 otherwise
//				// (to avoid long search of another randoms)
//				// P.S. this won't save from repeat chance though
//				secondParent = (firstParent == 0) ? 1 : 0;
//			}
//		}
//
//		gene_t* child1 = new gene_t[genomeSize];
//		gene_t* child2 = new gene_t[genomeSize];
//
//		Crossover(selection[firstParent]->getKeyGenome, selection[secondParent]->getGenome, child1, child2, genomeSize);
//		crossedPairs[i][0] = firstParent;
//		crossedPairs[i][1] = secondParent;
//
//		Mutate(child1, genomeSize);
//		Mutate(child2, genomeSize);
//
//		nextGeneration[i] = new TSPGenome(citiesNumber, child1);
//		nextGeneration[i + 1] = new TSPGenome(citiesNumber, child2);
//	}
//
//	return nextGeneration;
//}