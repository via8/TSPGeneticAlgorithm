#include "TSPGeneticAlgorithm.h"
#include <cstdlib>

#define TSP_GA_DEBUG

#ifdef TSP_GA_DEBUG
#include <iostream>
#include <algorithm>
#define TSP_GA_DEBUG_REPORT(file, line) (std::cout << "ERROR: file " << file << ", line " << line << std::endl)
#endif /* TSP_GA_DEBUG */

TSPGeneticAlgorithm::TSPGeneticAlgorithm(int generations, int generationSize, int crossoverRate, int tournamentSize, double mutationRate) :
	generations(generations),
	generationSize(generationSize),
	tournamentSize(tournamentSize),
	crossoverRate(crossoverRate),
	mutationRate(mutationRate) {
}

TSPGeneticAlgorithm* TSPGeneticAlgorithm::create(int generations, int generationSize, int crossoverRate, int tournamentSize, double mutationRate) {
	// make sure that generation size is even as crossover operation always builds pairs of children
	if (generationSize % 2 != 0)
		generationSize++;
	
	if (generations < 0 ||	// number of generations can't be negative
		generationSize - crossoverRate < tournamentSize - 1 ||	// otherwise it would be impossible to arrange tournament
		generationSize < MIN_GEN_SIZE  ||	// minimum possible size of generation to run the algorithm
		generationSize < crossoverRate ||	// impossible to cross more genomes then we have in generation
		tournamentSize < 0 ||	// tournament size can't be negative
		mutationRate   < 0.0)	// mutation rate cant' be negative
		return nullptr;

	return new TSPGeneticAlgorithm(generations, generationSize, crossoverRate, tournamentSize, mutationRate);
}

TSPGenome** TSPGeneticAlgorithm::createNextGeneration(int citiesNumber, TSPGenome** currentGeneration) {
	// create new array for next generation
	TSPGenome** nextGeneration = new TSPGenome* [generationSize];

	// creation of initial generation
	if (currentGeneration == nullptr) {
		for (int i = 0; i < generationSize; ++i) {
			nextGeneration[i] = new TSPGenome(citiesNumber);
			nextGeneration[i]->randomizeGenome();
		}
		return nextGeneration;
	}

	// get new array of selected genomes for crossing
	TSPGenome** selected = selection(currentGeneration);

	// each step of cycle adds 2 new children to the next generation
	for (int genomeSize = citiesNumber - 1, i = 0; i + 1 < generationSize; i += 2) {

		// get 2 random parents for crossing
		int parent1 = rand() % crossoverRate;
		int parent2 = rand() % crossoverRate;

		// if indexes matched then simply assign second to 0 if first isn't 0, 1 otherwise
		// (to avoid long search of another random)
		if (parent1 == parent2)
			parent2 = (parent1 == 0) ? 1 : 0;

		// create chidren genomes and fill them with new genes via crossover operation
		TSPGenome* child1 = new TSPGenome(citiesNumber);
		TSPGenome* child2 = new TSPGenome(citiesNumber);
		crossover(selected[parent1], selected[parent2], child1, child2, genomeSize);

		// apply mutation for minor changes (mutationRate supposed to be low)
		child1->mutate(mutationRate);
		child2->mutate(mutationRate);

		// add new children to the next generation
		nextGeneration[i] = child1;
		nextGeneration[i + 1] = child2;
	}

	// don't forget to free memory allocated for selected parents in Selection method
	delete[] selected;
	return nextGeneration;
}

TSPGenome** TSPGeneticAlgorithm::selection(TSPGenome** generation) {
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
		winner = tournament(generation, generationSize - i);
		selection[i] = generation[winner];

		int a = selection[i]->getFitness();
		// swap winner and the last genome in generation array to avoid winner's multiple selection
		swapGenomes(generation, winner, generationSize - (i + 1));
	}

	return selection;
}

int TSPGeneticAlgorithm::tournament(TSPGenome** generation, int remainingGenerationSize) {
	// set initial winner
	int winner  = 0;

#ifdef TSP_GA_DEBUG
	if (tournamentSize > remainingGenerationSize) {
		TSP_GA_DEBUG_REPORT(__FILE__, __LINE__);
		exit(-1);
	}

	int* compared = new int[tournamentSize + 1];
	compared[0] = generation[0]->getFitness();
#endif /* TSP_GA_DEBUG */

	// select one genome with minimum fitness from tournamentSize random genomes
	for (int i = 0; i < tournamentSize; ++i) {
		// get random genomes excluding those in the end that have already been compared
		int random = std::rand() % (remainingGenerationSize - i);

		// compare chosen genome's fitness with current winner's
		if (generation[winner]->getFitness() > generation[random]->getFitness())
			winner = random;

#ifdef TSP_GA_DEBUG
		compared[i + 1] = generation[random]->getFitness();
#endif /* TSP_GA_DEBUG */

		// swap compared random genome with the last not compared to avoid multiple comparation
		swapGenomes(generation, random, remainingGenerationSize - (i + 1));

		// if winner was swapped then assign it to the new index
		if (winner == random)
			winner = remainingGenerationSize - (i + 1);
	}

#ifdef TSP_GA_DEBUG
	int min = compared[0];
	for (int i = 1; i < tournamentSize + 1; ++i) {
		if (min > compared[i])
			min = compared[i];
	}
	delete[] compared;

	if (min != generation[winner]->getFitness()) {
		TSP_GA_DEBUG_REPORT(__FILE__, __LINE__);
		exit(-1);
	}
#endif /* TSP_GA_DEBUG */

	return winner;
}

void TSPGeneticAlgorithm::swapGenomes(TSPGenome** generation, int i, int j) {
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
void TSPGeneticAlgorithm::crossover(TSPGenome* parent1, TSPGenome* parent2, TSPGenome* child1, TSPGenome* child2, int genomeSize) {
	// create array of indexes of genes in parent genome, added to genome of first child
	bool* checked = new bool[genomeSize];

	// and initialize it
	for (int i = 0; i < genomeSize; ++i)
		checked[i] = false;

	// here algorithm starts with step 1 already accomplished with method call
	int first = 0;	// index of the first parent gene processed
	int  curr = 0;	// index of the current parent gene processed

	// step 2: select first bit of first offspring
	child1->setGene(0, parent2->getGene(curr));
	checked[0] = true;

	// step 3: select first bit of second offspring
	curr = parent1->getIndex(parent2->getGene(curr));
	curr = parent1->getIndex(parent2->getGene(curr));
	child2->setGene(0, parent2->getGene(curr));

	// repeat steps 3 - 4: select next bits for children
	for (int i = 1; i < genomeSize; ++i) {

		// step 5: start new cycle
		if (parent1->getIndex(child2->getGene(i - 1)) == first) {
			int j;
			for (j = 0; j < genomeSize; ++j) {
				// step 6: repeat steps 2 and 3 - 4 for remaining bits
				if (checked[j] == false) {
					// set first gene index of new cycle
					first = j;
					curr  = j;

#ifdef TSP_GA_DEBUG
					if (checked[curr] == true) {
						TSP_GA_DEBUG_REPORT(__FILE__, __LINE__);
						exit(-1);
					}
#endif /* TSP_GA_DEBUG */

					// step 2: select next bit of first offspring
					child1->setGene(i, parent2->getGene(curr));
					checked[curr] = true;

					// step 3: select next bit of second offspring
					curr = parent1->getIndex(parent2->getGene(curr));
					curr = parent1->getIndex(parent2->getGene(curr));
					child2->setGene(i, parent2->getGene(curr));
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
		curr = parent1->getIndex(parent2->getGene(curr));

#ifdef TSP_GA_DEBUG
		if (checked[curr] == true) {
			TSP_GA_DEBUG_REPORT(__FILE__, __LINE__);
			exit(-1);
		}
#endif /* TSP_GA_DEBUG */

		child1->setGene(i, parent2->getGene(curr));

		// mark current parent's bit as checked
		checked[curr] = true;

		// step 3
		curr = parent1->getIndex(parent2->getGene(curr));
		curr = parent1->getIndex(parent2->getGene(curr));
		child2->setGene(i, parent2->getGene(curr));
	}

#ifdef TSP_GA_DEBUG
	int* temp = new int[genomeSize];

	for (int i = 0; i < genomeSize; ++i)
		temp[i] = child1->getGene(i);
	std::sort(temp, temp + genomeSize);
	for (int i = 0; i < genomeSize - 1; ++i) {
		if (temp[i] + 1 != temp[i + 1]) {
			TSP_GA_DEBUG_REPORT(__FILE__, __LINE__);
			delete[] temp;
			exit(-1);
		}
	}

	for (int i = 0; i < genomeSize; ++i)
		temp[i] = child2->getGene(i);
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

int TSPGeneticAlgorithm::findFittestGenome(TSPGenome const* const* generation) {
	int min = 0;
	for (int i = 0; i < generationSize; ++i) {
		if (generation[min]->getFitness() > generation[i]->getFitness())
			min = i;
	}

	return min;
}

int TSPGeneticAlgorithm::run(int citiesNumber, int const** pricesMatrix, int* sequence) {
	// create initial (random) generation
	TSPGenome** currGeneration = createNextGeneration(citiesNumber);

	// calculate fitness of each genome in the generation
	for (int i = 0; i < generationSize; ++i)
		currGeneration[i]->recalculateFitness(pricesMatrix);

	// save inital fittest solution
	int min = findFittestGenome(currGeneration);
	int answer = currGeneration[min]->getFitness();
	for (int i = 0; i < citiesNumber - 1; ++i)
		sequence[i] = currGeneration[min]->getGene(i);

	// start reproducing next generations
	for (int i = 0; i < generations; ++i) {
		TSPGenome** nextGeneration = createNextGeneration(citiesNumber, currGeneration);

		// recalculate fitness of each genome in the next generation
		// and free memory allocated for previous one
		for (int j = 0; j < generationSize; ++j) {
			nextGeneration[j]->recalculateFitness(pricesMatrix);
			delete currGeneration[j];
		}
		delete[] currGeneration;

		// assign new pointer
		currGeneration = nextGeneration;

		// save current fittest solution if needed
		min = findFittestGenome(currGeneration);
		if (answer > currGeneration[min]->getFitness()) {
			answer = currGeneration[min]->getFitness();
			for (int i = 0; i < citiesNumber - 1; ++i)
				sequence[i] = currGeneration[min]->getGene(i);
		}
	}

	// free memory allocated for last generation
	for (int i = 0; i < generationSize; ++i)
		delete currGeneration[i];
	delete[] currGeneration;

	return answer;
}

#ifdef TSP_GA_DEBUG
#undef TSP_GA_DEBUG
#endif /* TSP_GA_DEBUG */
