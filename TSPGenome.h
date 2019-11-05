#ifndef TSP_GENOME_H
#define TSP_GENOME_H

class TSPGenome final {
private:
	int* genome;		// ordered array of numbers of visited cities
	int* indexes;		// ordered array of keys for genome array
	int genomeSize;		// size of genome which equals to citiesNumber
	int fitness;		// price of traversing cities according to current genome

	TSPGenome(TSPGenome const&)				 = delete;
	TSPGenome& operator = (TSPGenome const&) = delete;

public:
	TSPGenome(int citiesNumber);
	~TSPGenome() { delete[] genome; delete[] indexes; }
	void randomizeGenome();
	void setGene(int index, int gene) { genome[index] = gene; indexes[gene] = index; }
	int getGene(int index) const { return genome[index]; }
	int getIndex(int gene) const { return indexes[gene]; }
	int getFitness() const { return fitness; }
	void recalculateFitness(int const** pricesMatrix);
	void mutate(double mutationRate);
};

#endif /* TSP_GENOME_H */