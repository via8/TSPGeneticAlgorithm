#ifndef TSP_GENOME_H
#define TSP_GENOME_H

using gene_t = int;

class TSPGenome final {
private:
	int genomeSize;				// size of genome which equals to citiesNumber
	gene_t* genome;				// ordered array of numbers of visited cities (from 0 to citiesNumber - 1)
	int* keyGenome;				// ordered array of keys for genome array (from 0 to citiesNumber - 1)
	int fitness;				// price of traversing cities according to current genome

	void randomizeGenome();
	void fillKeyGenome();

	TSPGenome(TSPGenome const&)				 = delete;
	TSPGenome& operator = (TSPGenome const&) = delete;

public:
	TSPGenome(int citiesNumber, gene_t* genome = nullptr, int* keyGenome = nullptr);
	~TSPGenome() { delete[] genome; delete[] keyGenome; }
	int getFitness() const;
	void recalculateFitness(int const** pricesMatrix);
	gene_t const* getGenome() const { return	genome; }
	int const* getKeyGenome() const { return keyGenome; }
};

#endif /* TSP_GENOME_H */