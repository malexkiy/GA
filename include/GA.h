#ifndef _GA_H_
#define _GA_H_

#include <cstdint>
#include <vector>

class GA
{
private:
	struct Individual {
		double fitness;
		int32_t chrom;

		Individual(){}
		Individual(const GA&){}

		bool operator ()(const Individual& lhs, const Individual& rhs)
		{
			return lhs.fitness < rhs.fitness;
		}
	};

	struct Parents {
		uint32_t parentId[2];
	};

	struct Children {
		Individual child[2];
	};

	std::vector<double> _coefs;
	Individual *_population;
	Individual *_childPopulation;
	Individual *_newPopulation;
	uint32_t _popSize;
	uint32_t _childPopSize;
	uint32_t _maxIters;
	int32_t _minVal;
	int32_t _maxVal;
	uint32_t _steps;
	uint32_t _bits;
	float _mutationRate;

	double _coef;

	void GeneratePopulation();
	void Fitness(Individual & individ);
	Parents Selection();

	float CrossoverProb();
	Children* Crossover(Parents & parents, float prob);

	float MutationProb();
	void Mutation(Individual & indiv, float prob);

	int32_t Rand32();

	double calc(uint32_t chrom);
	bool FitnessSort(const Individual & indiv1, const Individual & indiv2);
	void SortByFitness(Individual *population, uint32_t count);

	void PopulationMerge();

public:
	GA(std::vector<double> coefs, uint32_t popSize, uint32_t maxIters, int32_t minVal, int32_t maxVal, uint32_t steps, float mutationRate);
	~GA();

	void Run();

	double BestFitness();
};

#endif
