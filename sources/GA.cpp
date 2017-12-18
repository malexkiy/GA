#include "GA.h"
#include <utility>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <functional>


GA::GA(std::vector<double> coefs, uint32_t popSize, uint32_t maxIters, int32_t minVal, int32_t maxVal, uint32_t steps, float mutationRate)
	:_coefs(coefs),
	_popSize(popSize),
	_maxIters(maxIters),
	_minVal(minVal),
	_maxVal(maxVal),
	_steps(steps),
	_mutationRate(mutationRate)
{
	srand(time(nullptr));

	_bits = 1;
	while (!((_steps >> _bits) & 1))
	{
		_bits++;
	}

	_coef = (_maxVal - _minVal + 1) / (float)_steps;

	_population = nullptr;
	_childPopulation = nullptr;
	_newPopulation = nullptr;

	_population = new Individual[_popSize];
	_childPopulation = new Individual[_popSize];
	_newPopulation = new Individual[_popSize];
}


GA::~GA()
{
	if (_population)
	{
		delete[] _population;
	}
	if (_childPopulation)
	{
		delete[] _childPopulation;
	}
	if (_newPopulation)
	{
		delete[] _newPopulation;
	}
}


double GA::calc(uint32_t chrom)
{
	return _minVal + _coef*chrom;
}


bool GA::FitnessSort(const Individual & lhs, const Individual & rhs)
{
	return lhs.fitness < rhs.fitness;
}


void GA::SortByFitness(Individual *population, uint32_t count)
{
	std::sort(population, population + count, Individual(*this));
}


void GA::PopulationMerge()
{
	uint32_t iPop = 0, iChPop = 0, iNewPop = 0;

	while (iNewPop < _popSize)
	{
		if (iPop == _popSize)
		{
			while (iChPop != _childPopSize)
			{
				_newPopulation[iNewPop++] = _childPopulation[iChPop++];
			}

			return;
		}

		if (iChPop == _childPopSize)
		{
			while (iPop != _popSize)
			{
				_newPopulation[iNewPop++] = _population[iPop++];
			}

			return;
		}

		_newPopulation[iNewPop++] = (_population[iPop].fitness < _childPopulation[iChPop].fitness) ? _population[iPop++] : _childPopulation[iChPop++];

	}
}


void GA::GeneratePopulation()
{
	for (uint32_t i = 0; i < _popSize; i++)
	{
		Individual individ;

		individ.chrom = Rand32() & (_steps - 1);

		_population[i] = individ;
	}
}


void GA::Fitness(GA::Individual & individ)
{
	double chrom = calc(individ.chrom);
	individ.fitness = 0;

	for (uint32_t i = 0; i < _coefs.size(); i++)
	{
		individ.fitness += _coefs[i] * pow(chrom, i);
	}
}


GA::Parents GA::Selection()
{
	Parents parents;

	parents.parentId[0] = Rand32() % _popSize;
	parents.parentId[1] = Rand32() % _popSize;

	return parents;
}


float GA::MutationProb()
{
	return rand() / (float)RAND_MAX;
}


float GA::CrossoverProb()
{
	return rand() / (float)RAND_MAX;
}


GA::Children* GA::Crossover(GA::Parents & parents, float prob)
{
	Children* children = nullptr;
	uint32_t chromP1, chromP2, maskH, maskL;
	uint32_t crossPoint;

	chromP1 = _population[parents.parentId[0]].chrom;
	chromP2 = _population[parents.parentId[1]].chrom;

	if (chromP1 != chromP2)
	{
		children = new Children();

		crossPoint = rand() % _bits;

		maskL = 0xFFFFFFFF << crossPoint;
		maskH = ~maskL;

		children->child[0].chrom = (chromP1 & maskH) | (chromP2 & maskL);

		children->child[1].chrom = (chromP2 & maskH) | (chromP1 & maskL);
	}

	return children;
}


int32_t GA::Rand32()
{
	return (rand() << 16 | rand());
}


void GA::Mutation(Individual& indiv, float prob)
{
	uint32_t mutPoint;

	if (prob <= _mutationRate)
	{
		mutPoint = rand() % _bits;
		indiv.chrom ^= (1 << mutPoint);
	}
}


double GA::BestFitness()
{
	return _population[0].fitness;
}


void GA::Run()
{
	using std::swap;
	uint32_t iters = 0, mates;

	GeneratePopulation();

	for (uint32_t i = 0; i < _popSize; i++)
	{
		Fitness(_population[i]);
	}

	while (iters < _maxIters)
	{
		_childPopSize = 0;
		mates = 0;

		while (_childPopSize < _popSize && mates < _popSize)
		{
			Parents parents = Selection();

			Children* children = Crossover(parents, CrossoverProb());
			if (children)
			{
				Mutation(children->child[0], MutationProb());
				Mutation(children->child[1], MutationProb());

				_childPopulation[_childPopSize++] = children->child[0];
				_childPopulation[_childPopSize++] = children->child[1];

				delete children;
			}

			mates += 2;
		}
		for (uint32_t i = 0; i < _childPopSize;i++)
		{
			Fitness(_childPopulation[i]);
		}

		SortByFitness(_population, _popSize);
		SortByFitness(_childPopulation, _childPopSize);

		PopulationMerge();

		swap(_newPopulation, _population);

		iters++;
	}
}
