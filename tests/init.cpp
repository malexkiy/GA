#include <GA.h>
#include <vector>
#include <catch.hpp>

SCENARIO("ga test", "[test]") {
	std::vector<double> coefs = { 5., -24., 17., -11. / 3, 1. / 4 };
	uint32_t popSize = 100;
	uint32_t maxIters = 1000;
	int32_t minVal = 0, maxVal = 7;
	uint32_t steps = 32; //MUST BE POWER OF 2
	float mutationRate = 0.25f;

	GA ga(coefs, popSize, maxIters, minVal, maxVal, steps, mutationRate);

	ga.Run();

	REQUIRE(ga.BestFitness() == -5.4166666666666661);
}
