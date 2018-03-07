#include <iostream>
#include <vector>
#include "dust.cpp"


extern "C" void simulate (int init_iterationsInput, int iterationsBInput, int numParticles, double *output) {

    iterationsB = iterationsBInput;
    init_iterations = init_iterationsInput;
    n_particles = numParticles;

	Dust dust0, dust1;
	std::vector<float> positiond0(3);
	positiond0={-3*lambdaD,0,0.0003825734};
	std::vector<float> positiond1(3);
	positiond1={3*lambdaD,3*lambdaD,0.0003825734};
	std::vector<float> posd0 = positiond0;
	std::vector<float> posd1 = positiond1;
	std::vector<float> vel1 = {0,0,0};
	std::vector<float> acc1 = {0,0,0};
	std::vector<float> vel = vel1;
	std::vector<float> acc = acc1;
	dust0.set_values(true, md, radd, lambdaD, phia, Zd*e, OMEGATAU, posd0, vel1, acc1, positiond0, vel, acc, multifields);
	// dust1.set_values(true, md, radd, lambdaD, phia, Zd, OMEGATAU, posd1, vel1, acc1, positiond1, vel, acc, multifields);
	// std::vector<float> test(3);
	// std::vector<float> testB(3);
	// testB = {0,0,0.014};
	// test=dust1.combinedrift(testB);
	// std::cout << test[0] << "," << test[1] << "," << test[2] << std::endl;


	// Testing analysis functions work
	DustAnalysis analysis;
	initpositionsx = {3*lambdaD,-3*lambdaD};
	initpositionsy = {0,0};
	initpositionsz = {0.0003825734,0.0003825734};
	analysis.set_values(dustlist, pairs, initpositionsx, initpositionsy, initpositionsz, positionx, positiony, positionz,t, n_particles, iterationsB, init_iterations);
	analysis.create_particles();
	std::vector<Dust> testname(n_particles);
	testname = analysis.getdustnames();
	analysis.create_pairs();
	analysis.interact_and_iterate();

    // Output
	std::vector<float> xpos;
	xpos = analysis.getallpositionx();

	std::copy(xpos.begin(), xpos.end(), output);

}
