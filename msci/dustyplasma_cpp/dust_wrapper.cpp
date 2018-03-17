#include <iostream>
#include <vector>
#include <string>
#include "dust.cpp"

extern "C" void simulate (int meth, int init_iterationsInput, int iterationsBInput, int numParticles, double *x, double *y, double *z, double* gr, double * gz, double * evalr, double * evalz, float maxr, float sheathsep, float sephor1, float sephor2, float pointone, double *output, double *output2, double *output3) {
    METH = meth;
    iterationsB = iterationsBInput;
    init_iterations = init_iterationsInput;
    n_particles = numParticles;

    std::vector<Dust> dustlist(n_particles);


//	Dust dust0, dust1;
//	std::vector<float> positiond0(3);
//	positiond0={-3*lambdaD,0,0.0003825734};
//	std::vector<float> positiond1(3);
//	positiond1={3*lambdaD,3*lambdaD,0.0003825734};
//	std::vector<float> posd0 = positiond0;
//	std::vector<float> posd1 = positiond1;
//	std::vector<float> vel1 = {0,0,0};
//	std::vector<float> acc1 = {0,0,0};
//	std::vector<float> vel = vel1;
//	std::vector<float> acc = acc1;
//	dust0.set_values(true, md, radd, lambdaD, phia, Zd*e, OMEGATAU, posd0, vel1, acc1, positiond0, vel, acc, multifields);
	// dust1.set_values(true, md, radd, lambdaD, phia, Zd, OMEGATAU, posd1, vel1, acc1, positiond1, vel, acc, multifields);
	// std::vector<float> test(3);
	// std::vector<float> testB(3);
	// testB = {0,0,0.014};
	// test=dust1.combinedrift(testB);
	// std::cout << test[0] << "," << test[1] << "," << test[2] << std::endl;


	// Testing analysis functions work
	DustAnalysis analysis;
	int arrlen = n_particles;
	std::vector<float> initpositionsx(x, x + arrlen);
	std::vector<float> initpositionsy(y, y + arrlen);
	std::vector<float> initpositionsz(z, z + arrlen);

    RMAX = maxr;
    SEPARATIONSHEATH = sheathsep;
    SEPARATIONHOR1 =sephor1;
    SEPARATIONHOR2 = sephor2;
    FIRSTPOINT = pointone;

    int flattenlen = int(15*126);
    std::vector<float> GRIDR(gr, gr + flattenlen);
    std::vector<float> GRIDZ(gz, gz + flattenlen);
    std::vector<float> EVALSR(evalr, evalr + flattenlen);
    std::vector<float> EVALSZ(evalz, evalz + flattenlen);

	analysis.set_values(dustlist, pairs, newEVALSr, newEVALSz, newGRIDr, newGRIDz, initpositionsx, initpositionsy, initpositionsz, positionx, positiony, positionz,EVALSR, EVALSZ, GRIDR, GRIDZ, t, FIRSTPOINT, RMAX, SEPARATIONSHEATH, SEPARATIONHOR1, SEPARATIONHOR2, n_particles, iterationsB, init_iterations);
	analysis.changefieldshape();
	analysis.create_particles();
	std::vector<Dust> testname(n_particles);
	testname = analysis.getdustnames();
	analysis.create_pairs();
	analysis.interact_and_iterate(METH);

    // Output
//	std::vector<float> xpos;
//	std::vector<float> ypos;
//	std::vector<float> zpos;
//
//	xpos = analysis.getallpositionx();
//	ypos = analysis.getallpositiony();
//	zpos = analysis.getallpositionz();
//
//	std::copy(xpos.begin(), xpos.end(), output);
//	std::copy(ypos.begin(), ypos.end(), output2);
//	std::copy(zpos.begin(), zpos.end(), output3);

}
