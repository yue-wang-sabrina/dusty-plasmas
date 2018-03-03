#include <iostream>  
#include <vector>
#include <math.h>
#include <string>
#include <fstream>
#include <cmath>
#include <utility>
#include <unistd.h>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <random>
#include <chrono>
#include <iomanip>
#include <thread>

const double PI = std::atan(1.0)*4;
float radd = 1.5 * pow(10,-6);
float mass = 1000. * (4 / 3) * PI * pow(radd,3);
float omega = (1.60217662 * pow(10,-19)) * 0.014 / (39.948 * 1.66053904 * pow(10,-27));
float mi = 39.948 * 1.66053904 * pow(10,-27);
float me = 9.10938356 * pow(10,-31);
float tau = 1*pow(10,-5);
float Babs = 0.014;
float Eabs = 1.0;
float dt = 1*pow(10,-7);
float kb = 1.38064852 * pow(10,-23);
float Ti = 310;
float vinit = sqrt(kb*Ti/mi);
float ioncharge = 1.60217662 * pow(10,-19);
float phia =  -9.0126016670216487;
float Zd =  -9388.3579633332938;
float gammadamp = 5000;
float Te = 46000;
float e0 = 8.85418782 * pow(10,-12);
float e = 1.60217662 * pow(10 ,-19);
float ne0 = 1. * pow(10, 15);
float ni0 = ne0;
float lambdade = kb * Te * e0 / pow((ne0 * pow(e,2)),0.5);
float lambdadi = kb * Ti * e0 / pow((ni0 * pow(e,2)), 0.5);
float lambdaD = 1. / pow(1. / pow(lambdadi, 2) + 1. / pow(lambdade,2), 0.5);
float sheathd = 10*lambdaD;
float boxr = 523*lambdaD;
float dt = 0.0001;
float g=-9.8;
float electrodeV = abs(kb * Te / (2 * e)) * log(2 * PI * me / mi);
float wallV = electrodeV; 
int iterations = 1*pow(10,7);
std::vector<float> pos = {0,0,0};
std::vector<float> vel = {0,0,0};
std::vector<float> acc = {0,0,0};
std::vector<float> pos1 = {0,0,0};
std::vector<float> vel1 = {0,0,0};
std::vector<float> acc1 = {0,0,0};
std::vector<float> positionx(iterations+1);
std::vector<float> positiony(iterations+1);
std::vector<float> positionz(iterations+1);

float Normalise(std::vector<float> v1, float output_float){
	output_float = sqrt(1*pow(v1[0],2) + 1*pow(v1[1],2) + 1*pow(v1[2],2));
	return output_float;
};

std::vector<float> Multiplyscalar(std::vector<float> v1, float scalar, std::vector<float> result_v){
std::transform(v1.begin(), v1.end(), result_v.begin(), std::bind1st(std::multiplies<float>(),scalar));
return result_v;
};

class Dust{
	float md, radius, lambdaD, phia, initcharge;
	std::vector<float> initpos, initvel, initacc, pos, vel, acc;

public:
	void set_values(float, float, float, float, float, std::vector<float>,std::vector<float>,std::vector<float>,std::vector<float>,std::vector<float>,std::vector<float>);

	std::vector<float> getselfpos(){
		return pos;
	}

	std::vector<float> getselfvel(){
		return vel;
	}

	std::vector<float> getselfacc(){
		return acc;
	}

	void updateEuler(){


	}

	std::vector<float> damping(){
		float magvel;
		std::vector<float> acctemp(3);
		magvel = Normalise(vel, magvel);
	if (magvel == float(0)){
		acctemp = {0,0,0};
		return acctemp;
	}
	else{
		float multiplyfactor;
		multiplyfactor = -1*gammadamp*pow(magvel,2);
		std::vector<float> acctemp(3);
		acctemp = Multiplyscalar(acctemp,1./magvel,acctemp);
		acctemp = Multiplyscalar(acctemp, multiplyfactor, acctemp);
		return acctemp;
	}
	}

	std::vector<float> verticalion(){
		if (pos[2]<sheathd){
			float mach = 5;
			float beta = abs(Zd*e/(Ti*lambdade));
			std::vector<float> acctemp(3)={0,0,0};
			acctemp[3] = (pow(Ti / e), 2) * log(lambdade * pow(mach, 2) / (beta * lambdadi)) * pow(beta, 2) / pow(mach,2) / md)
			return acctemp;
		}
		else{
			std::vector<float> acctemp(3) = {0,0,0};
			return acctemp;
		}

	}

	std::vector<float> selffield(Dust g2){

		return g2.getselfpos();
	}




};

void Dust::set_values(float a, float b, float c, float d, float e, std::vector<float> f, std::vector<float> g, std::vector<float> h, std::vector<float> i, std::vector<float> j, std::vector<float> k){
	md = a;
	radius = b;
	lambdaD = c;
	phia = d;
	initcharge = e;
	initpos = f;
	initvel = g;
	initacc = h;
	pos = i;
	vel = j;
	acc = k;
};




int main(){

	Dust dust0, dust1;
	dust0.set_values(md, radd, lambdaD, phia, Zd, pos1, vel1, acc1, std::vector<float> positiond0={0,0,0}, vel, acc);
	dust0.set_values(md, radd, lambdaD, phia, Zd, pos1, vel1, acc1, std::vector<float> positiond1={0,0,0}, vel, acc);

	std::vector<float> test(3);
	test = dust0.selffield(dust1);
	std::cout << test[0] << std::endl;


	return 0;
}














