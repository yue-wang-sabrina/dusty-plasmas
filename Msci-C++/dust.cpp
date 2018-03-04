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
float radd = 1.5 * 1*pow(10,-6);
float md = 1*pow(radd,3)*PI*1000*4/3;
float omega = (1.60217662 * pow(10,-19)) * 0.014 / (39.948 * 1.66053904 * pow(10,-27));
float OMEGATAU = 0.01; // From Konopka experimentalbasis paper
float mi = 39.948 * 1.66053904 * pow(10,-27);
float me = 9.10938356 * pow(10,-31);
// float tau = 1*pow(10,-5);
float Babs = 0.014;
float Eabs = 1.0;
float dt = 0.0001;
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
float lambdade = 1*pow(kb * Te * e0 / (ne0 * pow(e,2)),0.5);
float lambdadi = 1*pow(kb * Ti * e0 / (ni0 * pow(e,2)), 0.5);
float lambdaD = 1*pow(1. / (1. / pow(lambdadi, 2) + 1. / pow(lambdade,2)), 0.5);
float sheathd = 10*lambdaD;
float boxr = 523*lambdaD;
float g=-9.8;
float electrodeV = abs((kb * Te / (2 * e)) * log(2 * PI * me / mi));
float wallV = electrodeV; 
float radinfluence = 10*lambdaD;
int iterations = 1*pow(10,5);

std::vector<float> positionx(iterations+1);
std::vector<float> positiony(iterations+1);
std::vector<float> positionz(iterations+1);

float Magnitude(std::vector<float> v1, float output_float){
	output_float = sqrt(1*pow(v1[0],2) + 1*pow(v1[1],2) + 1*pow(v1[2],2));
	return output_float;
};

std::vector<float> Multiplyscalar(std::vector<float> v1, float scalar, std::vector<float> result_v){
std::transform(v1.begin(), v1.end(), result_v.begin(), std::bind1st(std::multiplies<float>(),scalar));
return result_v;
};

std::vector<float> Minusvectors(std::vector<float> v1, std::vector<float> v2, std::vector<float> minus_v){
std::transform (v1.begin(), v1.end(), v2.begin(), minus_v.begin(), std::minus<float>());
return minus_v;
}

std::vector<float> Addvectors(std::vector<float> v1, std::vector<float> v2, std::vector<float> minus_v){
std::transform (v1.begin(), v1.end(), v2.begin(), minus_v.begin(), std::plus<float>());
return minus_v;
}

std::vector<float> Crossproduct(std::vector<float> v1, std::vector<float> v2, std::vector<float> cross_p) {
cross_p[0] = (v1[1]*v2[2]) - (v1[2]*v2[1]);
cross_p[1] = -((v1[0]*v2[2]) - (v1[2]*v2[0]));
cross_p[2] = (v1[0]*v2[1]) - (v1[1]*v2[0]);
return cross_p;
}


class Dust{
	bool Bswitch;
	float md, radd, lambdaD, phia, initcharge, OMEGATAU;
	std::vector<float> initpos, initvel, initacc, pos, vel, acc, multifields;

public:
	void set_values(bool, float, float, float, float, float, float, std::vector<float>,std::vector<float>,std::vector<float>,std::vector<float>,std::vector<float>,std::vector<float>, std::vector<float>);

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
		magvel = Magnitude(vel, magvel);
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
			std::vector<float> acctemp(3);
			acctemp={0,0,0};
			acctemp[3] = (1*pow((Ti / e), 2) * log(lambdade * pow(mach, 2) / (beta * lambdadi)) * pow(beta, 2) / pow(mach,2)) / md;
			return acctemp;
		}
		else{
			std::vector<float> acctemp(3);
			acctemp = {0,0,0};
			return acctemp;
		}

	}

	std::vector<float> selffield(Dust g2){
		std::vector<float> g2pos(3);
		g2pos = g2.getselfpos();
		std::vector<float> diff(3);
		diff = Minusvectors(g2pos,pos,diff);
		float d = Magnitude(diff,d);
		std::vector<float> r(3);
		r = Minusvectors(pos,g2pos,r);
		std::vector<float> unitr(3);
		unitr = Multiplyscalar(r,1./d,unitr);
		float E = abs(phia/lambdaD)*exp(-1*(d-2*radd)/lambdaD);
		std::vector<float> acctemp(3);
		acctemp = Multiplyscalar(unitr, E,acctemp);
		return acctemp;
		
	}

	void selffieldmany(std::vector<float> E){
		multifields = Addvectors(multifields,E,multifields);
	}
 	
 	std::vector<float> sheathfield(){
 		std::vector<float> acctemp;
 		if (pos[2] >= sheathd & pos[2]>0){
 			acctemp = {0,0,0};
 			return acctemp;
 		}
 		else if(0< pos[2] & pos[2] < sheathd){
 			float field;
            field = abs(2 * (1. - pos[2] / sheathd) * (electrodeV / sheathd));
            acctemp = {0,0,field};
            return acctemp;
 		}
 		else{
 			throw std::invalid_argument( "Dust fell out of cylinder! Origin: sheathfield()" );
 		}
 		}

 	std::vector<float> radialfield(){
 		std::vector<float> acctemp(3);
 		std::vector<float> r(3);
 		r = {-pos[0], -pos[1], 0};
 		float magr;
 		magr = Magnitude(r,magr);

 		if (pos[2] < sheathd){
 			float magE = (wallV/(1*pow(boxr,2))) * 2 * sqrt(1*pow(pos[0],2) + 1*pow(pos[1],2));
 			if (float(1*pow(r[0],2)+1*pow(r[1],2)) == float(0)){
 				acctemp = {0,0,0};
 				return acctemp;
 			}
 			else{
 				acctemp = Multiplyscalar(r,magE/magr,acctemp);
 				return acctemp;
 			}
 		}
 		else if(pos[2] > sheathd){
 			acctemp = {0,0,0};
 			return acctemp;
 		}
 		else {
 			throw std::invalid_argument( "Dust fell out of cylinder! Origin: radialfield()");
 		}

 	}

 	bool intergraind(Dust g2){
 		std::vector<float> diff(3);
 		diff = Minusvectors(pos,g2.getselfpos(),diff);
 		float magdiff;
 		magdiff = Magnitude(diff,magdiff);
 		return magdiff <= radinfluence;
 	}

 	std::vector<float> EXBacchybrid(std::vector<float>B, std::string method){
 		std::vector<float> acctemp(3);
 		if (Bswitch && pos[2] < sheathd){
 			float magB;
 			magB = Magnitude(B,magB);
 			float vT;
 			vT = sqrt(kb * Ti / mi);
 			std::vector<float> vdrift(3);
 			std::vector<float> radialinstance(3);
 			std::vector<float> sheathinstance(3);
 			radialinstance = radialfield();
 			sheathinstance = sheathfield();
 	 		std::vector<float> E(3);
 			E = Addvectors(radialinstance, sheathinstance, E);
 			vdrift = Crossproduct(E, B, vdrift);
 			vdrift = Multiplyscalar(vdrift,1./(1.*pow(magB,2)),vdrift);
 			
 			float OMEGA = abs(e*magB/mi);
 			float TAU = OMEGATAU/OMEGA;

 			if (method == "factor"){
 				float factor = 1*pow(OMEGATAU,2)/(1. + 1*pow(OMEGATAU,2));
 				vdrift = Multiplyscalar(vdrift,factor,vdrift);
 			}

 			else if (method == "derivation"){
 				float Br = sqrt(1*pow(B[0],2) + 1*pow(B[1],2));
 				float Er = sqrt(1*pow(E[0],2) + 1*pow(E[1],2));
 				float k = 1*pow(e*mi,2) * (B[2]*Er - Br*E[2]);
 				float p = Er - B[2] * k / 1*pow(OMEGA,2);
 				float s = E[2] + Br * k / 1*pow(OMEGA,2);
 				float drift0 = -k  / (1*pow(OMEGA,2) + 1./ 1*pow(TAU,2));
 				float drift1 = (e / (mi * 1*pow(TAU,2))) * (1*pow(TAU, 3) * p - (B[2] * k / 1*pow(OMEGA,4)) * (-1*pow(TAU, 3) * 1*pow(OMEGA, 2) / (1. + 1*pow(TAU, 2) * 1*pow(OMEGA,2)) ));
 				float drift2 = (e / (1*pow(TAU, 2) * mi)) * (1*pow(TAU,3) * s + (Br * k / 1*pow(OMEGA,4) ) * (-1*pow(TAU, 3) * 1*pow(OMEGA, 3) / (1. + 1*pow(TAU, 2) * 1*pow(OMEGA, 2))));
                float r = sqrt(1*pow(pos[0],2) + 1*pow(pos[1],2));
                float theta = atan(pos[2]/pos[0]);
                vdrift[0] = abs((drift1) * r * sin(theta)) * 1*pow(-1,signbit(vdrift[0]));
 				vdrift[1] = abs((drift1) * r * cos(theta)) * 1*pow(-1, signbit(vdrift[1]));
 				vdrift[2] = abs((drift2) * 1*pow(-1, signbit(vdrift[2])));
 			}

 			else {
 				throw std::invalid_argument("Method does not exist, choose between factor and derivation");
 			}

 			std::vector<float> mach(3);
 			mach = Minusvectors(vdrift,vel, mach);
 			mach = Multiplyscalar(mach,1./vT,mach);
  			float machmag = Magnitude(mach,machmag);
 			float LAMBDA = sqrt(1. / (exp(-1*pow(machmag, 2) / 2.) * 1*pow(lambdadi, -2) + 1*pow(lambdade,-2) ));
 			float beta = abs(Zd*1*pow(e,2)/ (LAMBDA*Ti*kb) );
 			float u = Magnitude(vdrift,u);
 			u = u/vT;
 			float z = abs(Zd) * 1*pow(e,2) / (4 * PI * e0 * radd * Te * kb);
 			float temptau = Te / Ti ;
 			float coloumblog = 5;
 			float forcefactor = sqrt(2 * PI) * 1*pow(radd, 2) * ni0 * mi * 1*pow(vT,2) *
 				(
 				sqrt(PI/ 2) * erf(u / sqrt(2)) *
                (1 + 1*pow(u,2) + (1 - 1*pow(u,-2)) * (1 + 2 * z * temptau) + 4 * 1*pow(z,2) * 1*pow(temptau, 2) * 1*pow(u,-2) * log(
                coloumblog)) 
                +
                (1*pow(u,-1) * (1 + 2 * z * temptau + 1*pow(u,2) - 4 * 1*pow(z,2) * 1*pow(temptau, 2) * log(coloumblog)) 
                * exp(- 1*pow(u, 2) / 2.))
                ) / machmag ;
 			acctemp = Multiplyscalar(mach,forcefactor/md,acctemp);
 			return acctemp;
 		}
 		else{
 			acctemp = {0,0,0};
 			return acctemp; 
 		}
 	}


 	};


void Dust::set_values(bool SWITCH, float a, float b, float c, float d, float e, float OT, std::vector<float> f, std::vector<float> g, std::vector<float> h, std::vector<float> i, std::vector<float> j, std::vector<float> k, std::vector<float> l){
	Bswitch = SWITCH;
	md = a;
	radd = b;
	lambdaD = c;
	phia = d;
	initcharge = e;
	OMEGATAU = OT;
	initpos = f;
	initvel = g;
	initacc = h;
	pos = i;
	vel = j;
	acc = k;
	multifields = l;
};




int main(){
	Dust dust0, dust1;
	std::vector<float> positiond0(3);
	positiond0={-3*lambdaD,0,0.0003825734};
	std::vector<float> positiond1(3);
	positiond1={3*lambdaD,0,0.0003825734};
	std::vector<float> posd0 = positiond0;
	std::vector<float> posd1 = positiond1;
	std::vector<float> vel1 = {0,0,0};
	std::vector<float> acc1 = {0,0,0};
	std::vector<float> vel = vel1;
	std::vector<float> acc = acc1;
	std::vector<float> multifields = {0,0,0};
	dust0.set_values(true, md, radd, lambdaD, phia, Zd, OMEGATAU, posd0, vel1, acc1, positiond0, vel, acc, multifields);
	dust1.set_values(true, md, radd, lambdaD, phia, Zd, OMEGATAU, posd1, vel1, acc1, positiond1, vel, acc, multifields);
	std::vector<float> test(3);
	std::vector<float> Btest(3);
	Btest = {0,0,0.014};
	test=dust0.EXBacchybrid(Btest, "derivation");
	std::cout << test[0] << "," << test[1] << "," << test[2] << std::endl;



	return 0;
}














