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
#include <sstream>
#include <map>
#include "progbar.h"

const double PI = std::atan(1.0)*4;
float radd = 1.5 * 1*pow(10,-6);
float md = 1*pow(radd,3)*PI*1000*4/3;
float omega = (1.60217662 * pow(10,-19)) * 0.014 / (39.948 * 1.66053904 * pow(10,-27));
float OMEGATAU = 0.01; // From Konopka experimentalbasis paper
float mi = 39.948 * 1.66053904 * pow(10,-27);
float me = 9.10938356 * pow(10,-31);
float gravity = -9.8;
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
float mu0 = 4 * PI * 1*pow(10,-7);
float magBmom = (2 * PI * 1*pow(0.003,3) * 0.014 / mu0);
std::vector<float> Bmomhat = {0,0,1};
std::vector<float> Bmom = {0,0, magBmom};
std::vector<float> dipolepos = {0, 0, -0.0005};
int iterations = 1*pow(10,5);
std::vector<float> positionx;
std::vector<float> positiony;
std::vector<float> positionz;
std::vector<float> multifields = {0,0,0};

int n_particles = 2;
int iterationsB = 10000;
int init_iterations = 100;
std::vector<float> initpositionsx(n_particles); 
std::vector<float> initpositionsy(n_particles);
std::vector<float> initpositionsz(n_particles); 
float t = 1;
std::vector<std::vector<int>> pairs;



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
};

std::vector<float> Addvectors(std::vector<float> v1, std::vector<float> v2, std::vector<float> minus_v){
std::transform (v1.begin(), v1.end(), v2.begin(), minus_v.begin(), std::plus<float>());
return minus_v;
};

std::vector<float> Crossproduct(std::vector<float> v1, std::vector<float> v2, std::vector<float> cross_p) {
cross_p[0] = (v1[1]*v2[2]) - (v1[2]*v2[1]);
cross_p[1] = -((v1[0]*v2[2]) - (v1[2]*v2[0]));
cross_p[2] = (v1[0]*v2[1]) - (v1[1]*v2[0]);
return cross_p;
};

std::string Convert (int number){
    std::ostringstream buff;
    buff<<number;
    return buff.str();   
};

float Dotproduct(std::vector<float> v1, std::vector<float> v2, float result_dot){
	result_dot = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
	return result_dot;
};


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

	void changevelocity(std::vector<float> newvel){
		vel = newvel;
	}

	void changeposition(std::vector<float> newpos){
		pos = newpos;
	}

	void changeacc(std::vector<float> newacc){
		acc = newacc;
	}

	void changeBswitch(bool validity){
		if (validity){
			Bswitch = true;
		}
		else{
			Bswitch = false;
		}
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
		acctemp = Multiplyscalar(vel,1./magvel,acctemp);
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
 			throw std::invalid_argument( "Dust fell out of cylinder! Origin: sheathfield()" + std::to_string(pos[0]) + "," + std::to_string(pos[1]) + "," + std::to_string(pos[2]));
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
            float r = sqrt(1*pow(pos[0],2) + 1*pow(pos[1],2));
            float theta;
            float Factor;

            if (float(abs(pos[0])) == float(0)){
            	if (pos[1] > float(0)){
            		theta = PI/2;
            	}
            	else if (pos[1] < float(0)){
            		theta = 3*PI/2;
            	}
            	else{
            		throw std::invalid_argument("1Problem with theta calculation EXB drift");
            	}
            }
            else if (float(abs(pos[1])) == float(0)){
            	if (pos[0] > float(0)){
            		theta = 0;
            	}
            	else if (pos[0] < float(0)){
            		theta = PI;
            	}
            	else{
            		throw std::invalid_argument("2Problem with theta calculation EXB drift");
            	}

            }
            else{
            	if (pos[0]>float(0) and pos[1]>float(0)){
            		theta = atan(pos[1]/pos[0]);
            	}
            	else if (pos[0]>float(0) and pos[1]<float(0)){
            		theta = 2*PI - atan(abs(pos[1] /pos[0]));
            	}
            	else if (pos[0]<float(0) and pos[1]>float(0)){
            		theta = PI - atan(abs(pos[1] / pos[0]));
            	}
            	else if (pos[0]<float(0) and pos[1]<float(0)){
                    theta = PI + atan(abs(pos[1] / pos[0]));
            	}
            	else{
            		throw std::invalid_argument("3Problem with theta calculation EXB drift");
            	}
            }
 			if (method == "factor"){
 				float factor = 1*pow(OMEGATAU,2)/(1. + 1*pow(OMEGATAU,2));
 				vdrift = Multiplyscalar(vdrift,factor,vdrift);

 			}

 			else if (method == "derivation"){
 				float Br = sqrt(1*pow(B[0],2) + 1*pow(B[1],2));
 				float Er = sqrt(1*pow(E[0],2) + 1*pow(E[1],2));
 				// float k = 1*pow(e*mi,2) * (B[2]*Er - Br*E[2]);
 				// float p = Er - B[2] * k / 1*pow(OMEGA,2);
 				// float s = E[2] + Br * k / 1*pow(OMEGA,2);
 				// float drift0 = -k  / (1*pow(OMEGA,2) + 1./ 1*pow(TAU,2));
 				// float drift1 = (e / (mi * 1*pow(TAU,2))) * (1*pow(TAU, 3) * p - (B[2] * k / 1*pow(OMEGA,4)) * (-1*pow(TAU, 3) * 1*pow(OMEGA, 2) / (1. + 1*pow(TAU, 2) * 1*pow(OMEGA,2)) ));
 				// float drift2 = (e / (1*pow(TAU, 2) * mi)) * (1*pow(TAU,3) * s + (Br * k / 1*pow(OMEGA,4) ) * (-1*pow(TAU, 3) * 1*pow(OMEGA, 3) / (1. + 1*pow(TAU, 2) * 1*pow(OMEGA, 2))));
                Factor = (1*pow(e,2) * 1*pow(TAU,2) /(1*pow(mi,2) * (1+1*pow(OMEGATAU,2)))) * (Er * B[2] - E[2]*Br);
                vdrift[0] = Factor * -1*sin(theta);
 				vdrift[1] = Factor * cos(theta);
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
 		}
 		else{
 			acctemp = {0,0,0};
 		}

 		return acctemp; 
 	}

 	std::vector<float> dipoleB(std::vector<float> r){
 		std::vector<float> diffr(3);
 		diffr = Minusvectors(pos,r,diffr);
 		std::vector<float> Bfieldtemp(3);
 		if (Bswitch && pos[2] < sheathd) {
 			float magr = Magnitude(diffr,magr);
 			std::vector<float> rhat(3);
 			rhat = Multiplyscalar(diffr,1./magr,rhat);
 			float Bdotrhat;
 			Bdotrhat = Dotproduct(rhat, Bmom, Bdotrhat);
 			std::vector<float> temp1(3);
 			temp1 = Multiplyscalar(rhat,3*Bdotrhat,temp1);
 			std::vector<float> temp2(3);
 			temp2 = Minusvectors(temp1, Bmom,temp2);
 			Bfieldtemp = Multiplyscalar(temp2,mu0/(4*PI*1*pow(magr,3)),Bfieldtemp);
 		}
 		else{
 			Bfieldtemp = {0,0,0};
 		}
 		return Bfieldtemp;
 	}

 	std::vector<float> combinedrift(std::vector<float> B){
 		std::vector<float> totaldrift(3);

 		if (Bswitch){
 			float magB = Magnitude(B,magB);
 			float vperp = sqrt(2*kb*Ti/mi);
 			float vpa = vperp*sqrt(1./2.);
 			float r = Magnitude(pos,r);
 			float theta;
 			if (float(abs(pos[0])) == float(0)){
            	if (pos[1] > float(0)){
            		theta = PI/2;
            	}
            	else if (pos[1] < float(0)){
            		theta = 3*PI/2;
            	}
            	else{
            		throw std::invalid_argument("1Problem with theta calculation EXB drift");
            	}
            }
            else if (float(abs(pos[1])) == float(0)){
            	if (pos[0] > float(0)){
            		theta = 0;
            	}
            	else if (pos[0] < float(0)){
            		theta = PI;
            	}
            	else{
            		throw std::invalid_argument("2Problem with theta calculation EXB drift");
            	}

            }
            else{
            	if (pos[0]>float(0) and pos[1]>float(0)){
            		theta = atan(pos[1]/pos[0]);
            	}
            	else if (pos[0]>float(0) and pos[1]<float(0)){
            		theta = 2*PI - atan(abs(pos[1] /pos[0]));
            	}
            	else if (pos[0]<float(0) and pos[1]>float(0)){
            		theta = PI - atan(abs(pos[1] / pos[0]));
            	}
            	else if (pos[0]<float(0) and pos[1]<float(0)){
                    theta = PI + atan(abs(pos[1] / pos[0]));
            	}
            	else{
            		throw std::invalid_argument("3Problem with theta calculation EXB drift");
            	}
            }
            std::vector<float> gradB(3);
            std::vector<float> tempr(3);
            tempr = {cos(theta), sin(theta),0};
            float temp3 = (-3. * mu0 * sqrt(9 * 1*pow(pos[2],2) + 1*pow(magBmom, 2)) / (4 * PI * 1*pow(r,4)));
            std::vector<float> temp4(3);
            temp4 = Multiplyscalar(tempr,temp3,temp4);
            float temp5 = mu0 * 18 * pos[2] / (8 * PI * 1*pow(r, 3) * sqrt(9 * 1*pow(pos[2],2) + 1*pow(magBmom,2)));
            std::vector<float> temp6(3);
            std::vector<float> Bdirtemp = {0,0,1};
            temp6 = Multiplyscalar(Bdirtemp,temp5,temp6);
            gradB = Addvectors(temp4,temp6,gradB);
            float temp7 = (mi / ( e * magB )) * (1*pow(vpa,2) + 0.5 * 1 * pow(vperp,2)) / (1 * pow(magB,2));
            totaldrift = Crossproduct(B, gradB, totaldrift);
            totaldrift = Multiplyscalar(totaldrift,temp7,totaldrift);
 		}
 		else{
 			totaldrift = {0,0,0};
 		}
 		return totaldrift;
 	}

 	void updateEuler(){
 		acc = {0,0,gravity};
		std::vector<float> accsheath(3);
		std::vector<float> accradial(3);
		std::vector<float> accB(3);
		accsheath = sheathfield();
		accsheath = Multiplyscalar(accsheath, abs(initcharge/md),accsheath);
		accradial = radialfield();
		accradial = Multiplyscalar(accradial, abs(initcharge/md),accradial);
		accB = EXBacchybrid(dipoleB(dipolepos),"factor");
		acc = Addvectors(acc,accsheath,acc);
		std::vector<float> temp(3);
		temp = Multiplyscalar(multifields, abs(initcharge/md),temp);
		acc = Addvectors(acc,temp,acc);
		std::vector<float> damp(3);
		damp = damping();
		acc = Addvectors(acc,damp,acc);
		acc = Addvectors(acc,accradial,acc);
		acc = Addvectors(acc, accB, acc);
		std::vector<float> temp2(3);
		temp2 = Multiplyscalar(acc,dt,temp2);
		vel = Addvectors(vel,temp2,vel);
		std::vector<float> temp3(3);
		temp3 = Multiplyscalar(vel,dt,temp3);
		std::vector<float> temp4(3);
		temp4 = Multiplyscalar(acc,0.5*1*pow(dt,2),temp4);
		pos = Addvectors(pos,temp3,pos);
		pos = Addvectors(pos,temp4,pos);
		multifields = {0,0,0};
	}

};



class DustAnalysis{
	std::vector<Dust> dustlist;
	std::vector<std::vector<int>> pairs;
	std::vector<float> initx, inity, initz, positionsx, positionsy, positionsz;
	float t;
	int n_particles, iterationsB, init_iterations;
	
	public:
	void set_values(std::vector<Dust>, std::vector<std::vector<int>>,std::vector<float>,std::vector<float>,std::vector<float>,std::vector<float>, std::vector<float>,std::vector<float>,float, int, int, int);

	std::vector<Dust> getdustnames(){
		return dustlist;
	}

	std::vector<float> getallpositionx(){
		return positionx;
	}

	std::vector<float> getallpositiony(){
		return positiony;
	}

	std::vector<float> getallpositionz(){
		return positionz;
	}	

	void create_particles(){
		std::vector<std::string> names(n_particles);

		for (int int_name=0; int_name < n_particles; int_name ++){
			Dust dust;
			std::vector<float> initzeros(3);
			initzeros = {0,0,0};
			std::vector<float> INITPOS(3);
			INITPOS[0] = initx[int_name];
			INITPOS[1] = inity[int_name];
			INITPOS[2] = initz[int_name];
			dust.set_values(false, md, radd, lambdaD, phia, Zd*e, OMEGATAU, INITPOS, initzeros, initzeros, INITPOS, initzeros, initzeros, multifields);
			dustlist[int_name] = dust;

		}
	}

	void create_pairs(){
		for (int i=0; i<n_particles; i++){
			for (int j=i+1; j<n_particles;j++){
				std::vector<int> pair;
				pair = {i,j};
				pairs.push_back(pair);
			}
		}
	}

	std::vector<std::vector<int>> getselfpairs(){
		return pairs;
	}

	void interact_and_iterate(){
		progbar bar(std::cerr, 36, '=', ' ');

		for (int itone = 0; itone < init_iterations; itone ++){
			std::vector<std::vector<int>> pairsfinal;
			for (int b = 0; b < pairs.size(); b++) {
				if (dustlist[pairs[b][0]].intergraind(dustlist[pairs[b][1]])){
					pairsfinal.push_back(pairs[b]);
				}
			}
			for (int k = 0; k<pairsfinal.size(); k++){
				std::vector<float> interactfield(3);
				interactfield = dustlist[pairsfinal[k][0]].selffield(dustlist[pairsfinal[k][1]]);
				dustlist[pairsfinal[k][0]].selffieldmany(interactfield);
				dustlist[pairsfinal[k][1]].selffieldmany(Multiplyscalar(interactfield,-1,interactfield));
			}
			for (int D = 0; D < dustlist.size(); D++){
				dustlist[D].updateEuler();
				std::vector<float> selfpos(3);
				selfpos = dustlist[D].getselfpos();
				positionx.push_back(selfpos[0]);
				positiony.push_back(selfpos[1]);
				positionz.push_back(selfpos[2]);
				}
		}
			

		for (int l=0; l<dustlist.size(); l++){
			dustlist[l].changeBswitch(true);
		}

		for (int ittwo = 0; ittwo < iterationsB; ittwo ++ ){
			std::vector<std::vector<int>> pairsfinal;
			for (int b = 0; b < pairs.size(); b++) {
				if (dustlist[pairs[b][0]].intergraind(dustlist[pairs[b][1]])){
					pairsfinal.push_back(pairs[b]);
				}		
			}
			for (int k = 0; k<pairsfinal.size(); k++){
				std::vector<float> interactfield(3);
				interactfield = dustlist[pairsfinal[k][0]].selffield(dustlist[pairsfinal[k][1]]);
				dustlist[pairsfinal[k][0]].selffieldmany(interactfield);
				dustlist[pairsfinal[k][1]].selffieldmany(Multiplyscalar(interactfield,-1,interactfield));
			}
			for (int D=0; D<dustlist.size(); D++){
				dustlist[D].updateEuler();
				std::vector<float> selfpos(3);
				selfpos = dustlist[D].getselfpos();
				positionx.push_back(selfpos[0]);
				positiony.push_back(selfpos[1]);
				positionz.push_back(selfpos[2]);	
			}

			if (ittwo % 1000 == 0) {
				bar.update(ittwo, iterationsB);
			}

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

void DustAnalysis::set_values(std::vector<Dust> a, std::vector<std::vector<int>> b, std::vector<float> inx, std::vector<float> iny, std::vector<float> inz, std::vector<float> c, std::vector<float> y ,std::vector<float> z,float d, int f, int g, int h){
	dustlist = a;
	pairs = b;
	initx = inx;
	inity = iny;
	initz = inz;
	positionsx = c;
	positionsy = y;
	positionsz = z;
	t = d;
	n_particles = f;
	iterationsB = g;
	init_iterations = h;
}




int main(){
	// Testing dust functions work

	std::vector<Dust> dustlist(n_particles);

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
//	DustAnalysis analysis;
//	initpositionsx = {3*lambdaD,-3*lambdaD};
//	initpositionsy = {0,0};
//	initpositionsz = {0.0003825734,0.0003825734};
//	analysis.set_values(dustlist, pairs, initpositionsx, initpositionsy, initpositionsz, positionx, positiony, positionz,t, n_particles, iterationsB, init_iterations);
//	analysis.create_particles();
//	std::vector<Dust> testname(n_particles);
//	testname = analysis.getdustnames();
//	analysis.create_pairs();
//	analysis.interact_and_iterate();
//	std::vector<float> test;
//	test = analysis.getallpositionx();
//	 for (int i=0; i<positionx.size();i++){
//	 	std::cout << test[i] << "," << std::endl;
//	 }

    // std::vector<float> test;
	// dust0.updateEuler();
	// test = dust0.getselfacc();
	// std::cout << test[0] << "," << test[1] << "," << test[2] << std::endl;

	return 0;
}














