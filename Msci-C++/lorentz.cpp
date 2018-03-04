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

float omega = (1.60217662 * pow(10,-19)) * 0.014 / (39.948 * 1.66053904 * pow(10,-27));
float m = 39.948 * 1.66053904 * pow(10,-27);
float tau = 1*pow(10,-5);
float Babs = 0.014;
float Eabs = 1.0;
float dt = 1*pow(10,-7);
float kb = 1.38064852 * pow(10,-23);
float T = 310;
float vinit = sqrt(kb*T/m);
float charge = 1.60217662 * pow(10,-19);
float phia =  -9.0126016670216487;
float Zd =  -9388.3579633332938;
int iterations = 1*pow(10,7);
std::vector<float> pos = {0,0,0};
std::vector<float> vel = {0,vinit,0};
std::vector<float> acc = {0,0,0};
std::vector<float> pos1 = {0,0,0};
std::vector<float> vel1 = {0,vinit,0};
std::vector<float> acc1 = {0,0,0};
std::vector<float> positionx(iterations+1);
std::vector<float> positiony(iterations+1);
std::vector<float> positionz(iterations+1);
const double PI = std::atan(1.0)*4;



std::vector<float> MUltiplyscalar(std::vector<float> v1, float scalar, std::vector<float> result_v){
	std::transform(v1.begin(), v1.end(), result_v.begin(), std::bind1st(std::multiplies<float>(),scalar));
	return result_v;
};


double RANDD() {
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(0, 1);//uniform distribution between 0 and 1
	return dis(gen);
};


class Ion {
	float omega, m, tau, Babs, Eabs, dt, charge;
	std::vector<float> pos, vel, acc, pos1, vel1, acc1;

	public:
		void set_values(float,float,float,float,float, float, float, std::vector<float>,std::vector<float>,std::vector<float>,std::vector<float>,std::vector<float>,std::vector<float>);
		std::vector<float> constB() {std::vector<float> B = {0, 0.014, 0};
		return B;
		}
		std::vector<float> constE() {std::vector<float> E = {0, 0, 1};
		return E;
		}
		std::vector<float> crossproduct(std::vector<float> v1, std::vector<float> v2, std::vector<float> cross_p) {
		cross_p[0] = (v1[1]*v2[2]) - (v1[2]*v2[1]);
	 	cross_p[1] = -((v1[0]*v2[2]) - (v1[2]*v2[0]));
	 	cross_p[2] = (v1[0]*v2[1]) - (v1[1]*v2[0]);
	 	return cross_p;
		}
		std::vector<float> addvectors(std::vector<float> v1, std::vector<float> v2, std::vector<float> added_v){
		std::transform (v1.begin(), v1.end(), v2.begin(), added_v.begin(), std::plus<float>());
		return added_v;
		}
		std::vector<float> minusvectors(std::vector<float> v1, std::vector<float> v2, std::vector<float> minus_v){
		std::transform (v1.begin(), v1.end(), v2.begin(), minus_v.begin(), std::minus<float>());
		return minus_v;
		}	
		std::vector<float> multiplyscalar(std::vector<float> v1, float scalar, std::vector<float> result_v){
		std::transform(v1.begin(), v1.end(), result_v.begin(), std::bind1st(std::multiplies<float>(),scalar));
		return result_v;
		}
		void changevelocity(std::vector<float> newvel){
		vel1 = newvel;
		}
		std::vector<float> getvel(){
			return vel;
		}
		float dotproduct(std::vector<float> v1, std::vector<float> v2, float result_dot){
			result_dot = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
			return result_dot;
		}

		std::vector<float> elementwisemultiply(std::vector<float> v1, std::vector<float> v2, std::vector<float> output_v){
		output_v[0] = v1[0]*v2[0];
		output_v[1] = v1[1]*v2[1];
		output_v[2] = v1[2]*v2[2];
		return output_v;
		}

		float getEabs(){
			return Eabs;
		}

		float getBabs(){
			return Babs;
		}

		std::vector<float> updateRK4() {
		std::vector<float> fv1(3);
		std::vector<float> cross_p(3);
		std::vector<float> added_v(3);
		std::vector<float> result_v(3);
		fv1 = multiplyscalar(addvectors(crossproduct(vel1, constB(),cross_p),constE(), added_v),charge/m,result_v);
		std::vector<float> fy1(3);
		fy1 = vel1;
		std::vector<float> v1(3);	
		v1 = addvectors(vel1,multiplyscalar(fv1,0.5*dt,result_v),added_v);

		std::vector<float> fv2(3);	
		fv2 = multiplyscalar(addvectors(crossproduct(v1, constB(),cross_p),constE(), added_v),charge/m,result_v);
		std::vector<float> fy2(3);
		fy2 = v1;
		std::vector<float> v2(3);	
		v2 = addvectors(vel1,multiplyscalar(fv2,0.5*dt,result_v),added_v);	

		std::vector<float> fv3(3);
		fv3 = multiplyscalar(addvectors(crossproduct(v2, constB(),cross_p),constE(), added_v),charge/m,result_v);
		std::vector<float> fy3(3);
		fy3 = v2;
		std::vector<float> v3(3);	
		v3 = addvectors(vel1,multiplyscalar(fv3,dt,result_v),added_v);

		std::vector<float> fv4(3);

		fv4 = multiplyscalar(addvectors(crossproduct(v3, constB(),cross_p),constE(), added_v),charge/m,result_v);
		std::vector<float> fy4(3);
		fy4 = v3;

		std::vector<float> added1(3);
		std::vector<float> added2(3);
		std::vector<float> added3(3);
		std::vector<float> added4(3);
		std::vector<float> multiplyresult(3);
		std::vector<float> added5(3);
		std::vector<float> added6(3);
		std::vector<float> added7(3);
		std::vector<float> added1p(3);
		std::vector<float> added2p(3);
		std::vector<float> added3p(3);
		std::vector<float> added4p(3);
		std::vector<float> multiplyresultp(3);
		std::vector<float> added5p(3);
		std::vector<float> added6p(3);
		std::vector<float> added7p(3);

		vel = addvectors(multiplyscalar(addvectors(addvectors(addvectors(addvectors(fv3,fv3,added7),fv4,added1),addvectors(fv2,fv2,added6),added2),fv1,added3),dt/6,multiplyresult),vel1,added5);					
		pos = addvectors(multiplyscalar(addvectors(addvectors(addvectors(addvectors(fy3,fy3,added7p),fy4,added1p),addvectors(fy2,fy2,added6p),added2p),fy1,added3p),dt/6,multiplyresultp),pos1,added5p);					
		vel1 = vel;
		pos1 = pos;
		return pos;
		}
};


class IonAnalysis {
	int iterations;
	float tau, omega, dt, vinit;
	// std::vector<float> KEERK, KEEF, KEEB, drift, driftdistancecol;
	std::vector<float> positionx, positiony, positionz;
	Ion ion0;
	// std::string method;

public:
	void set_values(int, float, float, float, float, std::vector<float>,std::vector<float>,std::vector<float>);
	
	float Dotproduct(std::vector<float> v1, std::vector<float> v2, float result_dot){
			result_dot = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
			return result_dot;
	}

	std::vector<float> Crossproduct(std::vector<float> v1, std::vector<float> v2, std::vector<float> cross_p) {
	cross_p[0] = (v1[1]*v2[2]) - (v1[2]*v2[1]);
 	cross_p[1] = -((v1[0]*v2[2]) - (v1[2]*v2[0]));
 	cross_p[2] = (v1[0]*v2[1]) - (v1[1]*v2[0]);
 	return cross_p;
	}

	std::vector<float> Elementwisemultiply(std::vector<float> v1, std::vector<float> v2, std::vector<float> output_v){
	output_v[0] = v1[0]*v2[0];
	output_v[1] = v1[1]*v2[1];
	output_v[2] = v1[2]*v2[2];
	return output_v;
	}

	float Normalise(std::vector<float> v1, float output_float){
		output_float = sqrt(1*pow(v1[0],2) + 1*pow(v1[1],2) + 1*pow(v1[2],2));
		return output_float;
	}

	std::vector<float> Multiplyscalar(std::vector<float> v1, float scalar, std::vector<float> result_v){
	std::transform(v1.begin(), v1.end(), result_v.begin(), std::bind1st(std::multiplies<float>(),scalar));
	return result_v;
	}

	std::vector<float> ABS(std::vector<float> v1, std::vector<float> out_v){
		out_v = {abs(v1[0]), abs(v1[1]), abs(v1[2])};
		return out_v;
	}

	void runsim(){
		float vinit = sqrt(kb * T / m);
		ion0.set_values(omega, tau, Babs, Eabs, dt, charge, m, pos, vel, acc, pos1, vel1, acc1);

        for (int i=0; i<iterations; i++){
        std::vector<float> updated;
        updated = ion0.updateRK4();
        positionx[i] = updated[0];
        positiony[i] = updated[1];
		positionz[i] = updated[2];	
        }          

	}

	void savepositions(){
		std::ofstream testfile1;
		std::ofstream testfile2;
		std::ofstream testfile3;

  		testfile1.open ("testx.txt");
  		testfile2.open ("testy.txt");
  		testfile3.open ("testz.txt");

		for (int j=0; j<iterations; j++){
			testfile1 << positionx[j] << std::endl;
			testfile2 << positiony[j] << std::endl;
			testfile3 << positionz[j] << std::endl;
			// std::cout << positionx[j] << "," << positiony[j] << "," << positionz[j] << std::endl;
		}
		testfile1.close();
		testfile2.close();
		testfile3.close();

	}

	double randd() {
  	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(0, 1);//uniform distribution between 0 and 1
	return dis(gen);
	}

	std::vector<float> thermalkick() {
		int tkick = floor(tau / dt);
		ion0.set_values(omega, tau, Babs, Eabs, dt, charge, m, pos, vel, acc, pos1, vel1, acc1);
		std::vector<float> updated(3);
		int index = 0;
		std::vector<float> newvelocity(3);

		positionx[0] = pos1[0];
		positiony[0] = pos1[1];
		positionz[0] = pos1[2];

		for (int i=1; i<int(iterations/tkick)+2; i++){
			index = index + 1 ;
			float theta = PI*randd();
			float phi = 2*PI*randd();
			std::vector<float> rhat = {sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)};
			newvelocity[0] = vinit*rhat[0];
			newvelocity[1] = vinit*rhat[1];
			newvelocity[2] = vinit*rhat[2];
			ion0.changevelocity(newvelocity);
			updated = ion0.updateRK4();

			positionx[index] = updated[0];
	        positiony[index] = updated[1];
			positionz[index] = updated[2];

			for (int k=1; k< tkick; k++){
				index = index + 1;
				updated = ion0.updateRK4();
				positionx[index] = updated[0];
				positiony[index] = updated[1];
				positionz[index] = updated[2];
			}
		}

		std::vector<float> directionnonorm(3);
		float Norm;
		directionnonorm = Crossproduct(ion0.constE(),ion0.constB(), directionnonorm);
		Norm = Normalise(directionnonorm,Norm);
		std::vector<float> direction(3);
		direction = Multiplyscalar(directionnonorm,1./Norm,direction);
		direction = ABS(direction,direction);
		std::vector<float> drift(3);
		float driftabs;
		std::vector<float> diffdist;
		float finalposx;
		// float finalposy;
		// float finalposz;
		float initposx;
		// float initposy;
		// float initposz;
		finalposx = positionx[iterations];

		// finalposy = positiony.back();
		// finalposz = positionz.back();
		initposx = positionx[0];
		// initposy = positiony[0];
		// initposz = positionz[0];

		diffdist = {finalposx-initposx, 0,0}; //finalposy-initposy, finalposz-initposz};
		drift = Elementwisemultiply(direction,diffdist,drift);
		
		float driftnocol;
		driftnocol = (ion0.getEabs()/ion0.getBabs())*(iterations-1)*dt;
		std::vector<float> driftnocolvec;
		driftnocolvec = {-drift[0]/driftnocol, 0, 0};

		return driftnocolvec;

	}

	std::vector<float> thermalkickexponential() {
		int tkick = floor(tau / dt);
		ion0.set_values(omega, tau, Babs, Eabs, dt, charge, m, pos, vel, acc, pos1, vel1, acc1);
		std::vector<float> updated(3);
		int index = 0;
		std::vector<float> newvelocity(3);

		positionx[0] = pos1[0];
		positiony[0] = pos1[1];
		positionz[0] = pos1[2];

		for (int i=1; i<iterations+1; i++){
			index = index + 1 ;

			float probcol;
			probcol = 1.0 - 1*pow(exp(1),-dt/tau);
			float detcol;
			detcol=randd();
			
			if (detcol <= probcol){
			float theta = PI*randd();
			float phi = 2*PI*randd();
			std::vector<float> rhat = {sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)};
			newvelocity[0] = vinit*rhat[0];
			newvelocity[1] = vinit*rhat[1];
			newvelocity[2] = vinit*rhat[2];
			ion0.changevelocity(newvelocity);
			updated = ion0.updateRK4();
			}

			else{
			updated = ion0.updateRK4();
			}
		
			positionx[index] = updated[0];
	        positiony[index] = updated[1];
			positionz[index] = updated[2];
		}


		std::vector<float> directionnonorm(3);
		float Norm;
		directionnonorm = Crossproduct(ion0.constE(),ion0.constB(), directionnonorm);
		Norm = Normalise(directionnonorm,Norm);
		std::vector<float> direction(3);
		direction = Multiplyscalar(directionnonorm,1./Norm,direction);
		direction = ABS(direction,direction);
		std::vector<float> drift(3);
		float driftabs;
		std::vector<float> diffdist;
		float finalposx;
		// float finalposy;
		// float finalposz;
		float initposx;
		// float initposy;
		// float initposz;
		finalposx = positionx[iterations];

		// finalposy = positiony.back();
		// finalposz = positionz.back();
		initposx = positionx[0];
		// initposy = positiony[0];
		// initposz = positionz[0];

		diffdist = {finalposx-initposx, 0,0}; //finalposy-initposy, finalposz-initposz};
		drift = Elementwisemultiply(direction,diffdist,drift);
		
		float driftnocol;
		driftnocol = (ion0.getEabs()/ion0.getBabs())*(iterations-1)*dt;
		std::vector<float> driftnocolvec;
		driftnocolvec = {-drift[0]/driftnocol, 0, 0};

		return driftnocolvec;

	}


};


void Ion::set_values(float x, float y, float z, float h, float deltat, float Q, float mass, std::vector<float> a, std::vector<float> b, std::vector<float>c, std::vector<float>d, std::vector<float>e, std::vector<float>f){
	omega = x;
	tau = y;
	Babs = z;
	Eabs = h;
	dt = deltat;
	charge = Q;
	m = mass;
	pos = a;
	vel = b;
	acc = c;
	pos1 = d;
	vel1 = e;
	acc1 = f;
};


void IonAnalysis::set_values(int a, float b, float c, float d, float p, std::vector<float> e, std::vector<float> f, std::vector<float> g){
	iterations = a;
	tau = b;
	omega = c;
	dt = d;
	vinit = p;
	positionx = e;
	positiony = f;
	positionz = g;
};


int main(int argc, char* argv[]) {   

	int c_num =0;
	
	switch(c_num) {
	    case 1 : {
	    	// Test for one value of omega*tau with many trials
	    	std::cout << "Single tau value old thermal kick method" << std::endl;
		    IonAnalysis analysis;
			analysis.set_values(iterations, tau, omega, dt, vinit, positionx, positiony, positionz);
			std::vector<float> test;
			std::ofstream DRIFT;
		 	DRIFT.open ("DRIFTRATIOtauEm5.txt");
			for (int i=0; i<50; i++){
			test = analysis.thermalkick();
			DRIFT << test[0] << std::endl;
			}
			DRIFT.close();
			break;
		}
		case 2 : {
			//Test for changing omega*tau for fixed time
			float trials = 5;
			std::vector<float> tauvals(trials);	
			tauvals = {float(8.9*pow(10,-6)), float(9.5*pow(10,-6)), float(1*pow(10,-5)),float(1.5*pow(10,-5)),float(1.7*pow(10,-5))};
			std::vector<std::string> names(trials);
			names = {"2DRIFT1.txt", "2DRIFT2.txt","2DRIFT3.txt","2DRIFT4.txt", "2DRIFT5.txt"};
			std::vector<float> omegatauvals(trials);
			omegatauvals =  MUltiplyscalar(tauvals,omega,omegatauvals);
			for (int j=0; j < names.size(); j++){
				IonAnalysis analysis;
				analysis.set_values(iterations, tauvals[j], omega, dt, vinit, positionx, positiony, positionz);
				std::vector<float> test;
				std::ofstream DRIFT;
			 	DRIFT.open(names[j]);
				for (int i=0; i<20; i++){
					test = analysis.thermalkick();
					DRIFT << test[0] << std::endl;
					}
				DRIFT.close();
			}	
			break;
		}
		case 3: {
			// Test new thermal kick method with exponential probabilistic implementation for tau = 10^(-5) to check it works ok
			float probcol;
			probcol = 1.0 - 1*pow(exp(1),-dt/(1*pow(10,-5)));
			float detcol;
			for (int l = 0; l <1000; l++){
				detcol=RANDD();
				if (detcol <= probcol){
					std::cout << detcol << std::endl;
				}
				else {}
			}
			break;
		}
		
		case 4: {
			// Test single value of omega*tau using new thermal kick with exponential method
			IonAnalysis analysis;
			analysis.set_values(iterations, 1*pow(10,-5), omega, dt, vinit, positionx, positiony, positionz);
			std::vector<float> test;
			for (int i=0; i<10; i++){
				test = analysis.thermalkickexponential();
				std::cout << test[0] << std::endl;
				}
			break;	
		}
		case 5: {
			//Testing changing omega*tau using new thermal kick with exponential method 
			float trials = 8;
			std::vector<float> tauvals(trials);	
			tauvals = {float(1*pow(10,-6)), float(6*pow(10,-6)), float(1.1*pow(10,-5)), float(1.3*pow(10,-5)), float(1.4*pow(10,-5)), float(1.5*pow(10,-5)), float(1.6*pow(10,-5)), float(1.9*pow(10,-5))};
			std::vector<std::string> names(trials);
			names = {"2NEWDRIFT1.txt", "2NEWDRIFT2.txt","2NEWDRIFT3.txt","2NEWDRIFT4.txt", "2NEWDRIFT5.txt","2NEWDRIFT6.txt","2NEWDRIFT7.txt","2NEWDRIFT8.txt"};
			std::vector<float> omegatauvals(trials);
			omegatauvals =  MUltiplyscalar(tauvals,omega,omegatauvals);
			for (int j=0; j < names.size(); j++){
				IonAnalysis analysis;
				analysis.set_values(iterations, tauvals[j], omega, dt, vinit, positionx, positiony, positionz);
				std::vector<float> test;
				std::ofstream DRIFT;
			 	DRIFT.open(names[j]);
				for (int i=0; i<20; i++){
					test = analysis.thermalkickexponential();
					DRIFT << test[0] << std::endl;
				}
				DRIFT.close();
			}
			break;
			
		}
		}


	return 0; 
}
