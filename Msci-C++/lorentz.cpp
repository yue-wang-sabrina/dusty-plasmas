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
std::vector<float> pos = {0,0,0};
std::vector<float> vel = {vinit,0,0};
std::vector<float> acc = {0,0,0};
std::vector<float> pos1 = {0,0,0};
std::vector<float> vel1 = {vinit,0,0};
std::vector<float> acc1 = {0,0,0};
std::vector<float> positionx;
std::vector<float> positiony;
std::vector<float> positionz;
float iterations = 1000;
const double PI = std::atan(1.0)*4;

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
		vel = newvel;
		}
		std::vector<float> getvel(){
			return vel;
		}
		float dotproduct(std::vector<float> v1, std::vector<float> v2, float result_dot){
			result_dot = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
			return result_dot;
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
	float tau, omega, dt;
	// std::vector<float> KEERK, KEEF, KEEB, drift, driftdistancecol;
	std::vector<float> positionx , positiony, positionz ;
	Ion ion0;
	// std::string method;

public:
	void set_values(int,float,float,float, std::vector<float>,std::vector<float>,std::vector<float>);

	void runsim(){
		float vinit = sqrt(kb * T / m);
		ion0.set_values(omega, tau, Babs, Eabs, dt, charge, m, pos, vel, acc, pos1, vel1, acc1);

        for (int i=0; i<iterations; i++){
        std::vector<float> updated;
        updated = ion0.updateRK4();
        positionx.push_back(updated[0]);
        positiony.push_back(updated[1]);
		positionz.push_back(updated[2]);	
  		
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
  	return (double)rand() / ((double)RAND_MAX + 1);
	}

	float thermalkick(){
		int tkick = floor(tau / dt);
		ion0.set_values(omega, tau, Babs, Eabs, dt, charge, m, pos, vel, acc, pos1, vel1, acc1);
		for (int i=0; i<floor(iterations/tkick); i++){
			float theta = PI*randd();
			float phi = 2*PI*randd();
			std::vector<float> rhat = {sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)};
			std::vector<float> newvelocity;
			newvelocity = multiplyscalar(rhat,vinit)
			ion0.changevelocity(newvelocity);
			ion0.updateRK4();
			for (k=0; k< tkick-1; k++){
				ion0.updateRK4();
			}
		std::vector<float> direction;
		direction = crossproduct(ion0.constE(),ion0.constB());
		std::vector<float> drift;
		float driftabs;
		std::vector<float>diffdist;
		std::vector<float> finalposx;
		std::vector<float> finalposy;
		std::vector<float> finalposz;
		std::vector<float> initposx;
		std::vector<float> initposy;
		std::vector<float> initposz;
		finalposx = positionx.back();
		finalposy = positiony.back();
		finalposz = positionz.back();
		initposx = positionx[0];
		initposy = positiony[0];
		initposz = positionz[0];

		std::vector<float> diffdist = {finalposx-initposx, finalposy-initposy, finalposz-initposz};
		driftabs = dotproduct(direction, diffdist, driftabs);
		drift = multiplyscalar(direction,diffdist);
		return drift;

		}

        // for i in tqdm(numpy.arange(int(iterations / tkick))):
        //     theta = math.pi * numpy.random.random_sample()
        //     phi = 2 * math.pi * numpy.random.random_sample()
        //     rhat = numpy.array([numpy.sin(theta) * numpy.cos(phi), numpy.sin(theta) * numpy.sin(phi), numpy.cos(theta)])
        //     ion1.vel1 = rhat * numpy.sqrt(const.kb * const.Ti / const.mi)
        //     ion1.updateRK4(B=ion1.constB())
        //     position.append(list(ion1.getselfpos()))
        //     velocity.append(numpy.linalg.norm(ion1.getselfvel()))
        //     for j in numpy.arange(tkick - 1):
        //         ion1.updateRK4(B=ion1.constB())
        //         position.append(list(ion1.getselfpos()))
        //         velocity.append(numpy.linalg.norm(ion1.getselfvel()))
        // direction = numpy.cross(ion1.constE(), ion1.constB())
        // drift = numpy.dot(numpy.array(position[-1]) - numpy.array(position[0]), direction)
        // self.drift = drift
        // self.position = position

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
}

void IonAnalysis::set_values(int a, float b, float c, float d, std::vector<float> e, std::vector<float> f, std::vector<float> g){
	iterations = a;
	tau = b;
	omega = c;
	dt = d;
	positionx = e;
	positiony = f;
	positionz = g;

}


int main() {   
	IonAnalysis analysis;
	analysis.set_values(iterations, tau, omega, dt, positionx, positiony, positionz);
	analysis.runsim();
	analysis.savepositions();
	analysis.thermalkick();



	//auto lambda = [](auto x){ return x; 
	//};     
	//std::cout << lambda("Hello generic lambda!\n");  

	// Ion ion0;
	// ion0.set_values(omega, tau, Babs, Eabs, dt, charge, m, pos, vel, acc, pos1, vel1, acc1);

	// for (int i=0; i<100; i++) {
		// std::vector<float> test(3);
		// test = ion0.updateRK4();
		// std::vector<float> cross = {0,1,0};
		// std::transform (cross.begin(), cross.end(), cross.begin(), cross.begin(), std::plus<float>());
		// std::cout << test[0] <<","<< test[1]<<","<< test[2] << std::endl;
		// std::cout<< v1[0] << v1[1] <<v1[2] << std::endl;
		// std::cout << cross < std::endl;
		// std::cout<< ion0.getomega() << std::endl;
	// }


	return 0; 
}
