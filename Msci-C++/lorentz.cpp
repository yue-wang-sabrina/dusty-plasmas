#include <iostream>  


class Dust {
	int mass, radius, lambdaD;
	
	public:
		void set_values(int,int,int);
		int area() {return mass*lambdaD;}

};

void Dust::set_values(int x, int y, int z){
	mass = x;
	radius = y;
	lambdaD = z;
}


int main() {     
	auto lambda = [](auto x){ return x; 
	};     

	std::cout << lambda("Hello generic lambda!\n");    

	for (int i=0; i<1000000000; i++) {
		Dust dust0;
		dust0.set_values(1,2,3);
	}

	std::cout << "Done" << "\n";
	 
	return 0; 
}
