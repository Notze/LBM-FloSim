/**
 * LBM-FloSim's main file
 *
 * author: Tom Nordloh (nordloh@uni-bremen.de)
 */

#include <iostream>
#include <CL/cl.hpp>
#include <fstream>
#include <math.h>
#include <vector>
 
int main(int argc, char* argv[]) {
	using namespace std;
	
	// debugging option
	string arg1;
	if (argc>1) arg1=argv[1];
	bool debug = (arg1 == "debug") ? true : false;
	
	// declare values for simulation
	int    q    = 9;	// number of possible streaming directions in the 2DQ9 model
	int    re   = 10;	// choose a Reynoldsnumber (in the real world)
	double uMax = .05;	// the maximum velocity in LB methods must be <<< 1 (Mach expansion)
	int    tau  = 1;	// for this particular example, tau = 1, which results in \nu

	double nu    = 1./3.*(tau-.5);
	double omega = 1./tau;
	double ly1   = nu*re/uMax; // the width of the channel must match the realWorld Reynoldsnumber Re

	// scenario: consider a cylnder with length lx and height ly
	int lx = 16; // quite arbitrary, since periodic bc are implemented below
	int ly = (int) round(ly1);

	int    h1      = ly-1; // effective width of slit after half-way walls are considered
	int    a       = 1;
	double h       = 1.;
	double uMaxNew = nu*re/ly;
	
	// get gravity in LB units
	double uAvg = 2./3.*uMaxNew;
	double g    = 12.*uAvg*nu/pow(h1, 2);


	// debug
	if (debug) {
		cout << "\nMonitoring variables now:";
		cout << "\n q: \t\t" << q;
		cout << "\n re: \t\t" << re;
		cout << "\n uMax: \t\t" << uMax;
		cout << "\n tau: \t\t" << tau;
		
		cout << "\n nu: \t\t" << nu;
		cout << "\n omega: \t" << omega;
		cout << "\n ly1: \t\t" << ly1;
		
		cout << "\n lx: \t\t" << lx;
		cout << "\n ly: \t\t" << ly;
		
		cout << "\n h1: \t\t" << h1;
		cout << "\n a: \t\t" << a;
		cout << "\n h: \t\t" << h;
		cout << "\n uMaxNew: \t" << uMaxNew;
		
		cout << "\n uAvg: \t\t" << uAvg;
		cout << "\n g: \t\t" << g;
		cout << "\n\n";
	}
	
	// on a D2Q9 lattice the lattice velocities are c = (0,0), (1,0), (0,1) , ...
	pair<int,int> c[] = {{0,0}, {1,0}, {0,1}, {-1,0}, {0,-1}, {1,1}, {-1,1}, {-1,-1}, {1,-1}};
	//NotIndex = [c.tolist().index((-c[i]).tolist()) for i in range(Q)]
	//# the lattice weights w_i are in the same order as the velocities in c
	//w = 1./36. * ones(Q); w[0] = 4./9.; w[1:5] = 1./9.

	//# velocity initial condition: u(x,y,t=0) = 0
	//u = numpy.zeros((2,ly,lx)) # u = (ux,uy) = (ux(x,y),uy(x,y))
	//# rho initial condition: rho(x,y,t=0) = 1
	//rho = numpy.ones((ly,lx))
	
	return 0;
}
 
int stuff(){
    // get all platforms (drivers)
    std::vector<cl::Platform> all_platforms;
    cl::Platform::get(&all_platforms);
    if(all_platforms.size()==0){
        std::cout<<" No platforms found. Check OpenCL installation!\n";
        exit(1);
    }
    cl::Platform default_platform=all_platforms[0];
    std::cout << "Using platform: "<<default_platform.getInfo<CL_PLATFORM_NAME>()<<"\n";
 
    // get default device of the default platform
    std::vector<cl::Device> all_devices;
    default_platform.getDevices(CL_DEVICE_TYPE_ALL, &all_devices);
    if(all_devices.size()==0){
        std::cout<<" No devices found. Check OpenCL installation!\n";
        exit(1);
    }
    cl::Device default_device=all_devices[0];
    std::cout<< "Using device: "<<default_device.getInfo<CL_DEVICE_NAME>()<<"\n";
 
    cl::Context context({default_device});
 
	// read kernel source code from file
	std::ifstream code("simple_add.cl"); // open source file with kernel code
	code.seekg(0, std::ios::end);		 // jump to EOF
	int length = code.tellg();			 // length equals position of EOF
	code.seekg(0, std::ios::beg);		 // jump to beginning
	char kernel_code[length];			 // create char-array for code
	code.read(kernel_code, length);		 // read the source
	code.close();						 // close input file stream

	// build kernel program
    cl::Program program(context, kernel_code, false);
	if (program.build() != CL_SUCCESS) { // build(all_devices) bringt auch nix
		std::cerr << "failed to build program" << std::endl;
	}
 
    // create buffers on the device
    cl::Buffer buffer_A(context,CL_MEM_READ_WRITE,sizeof(int)*10);
    cl::Buffer buffer_B(context,CL_MEM_READ_WRITE,sizeof(int)*10);
    cl::Buffer buffer_C(context,CL_MEM_READ_WRITE,sizeof(int)*10);
    
	// build the kernel
    cl::Kernel simple_add(program, "simple_add");
    simple_add.setArg(0, buffer_A);
    simple_add.setArg(1, buffer_B);
    simple_add.setArg(2, buffer_C);
    
	//create queue to which we will push commands for the device.
    cl::CommandQueue queue(context, default_device);
    
    
    // declare input arrays as parameters for kernel
    int A[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    int B[] = {0, 1, 2, 0, 1, 2, 0, 1, 2, 0};
	
    // write arrays A and B to the device
    queue.enqueueWriteBuffer(buffer_A,CL_TRUE,0,sizeof(int)*10,A);
    queue.enqueueWriteBuffer(buffer_B,CL_TRUE,0,sizeof(int)*10,B);
    
    // run the kernel
    queue.enqueueNDRangeKernel(simple_add, cl::NullRange, cl::NDRange(10), cl::NullRange);
	queue.finish(); // wait until operations on device are finished

	// read the results from device
	int C[10];
	queue.enqueueReadBuffer(buffer_C,CL_TRUE,0,sizeof(int)*10,C);
 	
 	// monitor the arrays
	std::cout<<" Array A: ";
    for(int i=0;i<10;i++){
        std::cout<<A[i]<<" ";
    }
    std::cout<<std::endl;
    
    std::cout<<" Array B: ";
    for(int i=0;i<10;i++){
        std::cout<<B[i]<<" ";
    }
    std::cout<<std::endl;
    
    std::cout<<" Array C: ";
    for(int i=0;i<10;i++){
        std::cout<<C[i]<<" ";
    }
    std::cout<<std::endl;
    
    return 0;
}
