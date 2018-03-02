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
	
	// get all platforms (drivers)
	vector<cl::Platform> all_platforms;
	cl::Platform::get(&all_platforms);
	if(all_platforms.size()==0){
		cout<<" No platforms found. Check OpenCL installation!\n";
		exit(1);
	}
	cl::Platform default_platform=all_platforms[0];
	cout << "Using platform: "<<default_platform.getInfo<CL_PLATFORM_NAME>()<<"\n";
 
	// get default device of the default platform
	vector<cl::Device> all_devices;
	default_platform.getDevices(CL_DEVICE_TYPE_ALL, &all_devices);
	if(all_devices.size()==0){
		cout<<" No devices found. Check OpenCL installation!\n";
		exit(1);
	}
	cl::Device default_device=all_devices[0];
	cout<< "Using device: "<<default_device.getInfo<CL_DEVICE_NAME>()<<"\n";
 
	cl::Context context({default_device});
 
	// read kernel source code from file
	ifstream code("simple_add.cl"); // open source file with kernel code
	code.seekg(0, ios::end);		 // jump to EOF
	int length = code.tellg();			 // length equals position of EOF
	code.seekg(0, ios::beg);		 // jump to beginning
	char kernel_code[length];			 // create char-array for code
	code.read(kernel_code, length);		 // read the source
	code.close();						 // close input file stream

	// build kernel program
	cl::Program program(context, kernel_code, false);
	if (program.build() != CL_SUCCESS) { // build(all_devices) bringt auch nix
		cerr << "failed to build program" << endl;
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
	cout<<" Array A: ";
	for(int i=0;i<10;i++){
		cout<<A[i]<<" ";
	}
	cout<<endl;

	cout<<" Array B: ";
	for(int i=0;i<10;i++){
		cout<<B[i]<<" ";
	}
	cout<<endl;

	cout<<" Array C: ";
	for(int i=0;i<10;i++){
		cout<<C[i]<<" ";
	}
	cout<<endl;

	return 0;
}
