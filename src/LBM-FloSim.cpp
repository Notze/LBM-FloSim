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
#include "Util.cpp"

using namespace std;

// velocities according to D2Q9
vector<pair<int,int>> c = {{0,0}, {1,0}, {0,1}, {-1,0}, {0,-1}, {1,1}, {-1,1}, {-1,-1}, {1,-1}};
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

// TODO NotIndex = [c.tolist().index((-c[i]).tolist()) for i in range(Q)]

// the lattice weights w_i are in the same order as the velocities in c
// w = 1./36. * ones(Q); w[0] = 4./9.; w[1:5] = 1./9. TODO: w[1:4] = 1./9. (?)
vector<double> w = {4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36.};

/**
 *
 */
vector<int> equilibrium(vector<int> rho, vector<int> u){

	//    cu   = 3.0 * dot(c,u.transpose(1,0,2))
	//int cu = 3 * util::scalarProduct(c, u);
	// TODO wat macht dot() an dieser Stelle?

	//    usqr = 3./2.*(u[0]**2+u[1]**2)
	double usqr = (3./2.)*(u[0]^2 + u[1]^2);

	//    feq = zeros((Q,ly,lx))
	vector<int> feq(lx*ly*q, 0);

	// TODO ich kapier gar nix! :D
	//    for i in range(Q):
	//        feq[i,:,:] = rho*w[i]*(1.+cu[i]+0.5*cu[i]**2-usqr)
	//    return feq

	return feq;
}

/**
 * Grafikkarten-Shit
 */
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
	std::ifstream code("SimpleAdd.cl"); // open source file with kernel code
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


/**
 *
 */
int main(int argc, char* argv[]) {

	// debugging option
	string arg1;
	if (argc>1) arg1=argv[1];
	bool debug = (arg1 == "debug") ? true : false;
	
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

	// velocity initial condition: u(x,y,t=0) = 0
	//u = numpy.zeros((2,ly,lx)) // u = (ux,uy) = (ux(x,y),uy(x,y))
	vector<int> u(lx*ly*2, 0);

	//# rho initial condition: rho(x,y,t=0) = 1
	vector<int> rho(lx*ly, 1);
	
	//TODO wat?
	//sumpop = lambda fin: sum(fin,axis=0)

	// Initial populations TODO #rho != #u
	vector<int> f = equilibrium(rho, u);

//# change in velocity, due to gravity
//deltaU = tau/rho*g # ...is already an array, due to rho
//"""
//#MAIN LOOP
//"""
//#
//#
//#
//#
//t = 0
//tend = 3000
//ueq = u.copy()
//while t < tend:
//    # equilibrium function
//    feq = equilibrium(rho,ueq)
//
//
//    # collision step
//    fprime = f.copy()
//    fprime = f - omega * (f - feq)
//    f, fprime = fprime, f
//
//
//    # streaming step: f_i(r+c_i, t+1) = f_i(r, t+1)
//    fstreamed = f.copy()
//    for ii in range(Q):
//        fstreamed[ii] = numpy.roll( numpy.roll(f[ii],c[ii,0],axis=1) ,-c[ii,1], axis=0 )
//
//    # boundaries: half way bounce back at walls
//    northWallY = 0
//    southWallY = ly-1
//
//
//    for ii in [2,5,6]:
//        # half way bounce back
//        fstreamed[ii,southWallY,:] = fstreamed[NotIndex[ii],southWallY,:].copy()
//        fstreamed[NotIndex[ii],northWallY,:] = fstreamed[ii,northWallY,:].copy()
//
//        # density! adjust complementary populations
//        fstreamed[NotIndex[ii],southWallY,:] = fstreamed[ii,southWallY,:].copy()
//        fstreamed[ii,northWallY,:] = fstreamed[NotIndex[ii],northWallY,:].copy()
//
//
//    # swap the arrays
//    f, fstreamed = fstreamed, f
//
//    #rho = numpy.sum(f,axis=0)
//    rho = sumpop(f)
//    u = dot(c.transpose(), f.transpose((1,0,2)))/rho
//    ueq = u + deltaU
//    uprofile = numpy.zeros(ly)
//    uprofile[:] = u[0,:,lx/2]
//
//    # make plots at the t instants
//    tInstants = numpy.rint(numpy.linspace(0,tend-1,10))
//    for plotTime in tInstants:
//		if abs(t-plotTime) < 10e-03:
//			plt.plot(uprofile)
//
//
//
//
//#plt.show()
//
//
//    if (t%100==0):
//            ufield = numpy.zeros((ly,lx))
//            ufield = sqrt(u[0]**2 + u[1]**2)
//            plt.clf();
//            cmap = plt.cm.rainbow
//            norm = plt.Normalize(ufield.min(), ufield.max())
//            rgba = cmap(norm(ufield))
//            #rgba[obstacle, :3] = 0, 0, 0
//            plt.imshow(rgba, interpolation='nearest')
//            #plt.show()
//            #plt.imshow(obstacle2, cmap = cm.Greys)
//            #plt.imshow(sqrt(u[0]**2+u[1]**2),cmap=cm.Reds)
//            plt.colorbar()
//            #plt.show()
//            plt.savefig("t."+str(t/100).zfill(4)+".png")
//    t += 1
//    print numpy.around(float(t)/tend * 100.,decimals=0), "% \r",
//
	return 0;
}

