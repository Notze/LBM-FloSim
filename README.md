# LBM-FloSim
The Lattice-Boltzmann-Method Flow Simulator (LBM-FlowSim or short FlowSim) is a manipulatable real-time simulation of liquids according to the Lattice-Boltzmann-Methods.

* supports only 2D simulations
* calculation takes place on GPU

# requirements
* a working OpenCL installation
* a c++ compiler
* make

# manual compiling instructions
* download the source and step into the project directory
* _make_
* the executable _LBM-FloSim_ and the kernel source files are located in the _target_ folder

# setting up Eclipse IDE
* create a new project
  * File -> New -> Makefile Project With Existing Code (don't forget to check a toolchain for automatic includes)
* create run configurations
  * Run -> Profile Configurations... -> C/C++ Application -> New 
    * in the "Main"-tab: set "C/C++ Application:" to "target/LBM-FloSim" under 
    * in the "Arguments"-tab: add "debug" to "Program arguments:" (if you want to run the program with debug flag)
    * in the "Arguments"-tab: change "Working directory:" from "${workspace_loc:LBM-FloSim}" to "${workspace_loc:LBM-FloSim}/target"
* enable syntax highlighting for kernel source files
  * for OpenCL syntax highlighting add OpenCLKernel.hpp to the preporcessors include path
    * get OpenCLKernel.hpp from [here](https://gist.github.com/rjeschke/5701109)
    * Project Poperties -> C/C++ General -> Preprocessor Include Paths, Macros etc. -> Entries -> GNU C++ -> Add...
    * set content type for cl-files to c++ source syntax
      * Window -> General -> Content Types -> C Source File -> C++ Source File -> Add "*.cl"
* in case you have problems with ambiguous syntax errors while "using namespace std" just deactivate the ambiguous warnings
  * Project -> Properties -> C/C++ General -> Code Analysis -> Configure Workspace Settings -> uncheck "Ambiguous problem"
