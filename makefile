all:
	mkdir -p target
	g++ src/LBM-FloSim.cpp -lOpenCL -std=c++11 -o target/LBM-FloSim
	cp src/*.cl target/
clean:
	rm -rf target/
