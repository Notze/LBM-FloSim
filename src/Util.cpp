/**
 * Utility Class
 *
 * author: Tom Nordloh (nordloh@uni-bremen.de)
 */
#include "Util.hpp"
#include <vector>

using namespace std;

namespace util{

	/**
	 * Takes two arrays of the same length and builds the scalar product
	 *
	 * @param a: first input array
	 * @param b: second input array
	 * @return: the scalar product of a and b
	 */
	int scalarProduct(const vector<int> a, const vector<int> b){
		//check if arrays are of same length
		if (a.size() != b.size()) throw runtime_error("scalarProduct: length of input arrays doesn't match");

		int out = 0;
		for (int i=0; i<a.size(); i++){
			out += a[i]*b[i];
		}

		return out;
	}

}
