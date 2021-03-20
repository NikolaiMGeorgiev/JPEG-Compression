#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <algorithm>  
#include "HuffmanTree.h"
#include "Image.h"
#include <boost\dynamic_bitset.hpp>

using namespace std;

void testImageCompression(std::string path) {
	Image img(path);

	img.compress();
	img.printTestResults();
	img.printBits();
	img.printTestOfFirstTriplet();

	cout << endl;
}

int main() {
	testImageCompression("img/H16xW8.jpg");
	testImageCompression("img/H8xW16.jpg");
	testImageCompression("img/H8xW8.jpg");
	testImageCompression("img/H16xW12.jpg");
	testImageCompression("img/H12xW16.jpg");
	testImageCompression("img/H12xW12.jpg");
	testImageCompression("img/H7xW4.jpg");
	testImageCompression("img/H4xW7.jpg");
	testImageCompression("img/H7xW7.jpg");
	testImageCompression("img/H32xW16.jpg");
	testImageCompression("img/H16xW32.jpg");
	testImageCompression("img/H32xW32.jpg");
	testImageCompression("img/H32xW20.jpg");
	testImageCompression("img/H20xW32.jpg");
	testImageCompression("img/H36xW36.jpg");
	testImageCompression("img/nonexistent_image.jpg");
	
	system("pause");
	return 0;
}