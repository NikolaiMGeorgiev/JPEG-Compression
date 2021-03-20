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

extern "C" {
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
}

bool load_image(std::vector<unsigned char>& image, const std::string& filename, int& x, int&y)
{
	int n;
	unsigned char* data = stbi_load(filename.c_str(), &x, &y, &n, 3);

	if (data != nullptr) {
		image = std::vector<unsigned char>(data, data + x * y * 3);
	}

	stbi_image_free(data);
	return (data != nullptr);
}

Image::Image(string filename) {
	int width, height = 0;
	vector<unsigned char> image;
	bool success = load_image(image, filename, width, height);

	if (!success) {
		std::cout << "Error loading image\n";
		return;
	}

	//cut the pixel in the blocks with size less than 8x8
	this->height = height - height % 8;
	this->width = width - width % 8;

	vector<vector<double>> red(height, vector<double>(width, 0));
	vector<vector<double>> green(height, vector<double>(width, 0));
	vector<vector<double>> blue(height, vector<double>(width, 0));

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			//split the vector of chars produced by load_image() into 3 matrices with double values for 
			//the red, blue and red components of the pixel's color
			red[i][j] = static_cast<double>(image[(i * width + j) * 3]);
			green[i][j] = static_cast<double>(image[(i * width + j) * 3 + 1]);
			blue[i][j] = static_cast<double>(image[(i * width + j) * 3 + 2]);
		}
	}

	colorTransformAndInit(red, green, blue);
}

void Image::colorTransformAndInit(vector<vector<double>>& red, vector<vector<double>>& green,
	vector<vector<double>>& blue) {
	const float gamma = 1;

	vector<vector<double>> y(height, vector<double>(width, 0));
	vector<vector<double>> cr(height, vector<double>(width, 0));
	vector<vector<double>> cb(height, vector<double>(width, 0));

	for (int i = 0; i < this->height; i++) {
		for (int j = 0; j < this->width; j++)
		{
			//center around zero by subtracting 128; new values will be within [-128,127]
			y[i][j] = 0.299 * pow(red[i][j], gamma) + 0.587 * pow(green[i][j], gamma)
				+ 0.114 * pow(blue[i][j], gamma) - 128;
			cr[i][j] = 128 - 0.168736 * pow(red[i][j], gamma) - 0.331264 * pow(green[i][j], gamma)
				+ 0.5 * pow(blue[i][j], gamma) - 128;
			cb[i][j] = 128 + 0.5 * pow(red[i][j], gamma) - 0.418688 * pow(green[i][j], gamma)
				- 0.081312 * pow(blue[i][j], gamma) - 128;
		}
	}

	this->y = y;
	this->cr = cr;
	this->cb = cb;
}

bool comparePairs(pair<int, int> p1, pair<int, int> p2) {
	return p1.second < p2.second;
}

void Image::compress() {
	if (this->y.empty() || this->cr.empty() || this->cb.empty()) {
		return;
	}

	vector<vector<double>> yNew = this->y;
	vector<vector<double>> crNew = this->cr;
	vector<vector<double>> cbNew = this->cb;
	unordered_map<int, int> freq;

	//apply cos transform and quantization for every 8x8 block of the image
	for (int i = 0; i < this->height; i += 8) {
		for (int j = 0; j < this->width; j += 8) {
			cosTransform(i, j, yNew, crNew, cbNew);
			quantization(i, j);
		}
	}

	initFrequencies(this->y, freq);
	initFrequencies(this->cr, freq);
	initFrequencies(this->cb, freq);

	for (auto f : freq) {
		pair<int, int> p = make_pair(f.first, f.second);
		this->freq.push_back(p);
	}

	sort(this->freq.begin(), this->freq.end(), comparePairs);

	encodeHuffman(0, 0);
}

void Image::cosTransform(int x, int y, vector<vector<double>> yNew, vector<vector<double>> crNew,
	vector<vector<double>> cbNew) {
	vector<double> result;

	for (int i = x; i < x + 7; i++) {
		for (int j = y; j < y + 7; j++) {
			result = g(x, y, i - x, j - y, yNew, crNew, cbNew);
			this->y[i][j] = result[0];
			this->cr[i][j] = result[1];
			this->cb[i][j] = result[2];
		}
	}
}

vector<double> Image::g(int x, int y, int u, int v, vector<vector<double>> yNew, vector<vector<double>> crNew,
	vector<vector<double>> cbNew) {
	const double pi = 2 * asin(1.0);
	vector<double> result(3, 0);

	for (int i = 0; i < 8; i++) {
		for (int j = 0; j < 8; j++) {
			result[0] += yNew[x + i][y + j] * cos((2 * i + 1) * u * (pi / 16)) * cos((2 * j + 1) * v * (pi / 16));
			result[1] += crNew[x + i][y + j] * cos((2 * i + 1) * u * pi / 16) * cos((2 * j + 1) * v * pi / 16);
			result[2] += cbNew[x + i][y + j] * cos((2 * i + 1) * u * pi / 16) * cos((2 * j + 1) * v * pi / 16);
		}
	}

	for (int i = 0; i < 3; i++) {
		if (u == 0) {
			result[i] *= 1 / sqrt(2);
		}

		if (v == 0) {
			result[i] *= 1 / sqrt(2);
		}

		result[i] *= 0.25;
	}

	return result;
}

void Image::quantization(int x, int y) {
	vector<vector<short>> qMatrix = {
		{16, 11, 10, 16, 24, 40, 51, 61},
		{12, 12, 14, 19, 26, 58, 60, 55},
		{14, 13, 16, 24, 40, 57, 69, 56},
		{14, 17, 22, 29, 51, 87, 80, 62},
		{18, 22, 37, 56, 68, 109, 103, 77},
		{24, 35, 55, 64, 81, 104, 113, 92},
		{49, 64, 78, 87, 103, 121, 120, 101},
		{72, 92, 95, 98, 112, 100, 103, 99}
	};

	for (int i = x; i < x + 7; i++) {
		for (int j = y; j < y + 7; j++) {
			this->y[i][j] = round(this->y[i][j] / qMatrix[i - x][j - y]);
			this->cr[i][j] = round(this->cr[i][j] / qMatrix[i - x][j - y]);
			this->cb[i][j] = round(this->cb[i][j] / qMatrix[i - x][j - y]);
		}
	}
}

void Image::initFrequencies(const vector<vector<double>>& px, unordered_map<int, int>& freq) {
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			if (freq.find(px[i][j]) == freq.end()) {
				freq.insert(make_pair(px[i][j], 1));
			}
			else {
				freq[px[i][j]]++;
			}
		}
	}
}

void Image::encodeHuffman(int x, int y) {
	HuffmanTree tree(freq);
	vector<bool> bits;
	unordered_map<int, vector<bool>> table;
	tree.getBitRepresentations(bits, table, tree.root);

	this->bitsTable = table;
	this->bits = getBitRepresentation(table);
}

boost::dynamic_bitset<> Image::getBitRepresentation(unordered_map<int, vector<bool>>& table) {
	int size = 0;

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			size += table.at((int)y[i][j]).size();
			size += table.at((int)cb[i][j]).size();
			size += table.at((int)cr[i][j]).size();
		}
	}

	boost::dynamic_bitset<> bits(size);
	boost::dynamic_bitset<>::size_type k = 0;

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			vector<bool> bitsVec = table.at((int)y[i][j]);

			for (int l = 0; l < bitsVec.size(); l++)
			{
				bits[k] = bitsVec[l];
				k += 1;
			}

			bitsVec = table.at((int)cb[i][j]);

			for (int l = 0; l < bitsVec.size(); l++)
			{
				bits[k] = bitsVec[l];
				k += 1;
			}

			bitsVec = table.at((int)cr[i][j]);

			for (int l = 0; l < bitsVec.size(); l++)
			{
				bits[k] = bitsVec[l];
				k += 1;
			}
		}
	}

	return bits;
}

void Image::printBits() {
	cout << "Bits of image with height: " << height << "px  width: " << width << "px" << endl;

	for (boost::dynamic_bitset<>::size_type i = 0; i < this->bits.size(); ++i) {
		cout << this->bits[i];
	}

	cout << "\n\n";
}

void Image::printTestResults() {
	cout << endl;
	cout << "Test for image with height: " << height << "px  width: " << width << "px" << endl;
	cout << "Memory for storing with doubles: " << sizeof(double) * 3 * height * width << endl;
	cout << "Memory for storing with bits: " << bits.size() * 1 << endl;
	cout << endl;
}

double Image::getYPixel(int x, int y) {
	if (x >= 0 && x < this->width && y >= 0 && y < this->height) {
		return this->y[y][x];
	}
	else {
		cout << "Invalid coordinates!" << endl;
	}
}

double Image::getCrPixel(int x, int y) {
	if (x >= 0 && x < this->width && y >= 0 && y < this->height) {
		return this->cr[y][x];
	}
	else {
		cout << "Invalid coordinates!" << endl;
		return -130;
	}
}

double Image::getCbPixel(int x, int y) {
	if (x >= 0 && x < this->width && y >= 0 && y < this->height) {
		return this->cb[y][x];
	}
	else {
		cout << "Invalid coordinates!" << endl;
		return -130;
	}
}

void Image::printTestOfFirstTriplet() {
	if (height == 0 || width == 0) {
		return;
	}

	vector<bool> y = this->bitsTable.at(this->y[0][0]);
	vector<bool> cr = this->bitsTable.at(this->cr[0][0]);
	vector<bool> cb = this->bitsTable.at(this->cb[0][0]);

	cout << "Bits for y[0][0] are: ";
	for (int i = 0; i < y.size(); i++) {
		cout << y[i];
	}

	cout << "\nBits saved after compression are: ";
	for (int i = 0; i < y.size(); i++) {
		cout << this->bits[i];
	}

	cout << "\nBits for cb[0][0] are: ";
	for (int i = 0; i < cb.size(); i++) {
		cout << cb[i];
	}

	cout << "\nBits saved after compression are: ";
	for (int i = y.size(); i < y.size() + cb.size(); i++) {
		cout << this->bits[i];
	}

	cout << "\nBits for cr[0][0] are: ";
	for (int i = 0; i < cr.size(); i++) {
		cout << cr[i];
	}

	cout << "\nBits saved after compression are: ";
	for (int i = y.size() + cb.size(); i < y.size() + cb.size() + cr.size(); i++) {
		cout << this->bits[i];
	}

	cout << endl;
}