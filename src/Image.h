#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <boost\dynamic_bitset.hpp>

class Image {
private:
	vector<vector<double>> y;
	vector<vector<double>> cr;
	vector<vector<double>> cb;
	vector<pair<int, int>> freq;
	unordered_map<int, vector<bool>> bitsTable;

	void colorTransformAndInit(vector<vector<double>>& red, vector<vector<double>>& green,
		vector<vector<double>>& blue);
	void cosTransform(int x, int y, vector<vector<double>> yNew, vector<vector<double>> crNew,
		vector<vector<double>> cbNew);
	vector<double> g(int x, int y, int u, int v, vector<vector<double>> yNew, vector<vector<double>> crNew,
		vector<vector<double>> cbNew);
	void quantization(int x, int y);
	void initFrequencies(const vector<vector<double>>& px, unordered_map<int, int>& freq);
	void encodeHuffman(int x, int y);
	boost::dynamic_bitset<> getBitRepresentation(unordered_map<int, vector<bool>>& table);

public:
	int height;
	int width;
	boost::dynamic_bitset<> bits;

	Image(const string filename);

	void compress();
	void printBits();

	void printTestResults();
	void printTestOfFirstTriplet();

	double getYPixel(int x, int y);
	double getCrPixel(int x, int y);
	double getCbPixel(int x, int y);
};