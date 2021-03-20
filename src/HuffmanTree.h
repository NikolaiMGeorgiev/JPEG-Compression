#pragma once
#include <vector>
#include <unordered_map>

using namespace std;

struct Node {
	int num;
	int freq;
	Node* left;
	Node* right;

	Node();
	Node(int num, int freq);
};

struct HuffmanTree {
	Node* root = NULL;

	HuffmanTree(vector<pair<int, int>> freq);
	void getBitRepresentations(vector<bool>& bits, unordered_map<int,
		vector<bool>>&table, Node* node);
	void printPostorder(Node* node);
};

bool compareNodes(Node* n1, Node* n2);
