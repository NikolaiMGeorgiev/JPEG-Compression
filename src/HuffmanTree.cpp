#include <iostream>
#include <vector>
#include <algorithm>
#include <queue>
#include <unordered_map>
#include "HuffmanTree.h"

using namespace std;

Node::Node() {
	left = NULL;
	right = NULL;
}

Node::Node(int num, int freq) {
	this->num = num;
	this->freq = freq;
	left = NULL;
	right = NULL;
}

bool compareNodes(Node* n1, Node* n2) {
		return n1->freq > n2->freq;
}

HuffmanTree::HuffmanTree(vector<pair<int, int>> freq) {
	int size = freq.size() - 1;
	Node* current = new Node;
	Node* previous;
	Node* newNode = new Node;
	vector<Node*> nodes;

	for (int i = 0; i < freq.size(); i++) {
		nodes.push_back(new Node(freq[size - i].first, freq[size - i].second));
	}

	while (!nodes.empty()) {
		previous = nodes.back();
		nodes.pop_back();
		current = nodes.back();
		nodes.pop_back();
		newNode = new Node(-130, current->freq + previous->freq);

		if (previous->freq < current->freq) {
			newNode->left = previous;
			newNode->right = current;
		}
		else {
			newNode->left = current;
			newNode->right = previous;
		}

		nodes.push_back(newNode);

		if (nodes.size() > 1) {
			sort(nodes.begin(), nodes.end(), compareNodes);
		}
		else {
			root = newNode;
			break;
		}
	}
}

void HuffmanTree::getBitRepresentations(vector<bool>& bits, unordered_map<int, 
										vector<bool>>& table, Node* node) {
	if (node == NULL) {
		return;
	}

	if (node->num != -130) {
		if (table.find(node->num) == table.end()) {
			table.insert(make_pair(node->num, bits));
		}
	}

	bits.push_back(0);
	getBitRepresentations(bits, table, node->left);
	bits.pop_back();

	bits.push_back(1);
	getBitRepresentations(bits, table, node->right);
	bits.pop_back();
}

void HuffmanTree::printPostorder(Node* node) {
	if (node == NULL) {
		return;
	}

	printPostorder(node->left);
	printPostorder(node->right);
	cout << node->num << ", " << node->freq << endl;
}