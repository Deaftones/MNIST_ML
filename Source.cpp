#pragma once
//#include "Neuron.hpp"
#include <iostream>
#include <cstdint>
#include <fstream>
#include <ios>
#include <vector>
#include <string>
#include <unordered_map>

const std::string te_img = "C:\\Users\\serba\\source\\repos\\Deaftones\\MNIST_ML\\Data\\t10k_images.idx3-ubyte";
const std::string te_lab = "C:\\Users\\serba\\source\\repos\\Deaftones\\MNIST_ML\\Data\\t10k-labels.idx1-ubyte";
const std::string tr_img = "C:\\Users\\serba\\source\\repos\\Deaftones\\MNIST_ML\\Data\\train-images.idx3-ubyte";
const std::string tr_lab = "C:\\Users\\serba\\source\\repos\\Deaftones\\MNIST_ML\\Data\\train-labels-idx1-ubyte";

const uint16_t label_skip = 8;
const uint16_t img_skip = 16;

typedef double(*DoublePtr)(double, double);

double add(double x, double y) { return x + y; };
double subtract(double x, double y) { return x - y; };
double multiply(double x, double y) { return x * y; };
double divide(double x, double y) { return x / y; };
enum class formulae { add, subtract, multiply, divide };

class wtf
{
public:
	std::unordered_map<formulae, DoublePtr> mymap;
	inline void add_to_map()
	{
		double (*ptr_add)(double, double);
		ptr_add = &add;
		mymap.insert({ formulae::add, ptr_add });
	};
	inline double apply_formula(formulae fr, double vala, double valb)
	{
		double x{};
		x = mymap[fr](vala, valb);
		return x;
	}
};

int main()
{
	wtf lol;
	lol.add_to_map();
	std::cout << lol.apply_formula(formulae::add, 3, 5);
	

	return 0;
}