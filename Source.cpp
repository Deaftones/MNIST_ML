#pragma once
//#include "Neuron.hpp"
#include <iostream>
#include <cstdint>
#include <fstream>
#include <ios>
#include <vector>
#include <string>
#include <unordered_map>

/*const std::string te_img = "C:\\Users\\serba\\source\\repos\\Deaftones\\MNIST_ML\\Data\\t10k_images.idx3-ubyte";
const std::string te_lab = "C:\\Users\\serba\\source\\repos\\Deaftones\\MNIST_ML\\Data\\t10k-labels.idx1-ubyte";
const std::string tr_img = "C:\\Users\\serba\\source\\repos\\Deaftones\\MNIST_ML\\Data\\train-images.idx3-ubyte";
const std::string tr_lab = "C:\\Users\\serba\\source\\repos\\Deaftones\\MNIST_ML\\Data\\train-labels-idx1-ubyte";

const uint16_t label_skip = 8;
const uint16_t img_skip = 16;*/

typedef std::vector<std::vector<double>> matrix;

matrix dot_prod_test(const matrix& lhs, const matrix& rhs)
{
	matrix dot_product;
	std::vector<double> temp;

	//building null vector
	for (uint32_t i = 0; i < rhs.size(); i++)
	{
		std::vector<double> zeros;
		for (uint32_t j = 0; j < lhs[i].size(); j++)
		{
			zeros.push_back(0);
		};
		dot_product.push_back(zeros);
	};

	for (uint32_t i = 0; i < rhs.size(); i++)
	{
		for (uint32_t j = 0; j < lhs[i].size(); j++)
		{
			for (uint32_t x = 0; x < lhs.size(); x++)
			{
				dot_product[i][j] += lhs[x][j] * rhs[i][x];
			}
		}
	};

	return dot_product;
};

void print_vector(matrix& x)
{
	for (int i = 0; i < x[0].size(); i++)
	{
		for (int j = 0; j < x.size(); j++)
		{
			std::cout << x[j][i] << "...";
		};
		std::cout << std::endl;
	};
};



int main()
{
	matrix lhs; matrix rhs;
	for (int i = 1; i < 4; i++)
	{
		std::vector<double> temp;
		for (int j = 1; j < 5; j++)
		{
			temp.push_back(j);
		};
		lhs.push_back(temp);
	};

	for (int i = 0; i < 3; i++)
	{
		std::vector<double> temp;
		for (int j = 0; j < 3; j++)
		{
			temp.push_back(j);
		};
		rhs.push_back(temp);
	};

	matrix end = dot_prod_test(lhs, rhs);

	std::cout << "LHS" << std::endl; print_vector(lhs);
	std::cout << "RHS" << std::endl; print_vector(rhs);
	std::cout << "Their Dot product" << std::endl; print_vector(end);

	return 0;
}