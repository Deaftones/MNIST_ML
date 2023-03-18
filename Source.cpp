#pragma once
#include "Parser.hpp"
#include <iostream>
#include <cstdint>
#include <fstream>
#include <ios>
#include <vector>
#include <string>

const std::string tr_lab = "C:\\Users\\Serban\\source\\repos\\MNIST_ML\\Data\\train-labels-idx1-ubyte";
const std::string tr_img = "C:\\Users\\Serban\\source\\repos\\MNIST_ML\\Data\\train-images.idx3-ubyte";
const std::string te_lab = "C:\\Users\\Serban\\source\\repos\\MNIST_ML\\Data\\t10k-labels.idx1-ubyte";
const std::string te_img = "C:\\Users\\Serban\\source\\repos\\MNIST_ML\\Data\\t10k_images.idx3-ubyte";

const uint16_t label_skip = 8;
const uint16_t img_skip = 16;






int main()
{

	std::vector<char>* ttt;
	Parser pp(Parser::PARSE_MODE::IMAGE, 10000, 784, te_img);
	pp.read_file();
	ttt = pp.get_vec_main();
	std::cout << ttt->size();

	Imager img(&pp);

	img.make_image_vector();
	


	

	return 0;
}