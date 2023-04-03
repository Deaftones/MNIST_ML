#pragma once
//#include "Neuron.hpp"
#include <iostream>
#include <cstdint>
#include <fstream>
#include <ios>
#include <vector>
#include <string>
#include <utility>
#include <unordered_map>





class container
{
public:
	struct constant
	{
		char symbol;
		double value;
		constant(char s, double v) : symbol(s), value(v) {};
		~constant() {};
	};
private:
	std::vector<char> v_constants;
	std::vector<constant> v_pair;

public:
	container(std::initializer_list<constant> init) : v_pair(init) {};
	
	
	template <typename ... ARGS>
	void push(ARGS&&...args)
	{
		(void)std::initializer_list<char>{(v_constants.push_back(args), 0)...};
	};
	void print()
	{
		for (auto& i : v_constants) { std::cout << i << std::endl; };
	};
	void print_pair()
	{
		for (auto& i : v_pair)
		{
			std::cout << i.symbol << " = " << i.value << std::endl;
		};
	};
	
};


int main()
{
	container::constant a('a', 2.33);
	container::constant b('b', 35.3);
	container::constant c('c', 94.67);
	container contain{ a, b, c };
	
	contain.print_pair();
	

	return 0;
}