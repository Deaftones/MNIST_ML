#pragma once
//#include "Neuron.hpp"
#include <iostream>
#include <cstdint>
#include <fstream>
#include <ios>
#include <vector>
#include <string>
#include <unordered_map>



class Again
{
public:
	template <typename ... ARGS>
	void print(ARGS&& ... args)
	{
		(void)std::initializer_list<int>{ (std::cout << args, 0) ... };
	};
};

template < typename ... ARGS >
void println(ARGS&& ... args)
{
	// the (void) is just used to avoid a 'expression result unused' warning
	(void)std::initializer_list< int >{ (std::cout << args << std::endl, 0) ... };
	std::cout << '\n';
};

class container
{
private:
	std::vector<char> v_constants;

public:
	template <typename ... ARGS>
	void push(ARGS&&...args)
	{
		(void)std::initializer_list<char>{(v_constants.push_back(args), 0)...};
	};
	void print()
	{
		for (auto& i : v_constants) { std::cout << i << std::endl; };
	};
};


int main()
{
	container contain;
	contain.push('A', 'B', 'C', 'D');
	contain.print();

	return 0;
}