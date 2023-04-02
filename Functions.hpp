#pragma once
#include <vector>
#include <string>
#include <cmath>
#include <cstdint>
#include <utility>

class Fx
{
public:
	virtual ~Fx() { };
	virtual void derivative() = 0;
	virtual void integral() = 0;
	virtual void partial_derivative() = 0;
};


class Function
{
private:
	std::vector<char> v_constants;
	std::vector<char> v_variables;
public:
	Function() {};
	~Function() {};
	template <typename Arg, Arg...args>
	inline void add_constants(Arg args...)
	{
		v_constants.push_back(&args...);
	}
};