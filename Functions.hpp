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
public:
	struct constant
	{
		char _symbol; double _value;
		constant(char s, double d) : _symbol(s), _value(d) {};
		~constant() {};
	};
	struct variable
	{
		char _symbol; constant _coefficient; double _exponent; char _superscript; char _subscript;
		variable(char s, constant coeff) : _symbol(s), _coefficient(coeff)
		{
			_exponent = 1; _superscript = '%'; _subscript = '%';
		};
		variable(char s, constant coeff, double exp) : _symbol(s), _coefficient(coeff), _exponent(exp)
		{
			_superscript = '%'; _subscript = '%';
		};
		variable(char s, constant coeff, double exp, char sub) : _symbol(s), _coefficient(coeff), _exponent(exp), _subscript(sub)
		{
			_superscript = '%';
		};
		variable(char s, constant coeff, double exp, char sub, char super)
			: _symbol(s), _coefficient(coeff), _exponent(exp), _subscript(sub), _superscript(super) {};
		~variable() {};
	};
private:
	std::vector<constant> _v_constants;
	std::vector<variable> _v_variables;
public:
	Function(std::string function)
	{

	}


};