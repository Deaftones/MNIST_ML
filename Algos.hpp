#pragma once
#include "Parser.hpp"
#include <random>

class Algos
{
public:
	enum class NEURON_WEIGHTING { RANDOM, ZEROS, ONES };
	struct NEURON
	{
		std::vector<NEURON*> v_ptrs_prev_layer;
		std::vector<NEURON*> v_ptrs_next_layer;
		double weight_;
		double min_;
		double max_;
		int layer_;

		NEURON(int layer, double min, double max, NEURON_WEIGHTING nw) : layer_(layer), min_(min), max_(max)
		{
			switch (nw)
			{
			case NEURON_WEIGHTING::RANDOM:
			{

			} break;
			default: { } break;
			}
		};
		inline void add_to_prev_layer(NEURON& nr) { v_ptrs_prev_layer.push_back(&nr); };
		inline void add_to_next_layer(NEURON& nr) { v_ptrs_next_layer.push_back(&nr); };
		inline void remove_from_prev_layer(NEURON& nr)
		{
			try {
				if (std::remove(v_ptrs_prev_layer.begin(), v_ptrs_prev_layer.end(), nr) != v_ptrs_prev_layer.end())
				{
					std::cout << "Successfully removed neuron from previous layer." << std::endl;
				}
				else {
					throw(nr);
				}

			}
			catch (NEURON& neunr)
			{
				std::cout << "Neuron in layer " << neunr.layer_ << " with weight " << neunr.weight_ << " is not found in this neuron's previous layer vector." << std::endl;
			};
		};
		inline void remove_from_next_layer(NEURON& nr)
		{
			try {
				if (std::remove(v_ptrs_next_layer.begin(), v_ptrs_next_layer.end(), nr) != v_ptrs_next_layer.end())
				{
					std::cout << "Successfully removed neuron from next layer." << std::endl;
				}
				else {
					throw(nr);
				}

			}
			catch (NEURON& neunr)
			{
				std::cout << "Neuron in layer " << neunr.layer_ << " with weight " << neunr.weight_ << " is not found in this neuron's next layer vector." << std::endl;
			};
		};
	};
	
private:	
	Parser* parser_;
	Imager* imager_;
	std::vector<char>* v_training_labels_;
	std::vector<IMAGE>* v_training_images_;
	std::vector<char>* v_testing_labels_;
	std::vector<IMAGE>* v_testing_images_;

	std::vector<NEURON> v_neural_network_;

	inline void connect_neurons()
	{
		for (auto& i : v_neural_network_)
		{
			for ()
		}
	}

public:
	Algos(Parser* pp, Imager* ii) : parser_(pp), imager_(ii) {};
	~Algos() {};

	inline void set_training_data(std::vector<char>* labels, std::vector<IMAGE>* imgs)
	{
		v_training_labels_ = labels; v_training_images_ = imgs;
	};
	inline void set_testing_data(std::vector<char>* labels, std::vector<IMAGE>* imgs)
	{
		v_testing_labels_ = labels; v_testing_images_ = imgs;
	};
	inline void make_neural_network(uint16_t layers, uint16_t neurons_per_layer, NEURON_WEIGHTING nw, int range_1, int range_2) 
	{
		for (uint16_t i = 0; i < layers; i++)
		{
			for (uint16_t j = 0; j < neurons_per_layer; j++)
			{
				NEURON tt(i, range_1, range_2, nw);
				v_neural_network_.push_back(tt);
			};
		};
	};
	inline void make_neural_network(uint16_t layers, std::vector<uint16_t> neurons_per_layer, NEURON_WEIGHTING nw, int range_1, int range_2)
	{
		for (uint16_t i = 0; i < layers; i++) {
			for (uint16_t j = 0; j < neurons_per_layer[i]; j++)
			{
				NEURON tt(i, range_1, range_2, nw);
				v_neural_network_.push_back(tt);
			};
		};

	};

};