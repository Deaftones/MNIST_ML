#pragma once
#include "Parser.hpp"
#include <utility>
#include <xutility>
#include <thread>
#include <unordered_map>
#include "RNG.h"

class Neuron
{
private:
	double x_val_;
	uint32_t layer_;
	uint32_t number_;
	bool is_bias_;
	std::vector<std::pair<Neuron*, double>> v_fwd_weights_;

public:
	Neuron(uint32_t layer, uint32_t num, bool bias) : layer_(layer), number_(num), is_bias_(bias) {};
	~Neuron() {};
	inline void set_weight(double wt) {};
	inline void set_x_random(double min, double max) { x_val_ = RNG<double>::get_random_number(min, max); };
	inline double get_x() { return x_val_; };
	inline void add_fwd_weights(Neuron* nr, double wt)
	{
		v_fwd_weights_.push_back({ nr, wt });
	};
};

typedef std::vector<Neuron> Layer;
typedef std::vector<std::vector<double>> Weights;

class Network
{
public:
	struct NET_PROPERTIES
	{
		uint32_t num_input_layer;
		uint32_t num_output_layer;
		std::vector<uint32_t> num_neurons_hidden_layers;
		uint32_t num_layers;
		inline void make_hl_vector(const uint32_t layers, const uint32_t numPerLayer)
		{
			for (int i = 0; i < layers; i++)
			{
				num_neurons_hidden_layers.push_back(numPerLayer);
			};
		};
		inline void make_hl_vector(const std::vector<uint32_t>& neuronvec)
		{
			num_input_layer = neuronvec[0];
			num_output_layer = neuronvec[neuronvec.size() - 1];
			for (uint16_t i = 1; i < neuronvec.size() - 1; i++)
			{
				num_neurons_hidden_layers.push_back(neuronvec[i]);
			};
		};
	};
private:
	std::vector<Layer> v_neural_net_;
	std::unordered_map<std::string, Layer*> um_net_layers_;
	NET_PROPERTIES np_;
	double min_, max_;

	// === error flags ===
	bool is_neural_net_created;

	// === private functions ===
	inline void set_layer_weights_random(uint32_t target_layer, double min, double max)
	{
		//args Neuron*, double
		for (uint32_t i = 0; i < v_neural_net_[target_layer].size(); i++)
		{
			for (auto& j : v_neural_net_[target_layer - 1])
			{
				j.add_fwd_weights(&v_neural_net_[target_layer][i], RNG<double>::get_random_number(min, max));
			};
		};
	};

public:
	inline Network(uint32_t layers, uint32_t numPerLayer, uint32_t numOutput, double min, double max)
	{
		is_neural_net_created = false;
		np_.num_input_layer = numPerLayer;
		np_.num_output_layer = numOutput;
		np_.num_layers = layers;
		np_.make_hl_vector(layers, numPerLayer);
		min_ = min;
		max_ = max;
	};
	inline Network(uint32_t layers, std::vector<uint32_t>& neuronsByLayer, double min, double max)
	{
		np_.make_hl_vector(neuronsByLayer);
		is_neural_net_created = false;
		min_ = min;
		max_ = max;
	};
	inline ~Network() {};
	inline double const get_weight_min() { return min_; };
	inline double const get_weight_max() { return max_; };
	inline NET_PROPERTIES const get_net_properties() { return np_; };
	inline void make_neural_net()
	{
		// input layer
		Layer inputL;
		for (uint32_t i = 0; i < np_.num_input_layer; i++)
		{
			Neuron nn(0, i, false);
			inputL.push_back(nn);
		};
		Neuron bb(0, inputL.size(), true); inputL.push_back(bb);
		v_neural_net_.push_back(inputL);

		// hidden layers
		for (uint32_t i = 0; i < np_.num_neurons_hidden_layers.size(); i++)
		{
			Layer hiddenL;
			for (uint32_t j = 0; j < np_.num_neurons_hidden_layers[i]; j++)
			{
				Neuron nn(i + 1, j, false);
				hiddenL.push_back(nn);
			};
			Neuron bb(i + 1, np_.num_neurons_hidden_layers[i], true);
			hiddenL.push_back(bb);
			v_neural_net_.push_back(hiddenL);
		}

		//output layer
		Layer outputL;
		for (uint32_t i = 0; i < np_.num_output_layer; i++)
		{
			Neuron nn(np_.num_layers, i, false);
			outputL.push_back(nn);
		};
		v_neural_net_.push_back(outputL);

		is_neural_net_created = true;
	};
	inline void draw_layer_map()
	{
		if (is_neural_net_created == false) {
			std::cerr << "Unable to draw_layer_map() until neural net is made."
				<< std::endl; exit(55);
		};
		um_net_layers_.insert({ "INPUT", &v_neural_net_[0] });
		um_net_layers_.insert({ "OUTPUT", &v_neural_net_[v_neural_net_.size() - 1] });
		for (int i = 1; i < v_neural_net_.size() - 1; i++)
		{
			std::string aaa = "A" + std::to_string(i);
			um_net_layers_.insert({ aaa, &v_neural_net_[i] });
		};
	};
	inline void set_input_x_vals_random(double min, double max)
	{
		for (auto& i : v_neural_net_[0])
		{
			i.set_x_random(min, max);
		};
	};
	inline void set_input_x_vals_byVector(const std::vector<double>& input_xs)
	{
		if (input_xs.size() != v_neural_net_[0].size()) {
			std::cerr <<
				"ERROR: set_input_x_vals_byVector() has been given a vector argument with more parameters than input layer"
				<< std::endl; exit(44);
		};
		for (uint32_t i = 0; i < v_neural_net_[0].size(); i++)
		{
			v_neural_net_[0][i].set_weight(input_xs[i]);
		};
	};
	inline void set_all_layer_weights_random(double min, double max)
	{
		for (uint32_t i = 1; i < v_neural_net_.size(); i++)
		{
			set_layer_weights_random(i, min, max);
		};
	};
	inline std::vector<Layer>* get_ptr_neural_network() { return &v_neural_net_; };
};



//class Algos
//{
//public:
//	enum class NEURON_WEIGHTING { RANDOM, ZEROS, ONES };
//	struct NEURON
//	{
//		std::vector<NEURON*> v_ptrs_prev_layer;
//		std::vector<NEURON*> v_ptrs_next_layer;
//		std::vector<std::pair<NEURON*, double>> v_input_consts;
//		double weight_;
//		double min_;
//		double max_;
//		int layer_;
//		int number_;
//		bool is_bias_neuron;
//
//
//		NEURON(int layer, int number, double min, double max, bool is_bias) : layer_(layer), min_(min), max_(max)
//		{
//			layer_ = layer; number_ = number; min_ = min, max_ = max; is_bias_neuron = is_bias;
//			weight_ = RNG::get_rnd_double(min, max);
//		};
//		inline void add_to_prev_layer(NEURON& nr) { v_ptrs_prev_layer.push_back(&nr); };
//		inline void add_to_next_layer(NEURON& nr) { v_ptrs_next_layer.push_back(&nr); };
//		/*inline void remove_from_prev_layer(NEURON& nr)
//		{
//			try {
//				if (std::remove(v_ptrs_prev_layer.begin(), v_ptrs_prev_layer.end(), nr) != v_ptrs_prev_layer.end())
//				{
//					std::cout << "Successfully removed neuron from previous layer." << std::endl;
//				}
//				else {
//					throw(nr);
//				}
//
//			}
//			catch (NEURON& neunr)
//			{
//				std::cout << "Neuron in layer " << neunr.layer_ << " with weight " << neunr.weight_ << " is not found in this neuron's previous layer vector." << std::endl;
//			};
//		};
//		inline void remove_from_next_layer(NEURON& nr)
//		{
//			try {
//				if (std::remove(v_ptrs_next_layer.begin(), v_ptrs_next_layer.end(), nr) != v_ptrs_next_layer.end())
//				{
//					std::cout << "Successfully removed neuron from next layer." << std::endl;
//				}
//				else {
//					throw(nr);
//				}
//
//			}
//			catch (NEURON& neunr)
//			{
//				std::cout << "Neuron in layer " << neunr.layer_ << " with weight " << neunr.weight_ << " is not found in this neuron's next layer vector." << std::endl;
//			};
//		};
//		*/
//		inline void add_input_constants()
//		{
//			for (int i = 0; i < v_ptrs_prev_layer.size(); i++)
//			{
//				v_input_consts.push_back({ v_ptrs_prev_layer[i], RNG::get_rnd_double(0, 1) });
//			};
//		};
//		inline void set_weight(double wt) { weight_ = wt; };
//	};
//	
//private:	
//	Parser* parser_;
//	Imager* imager_;
//	std::vector<char>* v_training_labels_;
//	std::vector<IMAGE>* v_training_images_;
//	std::vector<char>* v_testing_labels_;
//	std::vector<IMAGE>* v_testing_images_;
//
//	std::vector<std::vector<NEURON>> v_neural_network_;
//
//	bool are_neurons_connected;
//	bool are_neuron_weightings_added;
//	static uint32_t num_fwd_props;
//	static uint32_t num_bwd_props;
//
//public:
//	Algos() 
//	{
//		parser_ = nullptr; imager_ = nullptr; v_training_images_ = nullptr; v_training_labels_ = nullptr;
//		v_testing_images_ = nullptr; v_testing_labels_ = nullptr;
//		are_neurons_connected = false; are_neuron_weightings_added = false;
//		num_fwd_props = 0; num_bwd_props = 0;
//	};
//	~Algos() {};
//
//	inline void set_training_data(std::vector<char>* labels, std::vector<IMAGE>* imgs)
//	{
//		v_training_labels_ = labels; v_training_images_ = imgs;
//	};
//	inline void set_testing_data(std::vector<char>* labels, std::vector<IMAGE>* imgs)
//	{
//		v_testing_labels_ = labels; v_testing_images_ = imgs;
//	};
//	inline void make_neural_network(uint16_t layers, uint16_t neurons_per_layer, uint16_t possible_outcomes, 
//		double range_1, double range_2) 
//	{
//		for (uint16_t i = 0; i < layers; i++)
//		{
//			std::vector<NEURON> temp_vec;
//			for (uint16_t j = 0; j <= neurons_per_layer; j++)
//			{
//				if (j < neurons_per_layer)
//				{
//					NEURON tt(i, j, range_1, range_2, false);
//					temp_vec.push_back(tt);
//				}
//				else if (j == neurons_per_layer) // making the bias neuron
//				{
//					NEURON tt(i, j, range_1, range_2, true);
//					temp_vec.push_back(tt);
//				};
//			};
//			v_neural_network_.push_back(temp_vec);
//		};
//		std::vector<NEURON> outcomes_vec;
//		for (uint16_t i = 0; i < possible_outcomes; i++)
//		{
//			NEURON tt(layers, i, range_1, range_2, false);
//		};
//		v_neural_network_.push_back(outcomes_vec);
//	};
//	inline void make_neural_network(uint16_t layers, std::vector<uint16_t> neurons_per_layer, NEURON_WEIGHTING nw, int range_1, int range_2)
//	{
//		for (uint16_t i = 0; i < layers; i++) {
//			std::vector<NEURON> temp_vec;
//			for (uint16_t j = 0; j <= neurons_per_layer[i]; j++)
//			{
//				if (j < neurons_per_layer[i])
//				{
//					NEURON tt(i, j, range_1, range_2, false);
//					temp_vec.push_back(tt);
//				}
//				else if (j == neurons_per_layer[i])
//				{
//					NEURON tt(i, j, range_1, range_2, true);
//					temp_vec.push_back(tt);
//				};
//			};
//			v_neural_network_.push_back(temp_vec);
//		};
//	};
//	inline void connect_neurons()
//	{
//		for (uint16_t q = 0; q < v_neural_network_.size(); q++)
//		{
//			for (auto& j : v_neural_network_[q])
//			{
//				if (q == 0)
//				{
//					for (uint16_t h = 0; h < v_neural_network_[1].size(); h++)
//					{
//						j.add_to_next_layer(v_neural_network_[1][h]);
//					};
//
//					for (uint16_t h = 0; h < v_neural_network_[v_neural_network_.size() - 1].size(); h++)
//					{
//						j.add_to_prev_layer(v_neural_network_[v_neural_network_.size() - 1][h]);
//					};
//				}
//				else if (q == (v_neural_network_.size() - 1))
//				{
//					for (uint16_t h = 0; h < v_neural_network_[v_neural_network_.size() - 1].size(); h++)
//					{
//						j.add_to_prev_layer(v_neural_network_[v_neural_network_.size() - 2][h]);
//					};
//
//					for (uint16_t h = 0; h < v_neural_network_[0].size(); h++)
//					{
//						j.add_to_next_layer(v_neural_network_[0][h]);
//					};
//				}
//				else
//				{
//					//add next
//					for (uint16_t n = 0; n < v_neural_network_[q + 1].size(); n++)
//					{
//						j.add_to_next_layer(v_neural_network_[q + 1][n]);
//					};
//					//add previous
//					for (uint16_t p = 0; p < v_neural_network_[q - 1].size(); p++)
//					{
//						j.add_to_prev_layer(v_neural_network_[q - 1][p]);
//					};
//				};
//			};
//		};
//		are_neurons_connected = true;
//	};
//
//	inline void add_weightings()
//	{
//		if (are_neurons_connected == false) {
//			std::cout <<
//				"Unable to add weightings before connect_neurons() is successfully called" << std::endl;
//		}
//		else
//		{
//			for (auto& x : v_neural_network_)
//			{
//				for (auto& y : x)
//				{
//					y.add_input_constants();
//
//					std::cout << "Layer: " << y.layer_ << "\tNumber: " << y.number_ << "\tBias Neuron: " << y.is_bias_neuron
//						<< "\tMy Wt. " << y.weight_
//						<< "\tPrev Wts: ";
//					for (auto& i : y.v_input_consts)
//					{
//						std::cout << i.second << '\t';
//					};
//					std::cout << std::endl;
//				};
//			};
//			are_neuron_weightings_added = true;
//		};
//	};
//
//	inline void forward_propagation()
//	{
//		if (are_neuron_weightings_added == false) {
//			std::cerr << "Unable to feed forward until weightings are added."
//				<< std::endl; exit(88);
//		};
//		for (uint16_t i = 1; i < v_neural_network_.size(); i++)
//		{
//			for (uint16_t n = 0; n < v_neural_network_[i].size() - 1; n++)
//			{
//				double val{};
//				for (uint16_t s = 0; s < v_neural_network_[i - 1].size(); s++)
//				{
//					val = val + (v_neural_network_[i - 1][s].weight_ * v_neural_network_[i][n].v_input_consts[s].second);
//				};
//				v_neural_network_[i][n].set_weight(val);
//			};
//		};
//		num_fwd_props += 1;
//	};
//
//	inline void backward_propagation()
//	{
//
//
//		num_bwd_props += 1;
//	};
//};