#pragma once
#include "Neuron.hpp"
#include <cmath>
#include "Functions.hpp"

typedef double (*PtrToDoubleFunction)(double);

class Matrix
{
public:
	enum class ACTIVATION_FUNCTION 
	{
		ReLU, ReLU_Leaky, ReLU_Parametric,
		GELU, SiLU, Softplus, ELU, Mish,
		Metallic_mean_function, Sigmoid, Binary_Step
	};
private:
	Network* nt_;
	Weights v_weights_;
	std::vector<Weights> v_all_weights_;
	std::vector<Weights> v_all_x_values_; // [0] = input layer. [.size() - 1] = output layer
	std::vector<Layer>* ptr_neural_network_;
	
	// === F L A G S === 
	bool is_input_x_set_in_vector;
	bool is_entire_x_set_in_vector;
	bool is_weight_vector_made;
	bool is_map_activation_fx_made;

	// === private functions ===
	inline void make_input_x_vector()
	{
		Weights twts;
		std::vector<double> temppp;
		for (auto& i : (*ptr_neural_network_)[0])
		{
			temppp.push_back(i.get_x());
		};
		twts.push_back(temppp);
		v_all_x_values_.push_back(twts);
		is_input_x_set_in_vector = true;
	};
	inline void make_remaining_x_vectors_skeleton()
	{
		Weights twts;
		std::vector<double> temppp;
		for (uint32_t i = 1; i < (*ptr_neural_network_).size() - 1; i++)
		{
			for (auto& j : (*ptr_neural_network_)[i])
			{
				temppp.push_back(0);
			};
			twts.push_back(temppp);
			v_all_x_values_.push_back(twts);
		};
		is_entire_x_set_in_vector = true;
	}
	inline void make_input_and_hidden_weight_vector(const uint32_t tar_layer_index)
	{
		Weights tempwt;
		for (uint32_t i = 0; i < (*ptr_neural_network_)[tar_layer_index - 1].size() - 1; i++)
		{
			std::vector<double> tempdub_I;
			for (uint32_t j = 0; j < (*ptr_neural_network_)[tar_layer_index].size() - 1; j++)
			{
				double rr = RNG<double>::get_random_number(nt_->get_weight_min(), nt_->get_weight_max());
				(*ptr_neural_network_)[tar_layer_index - 1][i].add_fwd_weights
				(&(*ptr_neural_network_)[tar_layer_index][j], rr);
				tempdub_I.push_back(rr);
			};
			tempwt.push_back(tempdub_I);
		};
		v_all_weights_.push_back(tempwt);
	}
	inline void make_output_weight_vector()
	{
		Weights tempwt;
		uint32_t layers = nt_->get_net_properties().num_layers;
		for (uint32_t i = 0; i < (*ptr_neural_network_)[layers - 2].size() - 1; i++)
		{
			std::vector<double> tempdub_I;
			for (uint32_t j = 0; j < (*ptr_neural_network_)[layers - 1].size(); j++)
			{
				double rr = RNG<double>::get_random_number(nt_->get_weight_min(), nt_->get_weight_max());
				(*ptr_neural_network_)[layers - 2][i].add_fwd_weights(&(*ptr_neural_network_)[layers - 1][j],
					rr);
				tempdub_I.push_back(rr);
			};
			tempwt.push_back(tempdub_I);
		};
		v_all_weights_.push_back(tempwt);
	};

	// === private calculations ===
	struct activation_functions
	{
		static inline double ReLU(double x)
		{
			double y;
			if (x > 0) { y = x; }
			else y = 0;
			return y;
		};
		static inline double Sigmoid(double x)
		{
			double y = 1 / (1 + pow(2.718, -x));
			return y;
		};
		static inline void Testy() { std::cout << "Testy" << std::endl; };
	};
	std::unordered_map <ACTIVATION_FUNCTION, PtrToDoubleFunction> m_activation_functions;
	inline std::vector<std::vector<double>> make_bias_matrix(uint32_t layer_pending_a_calculation)
	{
		if (is_entire_x_set_in_vector == false) {
			std::cerr << "ERROR: must set x value vector before making bias vector" << std::endl; exit(65);
		};
		std::vector<double> bias_v;
		double bias;
		bias = v_all_x_values_[layer_pending_a_calculation - 1][0]
							  [v_all_x_values_[layer_pending_a_calculation - 1][0].size() - 1];
		for (uint32_t i = 0; i < v_all_x_values_[layer_pending_a_calculation][0].size() - 1; i++)
		{
			bias_v.push_back(bias);
		};
		std::vector<std::vector<double>> shell; shell.push_back(bias_v);
		return shell;
	}
	inline std::vector<std::vector<double>> calc_dot_product(const std::vector<std::vector<double>>& lhs,
		const std::vector<std::vector<double>>& rhs)
	{
		//The number of columns of the 1st matrix must equal the number of rows of the 2nd matrix.
		//And the result will have the same number of rows as the 1st matrix, and the same number 
		//of columns as the 2nd matrix.
		if (lhs.size() != rhs[0].size()) {
			std::cerr << "Unable to calculate dot product due to nr of rows / columns"
				<< std::endl; exit(44);
		};
		std::vector<std::vector<double>> dot_product;

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
	inline std::vector<std::vector<double>> calc_matrix_addition(const std::vector<std::vector<double>>& lhs,
		const std::vector<std::vector<double>>& rhs)
	{
		if (lhs.size() != rhs.size() || lhs[0].size() != rhs[0].size())
		{
			std::cerr << "ERROR: Matrices must have the same dimensions" << std::endl; exit(67);
		};
		std::vector<std::vector<double>> summed;
		for (uint32_t i = 0; i < lhs.size(); i++)
		{
			for (uint32_t j = 0; j < lhs[0].size(); j++)
			{
				summed[i][j] = lhs[i][j] + rhs[i][j];
			};
		};
		return summed;
	};
	inline std::vector<double> calc_aNeuron_values_by_layer(uint32_t layer, ACTIVATION_FUNCTION af)
	{
		double temp;
		std::vector<double> bias_vector;
		std::vector<std::vector<double>> dot_result;
		dot_result = calc_dot_product(v_all_weights_[layer - 1], v_all_x_values_[layer - 1]);
		std::vector<std::vector<double>> after_bias;
		after_bias = calc_matrix_addition(dot_result, make_bias_matrix(layer));
		std::vector<double> a_values;
		for (auto& i : after_bias[0])
		{
			a_values.push_back(m_activation_functions[af](i));
		};
		return a_values;
	};
	inline double calc_MSE_single(const double output_val, const double label_val)
	{
		// E_p = 1/2 (SUM(d_k - O_k)^2)
		return ((output_val - label_val) * (output_val - label_val));
	};
	inline double calc_MSE_output(const std::vector<double>& label_vals)
	{
		double MSE{};
		uint32_t num_vals = v_all_x_values_[v_all_x_values_.size() - 1][0].size();
		for (uint32_t i = 0; i < num_vals; i++)
		{
			MSE += ((v_all_x_values_[v_all_x_values_.size() - 1][0][i] - label_vals[i]) *
				(v_all_x_values_[v_all_x_values_.size() - 1][0][i] - label_vals[i]));
		};
		MSE = MSE / 2;
		return MSE;
	}
	inline PtrToDoubleFunction calc_derivative(const double(*fxptr)(double))
	{
		
	}

public:
	inline Matrix(Network* nt) : nt_(nt) {
		ptr_neural_network_ = nt_->get_ptr_neural_network();
		is_input_x_set_in_vector = false;
		is_weight_vector_made = false;
		is_map_activation_fx_made = false;
		is_entire_x_set_in_vector = false;
	};
	inline ~Matrix() { nt_ = nullptr; ptr_neural_network_ = nullptr; };
	inline Weights const get_weights_vector() { return v_weights_; };
	inline Weights* get_ptr_to_weights_vector() { return &v_weights_; };
	inline void make_weight_vector()
	{
		for (uint32_t i = 0; i < nt_->get_net_properties().num_layers - 2; i++)
		{
			make_input_and_hidden_weight_vector(i);
		};
		make_output_weight_vector();
		is_weight_vector_made = true;
	};
	inline void make_activation_functions(double x)
	{
		double (*ptr_relu)(double); ptr_relu = &activation_functions::ReLU;
		double (*ptr_sigmoid)(double); ptr_sigmoid = &activation_functions::Sigmoid;
		m_activation_functions.insert({ACTIVATION_FUNCTION::ReLU, ptr_relu});
		m_activation_functions.insert({ ACTIVATION_FUNCTION::Sigmoid, ptr_sigmoid });
		is_map_activation_fx_made = true;
	}
	inline void make_x_values() { make_input_x_vector(); make_remaining_x_vectors_skeleton(); };
	
	inline void calc_neuron_values(ACTIVATION_FUNCTION af_all) 
	{
		if (is_map_activation_fx_made == false) {
			std::cerr << "ERROR: must set activation function map"
				<< std::endl; exit(69);
		};
		for (uint32_t i = 1; i < v_all_x_values_.size(); i++)
		{
			v_all_x_values_[i][0] = calc_aNeuron_values_by_layer(i, af_all);
		};
	};
	inline void calc_neuron_values(ACTIVATION_FUNCTION af_hidden, ACTIVATION_FUNCTION af_output)
	{
		if (is_map_activation_fx_made == false) {
			std::cerr << "ERROR: must set activation function map"
				<< std::endl; exit(69);
		};
		for (uint32_t i = 1; i < v_all_x_values_.size() - 1; i++)
		{
			v_all_x_values_[i][0] = calc_aNeuron_values_by_layer(i, af_hidden);
		};
		v_all_x_values_[v_all_x_values_.size()-1][0] = 
			calc_aNeuron_values_by_layer((v_all_x_values_.size() - 1), af_output);
	};
};