#pragma once
#include "Neuron.hpp"
#include <cmath>

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
	std::vector<std::vector<double>> v_all_x_values_;
	std::vector<Layer>* ptr_neural_network_;
	
	// === F L A G S === 
	bool is_input_x_set_in_vector;
	bool is_weight_vector_made;
	bool is_map_activation_fx_made;

	// === private functions ===
	inline void make_input_x_vector()
	{
		std::vector<double> temppp;
		for (auto& i : (*ptr_neural_network_)[0])
		{
			temppp.push_back(i.get_x());
		};
		v_all_x_values_.push_back(temppp);
		is_input_x_set_in_vector = true;
	};
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

public:
	inline Matrix(Network* nt) : nt_(nt) {
		ptr_neural_network_ = nt_->get_ptr_neural_network();
		is_input_x_set_in_vector = false;
		is_weight_vector_made = false;
		is_map_activation_fx_made = false;
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
	inline std::vector<std::vector<double>> calc_dot_product(std::vector<std::vector<double>> lhs,
		std::vector<std::vector<double>> rhs)
	{
		//The number of columns of the 1st matrix must equal the number of rows of the 2nd matrix.
		//And the result will have the same number of rows as the 1st matrix, and the same number 
		//of columns as the 2nd matrix.
		std::vector<std::vector<double>> dot;
		std::vector<double> temp;
		if (lhs.size() != rhs[0].size()) {
			std::cerr << "Unable to calculate dot product due to nr of rows / columns"
				<< std::endl; exit(44);
		};
		for (uint32_t i = 0; i < lhs.size(); i++)
		{
			double add{};
			for (uint32_t j = 0; j < lhs[i].size(); j++)
			{
				double d{};
				d = lhs[i][j] + rhs
			}
		}

	};
	inline void calc_neuron_values(ACTIVATION_FUNCTION af_all) 
	{
		if (is_map_activation_fx_made == false) {
			std::cerr << "ERROR: must set activation function map"
				<< std::endl; exit(69);
		};
		for (uint32_t i = 0; i < v_all_weights_.size(); i++)
		{
			
		}
	};
	inline void calc_neuron_values(ACTIVATION_FUNCTION af_hidden, ACTIVATION_FUNCTION af_output);
	

	

};