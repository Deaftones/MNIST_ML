#pragma once
#include <iostream>
#include <cstdint>
#include <string>
#include <vector>
#include <fstream>
#include <array>
#include <ios>
#include <stdlib.h>

class Imager;
class Parser;

struct IMAGE
{
	uint64_t id_;
	const uint8_t rows_ = 28;
	const uint8_t cols_ = 28;
	std::vector<char> image_by_rows;
	IMAGE(uint64_t id, std::vector<char>& raw_data)
	{
		id_ = id;
		image_by_rows = raw_data;
		
	};
	~IMAGE() 
	{
		
	};
};

class Parser
{
public:
	enum class PARSE_MODE { LABEL, IMAGE };
private:
	std::vector<char> vec_main;
	uint32_t num_items_, chunk_size_;
	std::string src_;
	PARSE_MODE pm_;
	uint16_t offset_;
	

public:
	Parser(PARSE_MODE pm, uint32_t num, uint32_t chunksize, std::string src) : pm_(pm), num_items_(num), chunk_size_(chunksize), src_(src) 
	{
		if (pm_ == PARSE_MODE::LABEL) { offset_ = 8; }
		else offset_ = 16;
		
	};

	inline void read_file()
	{
		std::ifstream file;
		file.open(src_, std::ios::binary);
		if (!file.is_open()) { std::cerr << "ERROR: Parser::read_file() is unable to read file." << std::endl; exit(44); };
		file.unsetf(std::ios::skipws);
		std::streampos file_size;
		file.seekg(0, std::ios::end);
		file_size = file.tellg();
		file.seekg(offset_, std::ios::beg);
;
		vec_main.reserve(((uint16_t)file_size - offset_));
		vec_main.insert(vec_main.begin(),
			std::istream_iterator<char>(file),
			std::istream_iterator<char>());
	};
	inline std::vector<char>* get_vec_main() { return &vec_main; };
	inline uint32_t get_num_items() { return num_items_; };
	inline uint32_t get_chunk_size() { return chunk_size_; };
};

class Imager
{
private:
	Parser* parser_;
	std::vector<IMAGE> v_images_;
	uint64_t counter_;

public:
	Imager(Parser* p) : parser_(p) { counter_ = 0; };

	inline void make_image_vector()
	{
		uint32_t nItems = parser_->get_num_items();
		uint32_t sChunk = parser_->get_chunk_size();
		std::vector<char>* v_ptr = parser_->get_vec_main();
		for (uint32_t i = 0; i < nItems; i++)
		{
			std::vector<char> v_temp;
			for (uint32_t j = 0; j < sChunk; j++)
			{
				v_temp.push_back(v_ptr->at(j + (i * sChunk)));
			};
			IMAGE iii(counter_, v_temp);
			v_images_.push_back(iii);
			counter_++;
			std::cout << "Created image number " << counter_ << " of " << nItems << std::endl;
		};
	};

	~Imager(){}
};

