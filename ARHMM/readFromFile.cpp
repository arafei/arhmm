#ifndef __MATRIX_6_H
#define __MATRIX_6_H
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
template <class T>//for generic programming, independent of any particular type.
class Matrix6 {
	public:
	std::vector< std::vector<T> > data;
	Matrix6(int nrow, int ncol, T val = 0) {
		data.resize(nrow); // make n rows
		for(int i=0; i < nrow; ++i)
			data[i].resize(ncol,val); // make n cols with default value val
		}
	int rowNums() { return (int)data.size(); }
	int colNums() { return ( data.size() == 0 ) ? 0 : (int)data[0].size(); }
//	void readFromFile(const char* fileName);
};
#endif
/*
void Matrix6::readFromFile(const char* fileName) {
	// open input file
	std::ifstream ifs(fileName);
	if ( ! ifs.is_open() ) {
		std::cerr << "Cannot open file " << fileName << std::endl;
		abort();
	}
	// set up the tokenizer
	std::string line;
	boost::char_separator <char> sep(" \t");
	// typedef is used to replace long type to a short alias
	typedef boost::tokenizer< boost::char_separator <char> > wsTokenizer;
	// clear the data first
	data.clear();
	int nr = 0, nc = 0;
	// read from file to fill the contents
	while( std::getline(ifs, line) ) {
		if ( line[0] == '#' ) continue; // skip meta-lines starting with #
		wsTokenizer t(line,sep);
		data.resize(nr+1);
		for(wsTokenizer::iterator i=t.begin(); i != t.end(); ++i) {
			data[nr].push_back(boost::lexical_cast <T>(i->c_str()));
			if ( nr == 0 ) ++nc; // count # of columns at the first row
		}
		if ( nc != (int)data[nr].size() ) {
			std::cerr << "The input file is not rectangle at line " << nr << std::endl;
			abort();
		}
		++nr;
	}
}


int main(int argc, char** argv){

Matrix6<double> m(1,1,0);
m.readFromFile(argv[1]);

return 0;
}

*/
