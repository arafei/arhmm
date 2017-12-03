#ifndef __MATRIX_6_H 
#define __MATRIX_6_H
#include <vector>

template <class T>
class Matrix6{
	public:
	std::vector< std::vector<T> > data;
	Matrix6(int nrow, int ncol, T val = 0){
		data.resize(nrow); // make n rows
		for(int i=0; i < nrow; ++i){
			data[i].resize(ncol,val);
		}
	}
	int rowNums() { return (int)data.size(); }
	int colNums() { return ( data.size() == 0 ) ? 0 : (int)data[0].size(); }
};
#endif