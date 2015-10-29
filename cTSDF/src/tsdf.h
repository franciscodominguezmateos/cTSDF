/*
 * tsdf.h
 *
 *  Created on: Oct 28, 2015
 *      Author: francisco
 */

#ifndef TSDF_H_
#define TSDF_H_
#include <cmath>

template <class T>
class Tsdf {
	int size;
	T *tsdf;
	float min,max;
public:
	Tsdf():size(256),tsdf(new T[size*size*size]),min(0),max(1){};
	inline int getCord(float f){
		int i;
		float d=max-min;
		i=(f-min)/d;
		i=fmax(fmin(i,1.0),0.0)*size;
		return i;
	}
	inline void setMinMax(float min,float max){
		this->min=min;
		this->max=max;
	}
	inline T getVoxel(int i,int j,int k){return tsdf[i+j*size+k*size*size];}
	inline T getVoxel(float x,float y,float z){
		int i,j,k;
		i=getCord(x);
		j=getCord(y);
		k=getCord(z);
		return getVoxel(i,j,k);
	}
	inline void setVoxel(int i,int j,int k,T v){tsdf[i+j*size+k*size*size]=v;}
	inline void setVoxel(float x,float y,float z,T v){
		int i,j,k;
		i=getCord(x);
		j=getCord(y);
		k=getCord(z);
		setVoxel(i,j,k,v);
	}
	inline void clear(T zeros){
	    for(int i=0;i<size;i++)
	    	for(int j=0;j<size;j++)
	    		for(int k=0;k<size;k++)
	    			setVoxel(i,j,k,zeros);
	}
	~Tsdf(){delete tsdf;}
};

#endif /* TSDF_H_ */
