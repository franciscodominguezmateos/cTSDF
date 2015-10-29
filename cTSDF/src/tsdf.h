/*
 * tsdf.h
 *
 *  Created on: Oct 28, 2015
 *      Author: Francisco Dominguez
 */

#ifndef TSDF_H_
#define TSDF_H_
#include <iostream>
#include <cmath>

using  namespace std;

template <class T>
class Tsdf {
	int size;
	T *tsdf;
	float min,max;
public:
	Tsdf():size(256),tsdf(new T[size*size*size]),min(0),max(1){};
	inline int getIdx(float f){
		int i;
		float d=max-min;
		float j;
		j=(f-min)/d;
		i=fmax(fmin(j,1.0),0.0)*(size-1);
		return i;
	}
	inline float i2f(int i){
		float j=(float)i/(size-1);
		float d=max-min;
		return min+d*j;
	}
	inline float voxelSize(){
		float d=max-min;
		return d/size;
	}
	inline void setMinMax(float min,float max){
		this->min=min;
		this->max=max;
	}
	inline int getSize(){return size;}
	inline T getVoxel(int i,int j,int k){return tsdf[i+j*size+k*size*size];}
	inline T getVoxel(float x,float y,float z){
		int i,j,k;
		i=getIdx(x);
		j=getIdx(y);
		k=getIdx(z);
		return getVoxel(i,j,k);
	}
	inline void setVoxel(int i,int j,int k,T v){
		//if out of range don't set
		if (i<0 || i>=size || j<0 || j>=size || k<0 || k>=size) return;
		//cout << "i"<<i<<j<<k<<endl;
		tsdf[i+j*size+k*size*size]=v;
	}
	inline void setVoxel(float x,float y,float z,T v){
		int i,j,k;
		i=getIdx(x);
		j=getIdx(y);
		k=getIdx(z);
		//cout << "x"<<x<<y<<z<<endl;
		setVoxel(i,j,k,v);
	}
	inline void clear(T zeros){
	    for(int i=0;i<size*size*size;i++)
	    	tsdf[i]=zeros;
	}
	~Tsdf(){delete tsdf;}
};

#endif /* TSDF_H_ */
