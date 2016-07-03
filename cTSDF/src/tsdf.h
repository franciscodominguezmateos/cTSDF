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
	vector<int>voxelsIdx;
	float min,max;
public:
	Tsdf(int size=128):size(size),tsdf(new T[size*size*size]),voxelsIdx(vector<int>()),min(0),max(1){};
	inline int getIdx(float f){
		int i;
		float d=max-min;
		float j;
		j=(f-min)/d;
		//i=fmax(fmin(j,1.0),0.0)*(size-1);
		i=j*(size-1);
		return i;
	}
	inline int getIdx(int &i,int &j,int &k){
		return i+j*size+k*size*size;
	}
	inline void getIJKfromIdx(int idx,int &i, int &j,int &k){
		i=idx % size;
		idx/=size;
		j=idx %size;
		idx/=size;
		k=idx;
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
	inline vector<int> &getVoxelsIdx(){
		return voxelsIdx;
	}
	inline void setVoxel(int i,int j,int k,T v){
		//if out of range don't set
		if (i<0 || i>=size || j<0 || j>=size || k<0 || k>=size) return;
		//cout << "i"<<i<<j<<k<<endl;
		int idx=i+j*size+k*size*size;
		tsdf[idx]=v;
		voxelsIdx.push_back(idx);
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
