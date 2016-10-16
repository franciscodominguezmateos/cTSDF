/*
 * grid.h
 *
 *  Created on: Jul 3, 2016
 *      Author: Francisco Dominguez
 */

#ifndef GRID_H_
#define GRID_H_
#include <vector>
using  namespace std;

template <class T>
class Grid {
	int sizeX,sizeY,sizeZ;
	T *grid;
	vector<int>voxelsIdx;
	float minX,maxX;
	float minY,maxY;
	float minZ,maxZ;
public:
	Grid(int sizeX=128,int sizeY=128,int sizeZ=128):
		sizeX(sizeZ),sizeY(sizeY),sizeZ(sizeZ),
		grid(new T[sizeX*sizeY*sizeZ]),
		voxelsIdx(vector<int>()),
		minX(0),maxX(1),
		minY(0),maxY(1),
		minZ(0),maxZ(1){}
	inline int getChildrenPos(int i,int j,int k,int level){
		int ibit=i>>level & 1;
		int jbit=j>>level & 1;
		int kbit=k>>level & 1;
		int ret= kbit<<2 | jbit <<1 | ibit;
		return ret;
	}
	inline int getXIdx(float f){
		int i;
		float d=maxX-minX;
		float j;
		j=(f-minX)/d;
		//i=fmax(fmin(j,1.0),0.0)*(size-1);
		i=j*(sizeX-1);
		return i;
	}
	inline int getYIdx(float f){
		int i;
		float d=maxY-minY;
		float j;
		j=(f-minY)/d;
		//i=fmax(fmin(j,1.0),0.0)*(size-1);
		i=j*(sizeY-1);
		return i;
	}
	inline int getZIdx(float f){
		int i;
		float d=maxZ-minZ;
		float j;
		j=(f-minZ)/d;
		//i=fmax(fmin(j,1.0),0.0)*(size-1);
		i=j*(sizeZ-1);
		return i;
	}
	inline int getIdx(int &i,int &j,int &k){
		return i+j*sizeX+k*sizeX*sizeY;
	}
	inline void getIJKfromIdx(int idx,int &i, int &j,int &k){
		i=idx % sizeX;
		idx/=sizeX;
		j=idx %sizeY;
		idx/=sizeY;
		k=idx;
	}
	inline float i2X(int i){
		float j=(float)i/(sizeX-1);
		float d=maxX-minX;
		return minX+d*j;
	}
	inline float j2Y(int i){
		float j=(float)i/(sizeY-1);
		float d=maxY-minY;
		return minY+d*j;
	}
	inline float k2Z(int i){
		float j=(float)i/(sizeZ-1);
		float d=maxZ-minZ;
		return minZ+d*j;
	}
	inline float voxelSizeX(){
		float d=maxX-minX;
		return d/sizeX;
	}
	inline float voxelSizeY(){
		float d=maxY-minY;
		return d/sizeY;
	}
	inline float voxelSizeZ(){
		float d=maxZ-minZ;
		return d/sizeZ;
	}
	inline void setMinMax(float minX,float maxX,float minY,float maxY,float minZ,float maxZ){
		this->minX=minX;
		this->maxX=maxX;
		this->minY=minY;
		this->maxY=maxY;
		this->minZ=minZ;
		this->maxZ=maxZ;
	}
	inline int getSizeX(){return sizeX;}
	inline int getSizeY(){return sizeY;}
	inline int getSizeZ(){return sizeX;}
	inline T getVoxel(int i,int j,int k){return grid[getIdx(i,j,k)];}
	inline T getVoxel(float x,float y,float z){
		int i,j,k;
		i=getXIdx(x);
		j=getYIdx(y);
		k=getYIdx(z);
		return getVoxel(i,j,k);
	}
	inline vector<int> &getVoxelsIdx(){
		return voxelsIdx;
	}
	inline void setVoxel(int i,int j,int k,T v){
		//if out of range don't set
		if (i<0 || i>=sizeX || j<0 || j>=sizeY || k<0 || k>=sizeZ) return;
		//cout << "i"<<i<<j<<k<<endl;
		int idx=getIdx(i,j,k);
		grid[idx]=v;
		voxelsIdx.push_back(idx);
	}
	inline void setVoxel(float x,float y,float z,T v){
		int i,j,k;
		i=getXIdx(x);
		j=getYIdx(y);
		k=getZIdx(z);
		//cout << "x"<<x<<y<<z<<endl;
		setVoxel(i,j,k,v);
	}
	inline void clear(T zeros){
	    for(int i=0;i<sizeX*sizeY*sizeZ;i++)
	    	grid[i]=zeros;
	}
	~Grid(){delete grid;}
};

#endif /* GRID_H_ */
