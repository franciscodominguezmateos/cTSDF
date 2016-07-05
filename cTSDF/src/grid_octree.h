/*
 * grid_octree.h
 *
 *  Created on: Jul 4, 2016
 *      Author: francisco
 */

#ifndef GRID_OCTREE_H_
#define GRID_OCTREE_H_


#include <vector>
#include "grid_octree_node.h"

using  namespace std;

template <class T>
class GridOctree {
	int sizeX,sizeY,sizeZ;
	vector<int>voxelsIdx;
	float minX,maxX;
	float minY,maxY;
	float minZ,maxZ;
	GridOctreeNode<T> *nodeRoot;
	int level;
public:
	GridOctree(int sizeX=128,int sizeY=128,int sizeZ=128):
		sizeX(sizeZ),sizeY(sizeY),sizeZ(sizeZ),
		voxelsIdx(vector<int>()),
		minX(0),maxX(1),
		minY(0),maxY(1),
		minZ(0),maxZ(1),
		nodeRoot(new GridOctreeNode<T>()),level(10){
	}
	inline int getChildrenPos(int i,int j,int k,int level){
		int ibit=i>>level & 1;
		int jbit=j>>level & 1;
		int kbit=k>>level & 1;
		int ret= kbit<<2 | jbit<<1 | ibit;
		return ret;
	}
	inline GridOctreeNode<T> *createNodes(GridOctreeNode<T> *nodeRoot,int i,int j,int k,int level){
		int nPos;
		GridOctreeNode<T> *node;
		for(int l=level;l>=0;l--){
			nPos=getChildrenPos(i,j,k,l);
			node=new GridOctreeNode<T>();
			nodeRoot->getChildren()[nPos]=node;
			nodeRoot=node;
		}
		return nodeRoot;
	}
	inline void insertNode(int i,int j,int k,T *value){
		int nPos;
		GridOctreeNode<T> *nodeRoot=this->nodeRoot;
		for(int l=level;l>=0;l--){
			nPos=getChildrenPos(i,j,k,l);
			if(nodeRoot->getChildren()[nPos]!=NULL){
				nodeRoot=nodeRoot->getChildren()[nPos];
			}
			else{
				nodeRoot=createNodes(nodeRoot,i,j,k,l);
				break;
			}
		}
		nodeRoot->value=value;
	}
	inline T getVoxel(int i,int j,int k){
		int nPos;
		T empty;
		GridOctreeNode<T> *nodeRoot=this->nodeRoot;
		for(int l=level;l>=0;l--){
			nPos=getChildrenPos(i,j,k,l);
			if(nodeRoot->getChildren()[nPos]!=NULL){
				nodeRoot=nodeRoot->getChildren()[nPos];
			}
			else{
				return empty;
			}
		}
		return *(nodeRoot->value);
	}
	inline void setVoxel(int i,int j,int k,T v){
		T *p=new T(v);//T need a copy constructor
		insertNode(i,j,k,p);
		int idx=getIdx(i,j,k);
		voxelsIdx.push_back(idx);
	}
	inline void setVoxel(float x,float y,float z,T v){
		int i,j,k;
		i=getXIdx(x);
		j=getYIdx(y);
		k=getZIdx(z);
		setVoxel(i,j,k,v);
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
	inline void clear(T zeros){
	}
	~GridOctree(){}
};




#endif /* GRID_OCTREE_H_ */
