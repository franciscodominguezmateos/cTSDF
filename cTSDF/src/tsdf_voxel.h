/*
 * tsdf_voxel.h
 *
 *  Created on: Jul 5, 2016
 *      Author: francisco
 */

#ifndef TSDF_VOXEL_H_
#define TSDF_VOXEL_H_

class TsdfVoxel {
public:
	float x,y,z,d;
	float wd;
	float r,g,b;
	float wc;
//public:
	TsdfVoxel(){
		d=1e32;
	}
	~TsdfVoxel(){
	}
	inline float getD(){
		return d;
	}
	inline float getR(){
		return r;
	}
	inline float getG(){
		return g;
	}
	inline float getB(){
		return b;
	}
	inline void setXYZ(float x,float y,float z){
		this->x=x;this->y=y;this->z=z;
	}
};



#endif /* TSDF_VOXEL_H_ */
