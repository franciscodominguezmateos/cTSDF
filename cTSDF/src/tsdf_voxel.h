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
	float wr,wg,wb;
//public:
	TsdfVoxel(){
		d=1e32;
	}
	~TsdfVoxel(){
	}
	inline float getD(){
		return d;
	}
};



#endif /* TSDF_VOXEL_H_ */
