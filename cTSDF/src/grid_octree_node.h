/*
 * grid_octree_node.h
 *
 *  Created on: Jul 4, 2016
 *      Author: francisco
 */

#ifndef GRID_OCTREE_NODE_H_
#define GRID_OCTREE_NODE_H_


template <class T>
class GridOctreeNode {
	T *value;
	GridOctreeNode *children[8];
public:
	GridOctreeNode(){
		for(int i=0;i<8;i++)
			children=NULL;
	}
	~GridOctreeNode(){
	}
	GridOctreeNode **getChildren(){
		return children;
	}
};


#endif /* GRID_OCTREE_NODE_H_ */
