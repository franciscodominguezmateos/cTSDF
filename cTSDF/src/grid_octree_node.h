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
public:
	//this two attibutes could be a unions since *value is only used in leave and children in *branches
	T *value;
	GridOctreeNode *children[8];// this could be **children to save memory
//public:
	GridOctreeNode(){
		for(int i=0;i<8;i++)
			children[i]=NULL;
	}
	~GridOctreeNode(){
	}
	GridOctreeNode **getChildren(){
		return children;
	}
};


#endif /* GRID_OCTREE_NODE_H_ */
