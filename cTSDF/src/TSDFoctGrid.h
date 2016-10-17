/*
 * TSDFoctGrid.h
 *
 *  Created on: Oct 17, 2016
 *      Author: francisco
 */

#ifndef TSDFOCTGRID_H_
#define TSDFOCTGRID_H_
#include <stdio.h>
#include <vector>
#include <utility>
#include <unistd.h>
#include <GL/glut.h>
#include <depthImage.h>
#include "grid_octree.h"
#include "tsdf_voxel.h"
extern "C" {
 #include "poligonise.h"
}

using namespace std;

typedef float S;

class TSDFoctGrid{
public:
	bool wires;
	bool boxes;
	int level;
	float l;//width/2 of the cube grid
	GridOctree<TsdfVoxel> g;
	vector<Point3f> vpts;
	vector<TRIANGLE> mesh;
	struct colorVertex{
		Point3f c[3];
	};
	vector<colorVertex> colors;
	TSDFoctGrid(int lev=10,float lp=2):level(lev),l(lp),g(GridOctree<TsdfVoxel>(1<<level,1<<level,1<<level)){
		g.setLevel(level);
	    g.setMinMax(-l,l,
	    		    -l,l,
					-l,l);
	    wires=true;
	    boxes=false;
	}
	~TSDFoctGrid();
	void computeNormals(TRIANGLE *t){
		XYZ v1,v2;
		v1.x=t->p[1].x-t->p[0].x;
		v1.y=t->p[1].y-t->p[0].y;
		v1.z=t->p[1].z-t->p[0].z;
		v2.x=t->p[2].x-t->p[1].x;
		v2.y=t->p[2].y-t->p[1].y;
		v2.z=t->p[2].z-t->p[1].z;
		XYZ n;
		n.x=v1.y*v2.z-v1.z*v2.y;
		n.y=v1.x*v2.z-v1.z*v2.x;
		n.z=v1.x*v2.y-v1.y*v2.x;
		double d=sqrt(n.x*n.x+n.y*n.y+n.z*n.z);
		n.x/=d;
		n.y/=d;
		n.z/=d;
		t->n[0].x=n.x;
		t->n[0].y=n.y;
		t->n[0].z=n.z;
		t->n[1].x=n.x;
		t->n[1].y=n.y;
		t->n[1].z=n.z;
		t->n[2].x=n.x;
		t->n[2].y=n.y;
		t->n[2].z=n.z;
	}
	float truncationDistance(float d){
		return ((int)d+1)*1.5;
	}
	inline S getTsdf(S x,S y,S z){
		TsdfVoxel *vxl;
		vxl=g.getVoxelPtr(x,y,z);
		if(vxl!=NULL){
			return vxl->getD();//this should be trilineal interpolated
		}
		else{
			return g.voxelSizeX()*2;//max distance TODO
		}
	}
	inline Point3f getTsdfColor(S x,S y,S z){
		TsdfVoxel *vxl;
		vxl=g.getVoxelPtr(x,y,z);
		if(vxl!=NULL){
			return Point3f(vxl->getR()/255.0,vxl->getG()/255.0,vxl->getB()/255.0);//this should be trilineal interpolated
		}
		else{
			return Point3f(1,0,0);//max distance TODO
		}
	}
	//not trilineal interpolated
	void computeVertexNormals(TRIANGLE *t,colorVertex *cv){
		int i,j,k;
		S x1,y1,z1;
		S d,dx,dy,dz;
		S x,y,z;
		XYZ n;
		for(int idx=0;idx<3;idx++){
			x=t->p[idx].x;
			y=t->p[idx].y;
			z=t->p[idx].z;
			g.XYZ2ijk(x,y,z,i,j,k);
			x1=g.i2X(i+1);
			y1=g.j2Y(j+1);
			z1=g.k2Z(k+1);
			d =getTsdf(x ,y ,z );
			dx=getTsdf(x1,y ,z );
			dy=getTsdf(x ,y1,z );
			dz=getTsdf(x ,y ,z1);
			n.x=dx-d;
			n.y=dy-d;
			n.z=dz-d;
			double m=sqrt(n.x*n.x+n.y*n.y+n.z*n.z);
			n.x/=m;
			n.y/=m;
			n.z/=m;
			t->n[idx].x=n.x;
			t->n[idx].y=n.y;
			t->n[idx].z=n.z;
			//Color calculation
			//Alway must exit voxel i,j,k
			Point3f c=getTsdfColor(x,y,z);
			cv->c[idx]=c;;
		}
	}
	void buildMesh(){
		TRIANGLE triangles[10];
	    GRIDCELL grid;
	    int i,j,k;
	    TsdfVoxel vxl,empty;
	    for(Idx idx:g.getVoxelsIdx()){
	       	g.getIJKfromIdx(idx,i,j,k);
	    			grid.p[0].x = g.i2X(i);
	    			grid.p[0].y = g.j2Y(j);
	    			grid.p[0].z = g.k2Z(k);
	    			vxl=g.getVoxel(i,j,k);
	    			if(vxl.d==empty.d) continue;
	    			grid.val[0] = vxl.d;
	                grid.p[1].x = g.i2X(i+1);
	                grid.p[1].y = g.j2Y(j);
	                grid.p[1].z = g.k2Z(k);
	    			vxl=g.getVoxel(i+1,j,k);
	    			if(vxl.d==empty.d) continue;
	    			grid.val[1] =  vxl.d;
	                grid.p[2].x = g.i2X(i+1);
	                grid.p[2].y = g.j2Y(j+1);
	                grid.p[2].z = g.k2Z(k);
	    			vxl=g.getVoxel(i+1,j+1,k);
	    			if(vxl.d==empty.d) continue;
	    			grid.val[2] =  vxl.d;
	                grid.p[3].x = g.i2X(i);
	                grid.p[3].y = g.j2Y(j+1);
	                grid.p[3].z = g.k2Z(k);
	    			vxl=g.getVoxel(i,j+1,k);
	    			if(vxl.d==empty.d) continue;
	    			grid.val[3] =  vxl.d;
	                grid.p[4].x = g.i2X(i);
	                grid.p[4].y = g.j2Y(j);
	                grid.p[4].z = g.k2Z(k+1);
	    			vxl=g.getVoxel(i,j,k+1);
	    			if(vxl.d==empty.d) continue;
	    			grid.val[4] =  vxl.d;
	                grid.p[5].x = g.i2X(i+1);
	                grid.p[5].y = g.j2Y(j);
	                grid.p[5].z = g.k2Z(k+1);
	    			vxl=g.getVoxel(i+1,j,k+1);
	    			if(vxl.d==empty.d) continue;
	    			grid.val[5] =  vxl.d;
	                grid.p[6].x = g.i2X(i+1);
	                grid.p[6].y = g.j2Y(j+1);
	                grid.p[6].z = g.k2Z(k+1);
	    			vxl=g.getVoxel(i+1,j+1,k+1);
	    			if(vxl.d==empty.d) continue;
	    			grid.val[6] =  vxl.d;
	                grid.p[7].x = g.i2X(i);
	                grid.p[7].y = g.j2Y(j+1);
	                grid.p[7].z = g.k2Z(k+1);
	    			vxl=g.getVoxel(i,j+1,k+1);
	    			if(vxl.d==empty.d) continue;
	    			grid.val[7] =  vxl.d;
	    			int	n = PolygoniseCube(grid,0.0,triangles);
					colorVertex color;
	    			for (int l=0;l<n;l++){
	    				//computeNormals(&triangles[l]);
	    				computeVertexNormals(&triangles[l],&color);
	    				mesh.push_back(triangles[l]);
	    				colors.push_back(color);
	    			}
	    }
	}

	// 30/7/2016 adding tau=truncation Distance
	void updateGrid(DepthImage &di1){
		Point3f p3D, p3Dg;
		TsdfVoxel vxl;
		TsdfVoxel *vxlPtr;
		float pd;
		float sz=g.voxelSizeZ();//grid Z size in m
		for(int u=0;u<di1.cols();u++){
	    	//if(u%300==0) cout << u << endl;
	    	for(int v=0;v<di1.rows();v++){
	    		//units metres
	    		if(di1.isGoodDepthPixel(u,v)){
					Point3f rp3D=di1.getPoint3D(u,v);
					Vec3b c=di1.getColor(u,v);
					float d=di1.getDepth(u,v);
					if(d>1.5) continue;
					if(g.isOut(rp3D.x,rp3D.y,rp3D.z)) continue;
					float tau=sz*3;//truncationDistance(d);
					float tau2=tau;
					float itau=sz/2.0;
					//cout << "d   =" << d << endl;
	//				cout << "sz  =" << sz << endl;
	//				cout << "tau=" << tau << endl;
	//				cout << "itau=" << itau  <<endl;
	    			for(float dt=-tau;dt<=sz;dt+=itau){
	    				//cout << "dt="<< dt <<endl;
	    				p3D=di1.getPoint3Ddeep(u,v,d+dt);
	    				pd=di1.projectiveDistance(p3D);
	    				vxl.x=rp3D.x;
	    				vxl.y=rp3D.y;
	    				vxl.z=rp3D.z;
	    				vxl.r=c[2];
	    				vxl.g=c[1];
	    				vxl.b=c[0];
	    				vxl.d=pd;
	    				vxl.wd=1;///tau2/(vxl.z*vxl.z);
	    				vxl.wr=1;///tau2/(vxl.z*vxl.z);
	    				vxl.wg=1;///tau2/(vxl.z*vxl.z);
	    				vxl.wb=1;///tau2/(vxl.z*vxl.z);
	    				p3Dg=di1.toGlobal(p3D);
	    				vxlPtr=g.getVoxelPtr(p3Dg.x,p3Dg.y,p3Dg.z);
	    				if(vxlPtr==NULL){
	        				g.setVoxel(p3Dg.x,p3Dg.y,p3Dg.z,vxl);
	    				}
	    				else{
	    					float &D=vxlPtr->d;
	    					float &W=vxlPtr->wd;
	    					float &R=vxlPtr->r;
	    					float &G=vxlPtr->g;
	    					float &B=vxlPtr->b;
	    					float &WR=vxlPtr->wr;
	    					float &WG=vxlPtr->wg;
	    					float &WB=vxlPtr->wb;
	    					float &d=vxl.d;
	    					float &w=vxl.wd;
	    					float &r=vxl.r;
	    					float &g=vxl.g;
	    					float &b=vxl.b;
	    					float &wr=vxl.wr;
	    					float &wg=vxl.wg;
	    					float &wb=vxl.wb;
	    					D=(W*D+w*d)/(W+w);
	    					W+=w;
	    					R=(WR*R+wr*r)/(WR+wr);
	    					G=(WR*G+wg*g)/(WG+wg);
	    					B=(WB*B+wb*b)/(WB+wb);
	    					WR+=wr;
	    					WG+=wg;
	    					WB+=wb;
	    				}
	    			}
	    		}
	    	}
	    }
	}
	void rayMarching(GridOctree<TsdfVoxel> &g,DepthImage &di1){
		float sz=g.voxelSizeZ();
		for(int u=0;u<di1.cols();u+=1){
			cout << "u="<<u<<endl;
			for(int v=0;v<di1.rows();v+=1){
				TsdfVoxel *vxlPtr=NULL;
				float d=10e32;
				float dt=0.5;
				float idt=sz*4;
				for(;dt<g.getSizeZ() && d>0.0;dt+=idt){
					Point3f p3D=di1.getPoint3Ddeep(u,v,dt*g.voxelSizeZ());
					Point3f p3Dg=di1.toGlobal(p3D);
					vxlPtr=g.getVoxelPtr(p3Dg.x,p3Dg.y,p3Dg.z);
					if(vxlPtr!=NULL){
						d=vxlPtr->d;
						idt=d;
					}
				}
				if(d<0)
					di1.setDepth(u,v,dt*g.voxelSizeZ());
				else
					di1.setDepth(u,v,0.0);
	 		}
		}
	}
	void glDrawMesh(){
	      glBegin(GL_TRIANGLES);
	      for(int i=0;i<mesh.size();i++){
	    	  TRIANGLE t=mesh[i];
	    	  colorVertex color=colors[i];
	          if(!wires)
	        	  glColor3f(1,1,1);
	          if(wires) glColor3f(color.c[0].x,color.c[0].y,color.c[0].z);
	          glNormal3f(t.n[0].x,t.n[0].y,t.n[0].z);
	    	  glVertex3f(t.p[0].x,t.p[0].y,t.p[0].z);
	    	  if(wires) glColor3f(color.c[1].x,color.c[1].y,color.c[1].z);
	          glNormal3f(t.n[1].x,t.n[1].y,t.n[1].z);
	    	  glVertex3f(t.p[1].x,t.p[1].y,t.p[1].z);
	    	  if(wires)glColor3f(color.c[2].x,color.c[2].y,color.c[2].z);
	          glNormal3f(t.n[2].x,t.n[2].y,t.n[2].z);
	    	  glVertex3f(t.p[2].x,t.p[2].y,t.p[2].z);
	      }
	      //Box Grid
	      glEnd();
	      glPushMatrix();
	      //glTranslatef(-l,-l,-l);
	      glColor3f(1,0,1);
	      glutWireCube(l*2);
	      //Axes
	      glBegin(GL_LINES);
	       glColor3f(1,0,0);
	       glVertex3f(0,0,0);
	       glVertex3f(l,0,0);
	       glColor3f(0,1,0);
	       glVertex3f(0,0,0);
	       glVertex3f(0,l,0);
	       glColor3f(0,0,1);
	       glVertex3f(0,0,0);
	       glVertex3f(0,0,l);
	      glEnd();
	      glPopMatrix();
	      //boxes of voxels
	      if(boxes){
	      int i,j,k;
	      float x,y,z;
	      float s=g.voxelSizeX();
	      for(Idx idx:g.getVoxelsIdx()){
	     	 glPushMatrix();
	    	 glColor3f(0,0,1);
	         g.getIJKfromIdx(idx,i,j,k);
	         g.ijk2XYZ(i,j,k,x,y,z);
	         glTranslatef(x,y,z);
	         glutWireCube(s);
	         glPopMatrix();
	      }
	      }
	}

};

#endif /* TSDFOCTGRID_H_ */
