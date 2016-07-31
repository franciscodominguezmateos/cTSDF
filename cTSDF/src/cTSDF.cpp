//============================================================================
// Name        : cTSDF.cpp
// Author      : Francico Dominguez
// Version     :
// Copyright   : Your copyright notice
// Description : easy TSDF in C++, Ansi-style
//============================================================================
#include <stdio.h>
#include <vector>
#include <utility>
#include <unistd.h>
#include <GL/glut.h>
#include <depthImage.h>
//#include "tsdf.h"
//#include "grid.h"
#include "grid_octree.h"
#include "tsdf_voxel.h"
extern "C" {
 #include "poligonise.h"
}

using namespace std;

GLint ancho=400;
GLint alto=400;
int hazPerspectiva = 1;
GLfloat angle = 0.0;
GLfloat yaw = 0.0;
GLfloat roll = 0.0;
GLfloat pitch = 0.0;
GLfloat t=-3.0f;

DepthImage di1,di2;
int level=9;
GridOctree<TsdfVoxel> g(1<<level,1<<level,1<<level);
vector<Point3f> vpts;
vector<TRIANGLE> mesh;
vector<Point3f> colors;

bool wires=true;
bool friccion=true;

GLfloat light0_ambient[] ={0.2, 0.2, 0.2, 1.0};
GLfloat light0_diffuse[] ={0.8, 0.8, 0.8, 1.0};

void displayMe(void)
{
    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
//    gluLookAt (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    GLfloat lightpos[] = {0.0, 15., 5., 0.0};
    glLightfv(GL_LIGHT0, GL_POSITION, lightpos);
    glTranslatef(0.0f, 0.0f, t);
    glRotatef(yaw  ,0.0,1.0,0.0);
    glRotatef(pitch,1.0,0.0,0.0);
    glRotatef(roll ,0.0,0.0,1.0);
    glPushMatrix();
     //glTranslatef(-di1.getCentroid().x, di1.getCentroid().y, di1.getCentroid().z+1);

      //if(wires) di1.glRender();
      //di2.glRender();

      glBegin(GL_TRIANGLES);
      for(int i=0;i<mesh.size();i++){
    	  TRIANGLE t=mesh[i];
    	  Point3f color=colors[i];
          if(wires) glColor4f(color.x,color.y,color.z,0.0);
          //glNormal3f(t.n[0].x,t.n[0].y,t.n[0].z);
    	  glVertex3f(t.p[0].x,t.p[0].y,t.p[0].z);
          //glColor4f(0,1,0,0.1);
          //glNormal3f(t.n[1].x,t.n[1].y,t.n[1].z);
    	  glVertex3f(t.p[1].x,t.p[1].y,t.p[1].z);
          //glColor4f(0,0,1,0.1);
          //glNormal3f(t.n[2].x,t.n[2].y,t.n[2].z);
    	  glVertex3f(t.p[2].x,t.p[2].y,t.p[2].z);
      }
      glEnd();

//      glColor3f(1,1,0);
//	  glBegin(GL_POINTS);
//      for(Point3f p:vpts)
//    	  glVertex3f(p.x,-p.y,-p.z);
//      glEnd();
    glPopMatrix();
    glutSwapBuffers();
    angle++;
}

void init (void) {
    glEnable (GL_DEPTH_TEST);
    //glDisable(GL_CULL_FACE);
    glEnable (GL_LIGHTING);
    glEnable(GL_COLOR_MATERIAL);
    glEnable (GL_LIGHT0);
    glLightfv(GL_LIGHT0, GL_AMBIENT, light0_ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
    glClearColor (0.0,0.0,0.0,0.0);
}

void reshape(int width, int height)
{
    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60.0f, (GLfloat)width/(GLfloat)height, 0.5f, 200.0f);

    glMatrixMode(GL_MODELVIEW);
    ancho = width;
    alto = height;
}
void idle()
{
    displayMe();
    usleep(100000);
}
void keyPressed (unsigned char key, int x, int y) {
	x++;
	y++;
	//printf("%c %d,",key,key);
    switch(key)
    {
    case 'z':
    case 'Z':
      t+=0.1;
      break;
    case 'x':
    case 'X':
      t-=0.1;
      break;
    case 'p':
    case 'P':
      yaw++;
      break;

    case 'o':
    case 'O':
      yaw--;
      break;
    case 'q':
    case 'Q':
      pitch++;
      break;

    case 'a':
    case 'A':
      pitch--;
      break;
    case 'w':
    case 'W':
      roll++;
      break;

    case 's':
    case 'S':
      roll--;
      break;
    case 'v':
    case 'V':
      wires=true;
      break;
    case 'b':
    case 'B':
      wires=false;
      break;
    case 'f':
    case 'F':
      friccion=true;
      break;
    case 'g':
    case 'G':
      friccion=false;
      break;

    case 27:   // escape
      exit(0);
      break;
    }
}
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
void buildMesh(GridOctree<TsdfVoxel> &g,vector<TRIANGLE> &mesh){
	TRIANGLE triangles[10];
    GRIDCELL grid;
//    for(int i=0;i<t.getSize()-1;i++){
//    	cout << i << endl;
//    	for(int j=0 ;j<t.getSize()-1;j++)
//    		for(int k=0;k<t.getSize()-1;k++)	{
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
    			for (int l=0;l<n;l++){
    				computeNormals(&triangles[l]);
    				mesh.push_back(triangles[l]);
    				Point3f color;
    				color.x=vxl.r/255.0;
    				color.y=vxl.g/255.0;
    				color.z=vxl.b/255.0;
    				colors.push_back(color);
    			}
    	//	}
    }
}
float truncationDistance(float d){
	return ((int)d+1)*1.5;
}
// 30/7/2016 adding tau=truncation Distance
void updateGrid(GridOctree<TsdfVoxel> &g,DepthImage &di1){
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
				//if(d>1) continue;
				float tau=sz*truncationDistance(d);
				float tau2=tau*2;
				float itau=sz;
				//cout << "d   =" << d << endl;
//				cout << "sz  =" << sz << endl;
//				cout << "tau=" << tau << endl;
//				cout << "itau=" << itau  <<endl;
    			for(float dt=-tau;dt<=tau;dt+=itau){
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
    				vxl.wd=1/tau2/(vxl.z*vxl.z);
    				vxl.wr=1/tau2/(vxl.z*vxl.z);
    				vxl.wg=1/tau2/(vxl.z*vxl.z);
    				vxl.wb=1/tau2/(vxl.z*vxl.z);
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
//void rayMarching(g,t){
//	for(int u=0;u<640;u++){
//		for(int v=0;v<480;v++){
//
//		}
//	}
//}
#define sqr(x) ((x)*(x))
int main(int argc, char** argv)
{
	string basepath;
	if ( argc != 2 )
		{
	        printf("usage: cRGB pathToRGBDdataset \n");
	        return -1;
	}
	basepath=argv[1];
	int posI=0;
	DepthImage dImg1(basepath,posI);
	DepthImage dImg2(basepath,posI+1);
    di1=dImg1;//.sparse();
    di1.bilateralDepthFilter();
    cout << "R" << di1.getR() << endl;
    cout << "t" << di1.getT() << endl;
    di2=dImg2.sparse();
    //di=dImg1;
    cout << di1.getCentroid()<< " centroid"<<endl;
    cout << di1.getPoints3D().size()/1000 << "mil filtered points" <<endl;
    //vector<Point3f> pts=di1.getPoints3DCentered();
    //cout << "pts.size()" << pts.size() <<endl;

    //g.clear(1e32);
	g.setLevel(level);
    g.setMinMax(-3.0,3.0,
    		    -3.0,3.0,
				-3.0,3.0);
    //for(Point3f p:pts){
    //	t.setVoxel(p.x,p.y,p.z,0.0);
    //}
    for(int i=1;i<400;i+=1){
    	cout<<"i="<<i<<endl;
        cout <<"voxels="<< g.getVoxelsIdx().size() <<endl;
    	di1=DepthImage(basepath,i);
        di1.bilateralDepthFilter();
    	updateGrid(g,di1);
    }

    //build sphere and projectiveDistance
//    for(int i=0;i<g.getSizeX();i++){
//    	cout << i << endl;
//    	for(int j=0;j<g.getSizeY();j++)
//    		for(int k=0;k<g.getSizeZ();k++){
//    			float x=g.i2X(i);
//    			float y=g.j2Y(j);
//    			float z=g.k2Z(k);
//    			float p=0.20;
//    			float d0=sqrt(sqr(x)+sqr(y)+sqr(z))-0.25;
//    			float d1=sqrt(sqr(x-p)+sqr(y-p)+sqr(z-p))-0.125;
//    			float d2=sqrt(sqr(x+p)+sqr(y-p)+sqr(z-p))-0.125;
//    			float db=fmin(fmin(d0,d1),d2);
//    			Point3f p3D=Point3f(x,y,z);
//    			float pd=di1.projectiveDistance(Point3f(x,y,z));
//    			float d=fmin(db,pd);
//    			//float d=d0+d1;
//    			if(abs(d)<0.035){
//    				TsdfVoxel vxl;
//        			Point2f uv=di1.project(p3D);
//        			Point3f rp3D=di1.getPoint3D(uv);
//        			Vec3b c=di1.getColor(uv);
//        			vxl.x=rp3D.x;
//        			vxl.y=rp3D.y;
//        			vxl.z=rp3D.z;
//        			vxl.r=c[2];
//        			vxl.g=c[1];
//        			vxl.b=c[0];
//    				vxl.d=d;
//    				g.setVoxel(i,j,k,vxl);
//    				//cout << d<<":"<< g.getVoxel(i,j,k)<<endl;
//    			}
//    		}
//    }
    cout << g.getVoxelsIdx().size() <<endl;
//    cout << g.getChildrenPos(0b100,0b010,0b100,2) << endl;
//    cout << g.getChildrenPos(0b100,0b010,0b101,1) << endl;
//    GridOctree<float> goct;
//    float v=5;
//    goct.setVoxel(0b100,0b010,0b101,v);
//    float v1=8;
//    goct.setVoxel(0b100,0b010,0b100,v1);
//    cout << goct.getVoxel(0b100,0b010,0b100) << endl;
//    cout << goct.getVoxel(0b100,0b010,0b101) << endl;
//    cout << goct.getVoxel(0b100,0b011,0b101) << endl;
//    cout << *goct.getVoxelPtr(0b100,0b010,0b100) << endl;
    buildMesh(g,mesh);
    glutInit(&argc, argv);
    glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
    //glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB);
    // Enable blending
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glutInitWindowSize(700, 700);
    glutInitWindowPosition(250, 250);
    glutCreateWindow("cTSDF");
    init();
    glutDisplayFunc(displayMe);
    glutIdleFunc(idle);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyPressed); // Tell GLUT to use the method "keyPressed" for key presses
    glutMainLoop();
    return 0;
}

