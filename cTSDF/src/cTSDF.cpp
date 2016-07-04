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
#include "tsdf.h"
#include "grid.h"
#include "grid_octree.h"
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

DepthImage di1,di2;
Grid<float> g(256/2,256/2,256/2);
vector<Point3f> vpts;
vector<TRIANGLE> mesh;

bool wires=true;
bool friccion=true;

GLfloat light0_ambient[] ={0.2, 0.2, 0.2, 1.0};
GLfloat light0_diffuse[] ={0.8, 0.8, 0.8, 1.0};

void displayMe(void)
{
    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    gluLookAt (0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    GLfloat lightpos[] = {0.0, 15., 5., 0.0};
    glLightfv(GL_LIGHT0, GL_POSITION, lightpos);
//    glTranslatef(0.0f, 0.0f, -3.0f);
    glRotatef(yaw  ,0.0,1.0,0.0);
    glRotatef(pitch,1.0,0.0,0.0);
    glRotatef(roll ,0.0,0.0,1.0);
    glPushMatrix();
     //glTranslatef(-di1.getCentroid().x, di1.getCentroid().y, di1.getCentroid().z+1);

     if(wires) di1.glRender();
      //di2.glRender();
      glBegin(GL_TRIANGLES);
      for(TRIANGLE t:mesh){
          glColor4f(1,0,0.0,0.1);
          glNormal3f(t.n[0].x,-t.n[0].y,-t.n[0].z);
    	  glVertex3f(t.p[0].x,-t.p[0].y,-t.p[0].z);
          //glColor4f(0,1,0,0.1);
          glNormal3f(t.n[1].x,-t.n[1].y,-t.n[1].z);
    	  glVertex3f(t.p[1].x,-t.p[1].y,-t.p[1].z);
          //glColor4f(0,0,1,0.1);
          glNormal3f(t.n[2].x,-t.n[2].y,-t.n[2].z);
    	  glVertex3f(t.p[2].x,-t.p[2].y,-t.p[2].z);
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
    gluPerspective(30.0f, (GLfloat)width/(GLfloat)height, 0.5f, 200.0f);

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
void buildMesh(Grid<float> &g,vector<TRIANGLE> &mesh){
	TRIANGLE triangles[10];
    GRIDCELL grid;
//    for(int i=0;i<t.getSize()-1;i++){
//    	cout << i << endl;
//    	for(int j=0 ;j<t.getSize()-1;j++)
//    		for(int k=0;k<t.getSize()-1;k++)	{
    int i,j,k;
    for(int idx:g.getVoxelsIdx()){
       	g.getIJKfromIdx(idx,i,j,k);
    			grid.p[0].x = g.i2X(i);
    			grid.p[0].y = g.j2Y(j);
    			grid.p[0].z = g.k2Z(k);
    			grid.val[0] = g.getVoxel(i,j,k);
                grid.p[1].x = g.i2X(i+1);
                grid.p[1].y = g.j2Y(j);
                grid.p[1].z = g.k2Z(k);
    			grid.val[1] =  g.getVoxel(i+1,j,k);
                grid.p[2].x = g.i2X(i+1);
                grid.p[2].y = g.j2Y(j+1);
                grid.p[2].z = g.k2Z(k);
    			grid.val[2] =  g.getVoxel(i+1,j+1,k);
                grid.p[3].x = g.i2X(i);
                grid.p[3].y = g.j2Y(j+1);
                grid.p[3].z = g.k2Z(k);
    			grid.val[3] =  g.getVoxel(i,j+1,k);
                grid.p[4].x = g.i2X(i);
                grid.p[4].y = g.j2Y(j);
                grid.p[4].z = g.k2Z(k+1);
    			grid.val[4] =  g.getVoxel(i,j,k+1);
                grid.p[5].x = g.i2X(i+1);
                grid.p[5].y = g.j2Y(j);
                grid.p[5].z = g.k2Z(k+1);
    			grid.val[5] =  g.getVoxel(i+1,j,k+1);
                grid.p[6].x = g.i2X(i+1);
                grid.p[6].y = g.j2Y(j+1);
                grid.p[6].z = g.k2Z(k+1);
    			grid.val[6] =  g.getVoxel(i+1,j+1,k+1);
                grid.p[7].x = g.i2X(i);
                grid.p[7].y = g.j2Y(j+1);
                grid.p[7].z = g.k2Z(k+1);
    			grid.val[7] =  g.getVoxel(i,j+1,k+1);
    			int	n = PolygoniseCube(grid,0.0,triangles);
    			for (int l=0;l<n;l++){
    				computeNormals(&triangles[l]);
    				mesh.push_back(triangles[l]);
    			}
    	//	}
    }
}
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
    di2=dImg2.sparse();
    //di=dImg1;
    cout << di1.getCentroid()<< " centroid"<<endl;
    cout << di1.getPoints3D().size()/1000 << "mil filtered points" <<endl;
    //vector<Point3f> pts=di1.getPoints3DCentered();
    //cout << "pts.size()" << pts.size() <<endl;

    g.clear(1e32);
    //t.setMinMax(-0.35,0.35);
    g.setMinMax(-1.0,1.0,
    		    -1.0,1.0,
				 0.9,2.0);
    //for(Point3f p:pts){
    //	t.setVoxel(p.x,p.y,p.z,0.0);
    //}
    //build sphere and projectiveDistance
    for(int i=0;i<g.getSizeX();i++){
    	cout << i << endl;
    	for(int j=0;j<g.getSizeY();j++)
    		for(int k=0;k<g.getSizeZ();k++){
    			float x=g.i2X(i);
    			float y=g.j2Y(j);
    			float z=g.k2Z(k);
    			float p=0.20;
    			float d0=sqrt(sqr(x)+sqr(y)+sqr(z))-0.25;
    			float d1=sqrt(sqr(x-p)+sqr(y-p)+sqr(z-p))-0.125;
    			float d2=sqrt(sqr(x+p)+sqr(y-p)+sqr(z-p))-0.125;
    			float db=fmin(fmin(d0,d1),d2);
    			float pd=di1.projectiveDistance(Point3f(x,y,z));
    			float d=fmin(db,pd);
    			//float d=d0+d1;
    			if(abs(d)<0.025)
    				g.setVoxel(i,j,k,d);
    		}
    }
    cout << g.getVoxelsIdx().size() <<endl;
    cout << g.getChildrenPos(0b100,0b010,0b100,2) << endl;
    cout << g.getChildrenPos(0b100,0b010,0b100,1) << endl;
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

