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
Tsdf<float> t;
vector<Point3f> vpts;
vector<TRIANGLE> tiangles;

bool wires=true;
bool friccion=true;

void displayMe(void)
{
    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    gluLookAt (0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    GLfloat lightpos[] = {5.0, 15., 5., 0.};
    glLightfv(GL_LIGHT0, GL_POSITION, lightpos);
//    glTranslatef(0.0f, 0.0f, -3.0f);
    glRotatef(yaw  ,0.0,1.0,0.0);
    glRotatef(pitch,1.0,0.0,0.0);
    glRotatef(roll ,0.0,0.0,1.0);

    glPushMatrix();
      glTranslatef(-di1.getCentroid().x, di1.getCentroid().y, di1.getCentroid().z);
      di1.glRender();
      //di2.glRender();
      glColor3f(1,1,0);
	  glBegin(GL_POINTS);
      for(Point3f p:vpts)
    	  glVertex3f(p.x,-p.y,-p.z);
      glEnd();
    glPopMatrix();
    glutSwapBuffers();
    angle++;
}

void init (void) {
    glEnable (GL_DEPTH_TEST);
    //glEnable (GL_LIGHTING);
    //glEnable (GL_LIGHT0);
    //glEnable(GL_COLOR_MATERIAL);
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

int main(int argc, char** argv)
{
	string basepath;
	if ( argc != 2 )
		{
	        printf("usage: cRGB pathToRGBDdataset \n");
	        return -1;
	}
	basepath=argv[1];
	int posI=60;
	DepthImage dImg1(basepath,posI);
	DepthImage dImg2(basepath,posI+1);
    di1=dImg1;//.sparse();
    di2=dImg2.sparse();
    //di=dImg1;
    cout << di1.getCentroid()<< " centroid"<<endl;
    cout << di1.getPoints3D().size()/1000 << "mil filtered points" <<endl;
    t.clear(0.0);
    t.setMinMax(-2,2);
    vector<Point3f> pts=di1.getPoints3D();
    cout << "pts.size()" << pts.size() <<endl;
    for(Point3f p:pts){
    	t.setVoxel(p.x,p.y,p.z,1.0);
    }
    for(int i=0;i<t.getSize();i++)
    	for(int j=0;j<t.getSize();j++)
    		for(int k=0;k<t.getSize();k++)
    			if(t.getVoxel(i,j,k)==1.0){
    				float x=t.i2f(i);
    				float y=t.i2f(j);
    				float z=t.i2f(k);
    				Point3f p(x,y,z);
    				vpts.push_back(p);
    			}
    cout << vpts.size() <<endl;
	TRIANGLE triangles[10];
	TRIANGLE *tri = NULL;
	int ntri = 0;
    GRIDCELL grid;
    for(int i=0;i<t.getSize()-1;i++)
    	for(int j=0 ;j<t.getSize()-1;j++)
    		for(int k=0;k<t.getSize()-1;k++){
    			grid.p[0].x = t.i2f(i);
    			grid.p[0].y = t.i2f(j);
    			grid.p[0].z = t.i2f(k);
    			grid.val[0] = t.getVoxel(i,j,k);
                grid.p[1].x = t.i2f(i+1);
                grid.p[1].y = t.i2f(j);
                grid.p[1].z = t.i2f(k);
    			grid.val[1] =  t.getVoxel(i+1,j,k);
                grid.p[2].x = t.i2f(i+1);
                grid.p[2].y = t.i2f(j+1);
                grid.p[2].z = t.i2f(k);
    			grid.val[2] =  t.getVoxel(i+1,j+1,k);
                grid.p[3].x = t.i2f(i);
                grid.p[3].y = t.i2f(j+1);
                grid.p[3].z = t.i2f(k);
    			grid.val[3] =  t.getVoxel(i,j+1,k);
                grid.p[4].x = t.i2f(i);
                grid.p[4].y = t.i2f(j);
                grid.p[4].z = t.i2f(k+1);
    			grid.val[4] =  t.getVoxel(i,j,k+1);
                grid.p[5].x = t.i2f(i+1);
                grid.p[5].y = t.i2f(j);
                grid.p[5].z = t.i2f(k+1);
    			grid.val[5] =  t.getVoxel(i+1,j,k+1);
                grid.p[6].x = t.i2f(i+1);
                grid.p[6].y = t.i2f(j+1);
                grid.p[6].z = t.i2f(k+1);
    			grid.val[6] =  t.getVoxel(i+1,j+1,k+1);
                grid.p[7].x = t.i2f(i);
                grid.p[7].y = t.i2f(j+1);
                grid.p[7].z = t.i2f(k+1);
    			grid.val[7] =  t.getVoxel(i,j+1,k+1);
    			int	n = PolygoniseCube(grid,0.5,triangles);
    			tri = (TRIANGLE *)realloc(tri,(ntri+n)*sizeof(TRIANGLE));
    			for (int l=0;l<n;l++)
    				tri[ntri+l] = triangles[l];
    			ntri += n;
    			//triangles.push_back(p);
    			}

    glutInit(&argc, argv);
    glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
    //glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB);
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

