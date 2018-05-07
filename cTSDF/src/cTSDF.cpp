//============================================================================
// Name        : cTSDF.cpp
// Author      : Francico Dominguez
// Version     :
// Copyright   : Your copyright notice
// Description : easy TSDF in C++, Ansi-style
//============================================================================
#include "TSDFoctGrid.h"

using namespace std;

typedef float S;

GLint ancho=400;
GLint alto=400;
int hazPerspectiva = 1;
GLfloat tx=0.8,ty=1.4,tz=0.8;
GLfloat angle = 0.0;
GLfloat yaw = 0.0;
GLfloat roll = 0.0;
GLfloat pitch = 0.0;
GLfloat t=-3.0f;
int mx=-1,my=-1;        // Prevous mouse coordinates
int rotangles[2] = {0}; // Panning angles
float zoom = 1;         // zoom factor

DepthImage di1,di2;

bool wires=true;
bool boxes=false;
bool friccion=true;
bool digl=false;
int i=0;
string basepath;

GLfloat light0_ambient[] ={0.2, 0.2, 0.2, 1.0};
GLfloat light0_diffuse[] ={0.8, 0.8, 0.8, 1.0};

TSDFoctGrid tog(9);
GridOctree<TsdfVoxel> &g=tog.g;

void displayMe(void)
{
    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
//    gluLookAt (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    GLfloat lightpos[] = {30.0, 30.0, 30.0, 0.0};
    glScalef(zoom,zoom,zoom);
    glRotatef(-90,1,0,0);//put Z up
    glTranslatef(0,0,0);
    glTranslatef(0.0f, -t, 0);
    glRotatef(rotangles[0], 1,0,0);
    glRotatef(rotangles[1], 0,0,1);
    //Axes
    int l=10;
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
    glTranslatef(-tx,-ty,-tz);
    glRotatef(yaw  ,0.0,1.0,0.0);
    glRotatef(pitch,1.0,0.0,0.0);
    glRotatef(roll ,0.0,0.0,1.0);
    glLightfv(GL_LIGHT0, GL_POSITION, lightpos);
    glPushMatrix();
    tog.glDrawMesh();
    if(digl) di1.glRender();
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
    gluPerspective(60.0f, (GLfloat)width/(GLfloat)height, 0.2f, 200.0f);

    glMatrixMode(GL_MODELVIEW);
    ancho = width;
    alto = height;
}
void idle()
{
    displayMe();
    usleep(100000);
    cv::waitKey(1);
}
void mouseMoved(int x, int y)
{
    if (mx>=0 && my>=0) {
        rotangles[0] += (float)(y-my)/2.0;
        rotangles[1] += (float)(x-mx)*0.5;
    }
    mx = x;
    my = y;
}

void mousePress(int button, int state, int x, int y)
{
    if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
        mx = x;
        my = y;
    }
    if (button == GLUT_LEFT_BUTTON && state == GLUT_UP) {
        mx = -1;
        my = -1;
    }
}

void keyPressed (unsigned char key, int x, int y) {
	x++;
	y++;
	//printf("%c %d,",key,key);
    switch(key)
    {
    case 'd':
    case 'D':
      tog.glPoints=!tog.glPoints;
      break;
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
      //yaw++;
      tx+=0.1;
      break;

    case 'o':
    case 'O':
      //yaw--;
    	tx-=0.1;
      break;
    case 'q':
    case 'Q':
      //pitch++;
    	tz+=0.1;
      break;

    case 'a':
    case 'A':
      //pitch--;
    	tz-=0.1;
      break;
    case 'w':
    case 'W':
      //roll++;
      ty+=0.1;
      break;

    case 's':
    case 'S':
      //roll--;
      ty-=0.1;
      cout<<"tx="<<tx<<"ty="<<ty<<"tz="<<tz<<endl;
      break;
    case 'v':
    case 'V':
      tog.wires=true;
      break;
    case 'b':
    case 'B':
      tog.wires=false;
      break;
    case 'f':
    case 'F':
      tog.iBoxes+=0.01;
      break;
    case 'g':
    case 'G':
      tog.iBoxes-=0.01;
      break;
    case 'r':
    case 'R':
      tog.iBoxesW+=0.01;
      break;
    case 't':
    case 'T':
      tog.iBoxesW-=0.01;
      break;
    case 'c':
    case 'C':
      tog.boxes=!tog.boxes;
      break;
    case 'm':
    case 'M':
    	digl=!digl;
      break;
    case 'i':
    case 'I':
    	cout<<"i="<<i<<endl;
    	di1=DepthImage(basepath,i);
        di1.bilateralDepthFilter();
        di2=di1.pyrDown(0.5);
    	tog.updateGrid3(di2);
    	//imshow("di1",di1.getImg());
    	imshow("di2",di2.getImg());
    	imshow("nor",di2.getNormals()*255);
    	//imshow("ddw",di2.getDepth()/2.5);
    	//imshow("dd1",di1.getDepth()/2.5);
        cout <<"voxels="<< g.getVoxelsIdx().size() <<endl;
        tog.mesh.clear();
        tog.colors.clear();
        tog.buildMesh();
        cout <<"mesh="<< tog.mesh.size() <<endl;
        i+=1;
      break;
    case 27:   // escape
      exit(0);
      break;
    }
}
#define sqr(x) ((x)*(x))
int main(int argc, char** argv)
{
	if ( argc != 2 )
		{
	        printf("usage: cTSDF pathToRGBDdataset \n");
	        return -1;
	}
	basepath=argv[1];
	DepthImage di(basepath,0);
	di.computeGrad();
	di2=di.pyrDown(0.5);
	imshow("di2",di2.getImg());
	/*
	Mat gx=di2.getGradXImg();
	double minGX,maxGX;
	minMaxLoc(gx,&minGX,&maxGX);
	float d=maxGX-minGX;
	gx=gx-minGX;
	gx=gx/d;
	imshow("gX",gx);
	*/

    for(i=0;i<1;i+=1){
    	cout<<"i="<<i<<endl;
    	di1=DepthImage(basepath,i);
        di1.bilateralDepthFilter();
        di2=di1.pyrDown(0.5);
    	tog.updateGrid3(di2);
        cout <<"voxels="<< g.getVoxelsIdx().size() <<endl;
    }
    //cout << "Comienzo RayMarching"<<endl;
    //di1=DepthImage(basepath,30);
    //rayMarching(g,di1);
    //cv::imshow("hola",di1.getDepth());
    //cout << "fin"<<endl;
    cout << g.getVoxelsIdx().size() <<endl;
    tog.buildMesh();
    cout << tog.mesh.size() <<endl;
    glutInit(&argc, argv);
    glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
    //glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB);
    // Enable blending
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glutInitWindowSize(700, 400);
    glutInitWindowPosition(250, 250);
    glutCreateWindow("cTSDF");
    init();
    glutDisplayFunc(displayMe);
    glutIdleFunc(idle);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyPressed); // Tell GLUT to use the method "keyPressed" for key presses
    glutMotionFunc(&mouseMoved);
    glutMouseFunc(&mousePress);
    glutMainLoop();
    return 0;
}

