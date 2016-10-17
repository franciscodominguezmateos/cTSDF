//============================================================================
// Name        : cTSDF.cpp
// Author      : Francico Dominguez
// Version     :
// Copyright   : Your copyright notice
// Description : easy TSDF in C++, Ansi-style
//============================================================================
#include "TSDFoctGrid.h""

using namespace std;

typedef float S;

GLint ancho=400;
GLint alto=400;
int hazPerspectiva = 1;
GLfloat angle = 0.0;
GLfloat yaw = 0.0;
GLfloat roll = 0.0;
GLfloat pitch = 0.0;
GLfloat t=-3.0f;

DepthImage di1,di2;

bool wires=true;
bool boxes=false;
bool friccion=true;
int i=0;
string basepath;

GLfloat light0_ambient[] ={0.2, 0.2, 0.2, 1.0};
GLfloat light0_diffuse[] ={0.8, 0.8, 0.8, 1.0};

TSDFoctGrid tog;
GridOctree<TsdfVoxel> &g=tog.g;

void displayMe(void)
{
    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
//    gluLookAt (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    GLfloat lightpos[] = {3.0, 3.0, 3.0, 0.0};
    glTranslatef(0.0f, 0.0f, t);
    glRotatef(yaw  ,0.0,1.0,0.0);
    glRotatef(pitch,1.0,0.0,0.0);
    glRotatef(roll ,0.0,0.0,1.0);
    glLightfv(GL_LIGHT0, GL_POSITION, lightpos);
    glPushMatrix();
    tog.glDrawMesh();
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
    cv::waitKey(1);
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
      tog.wires=true;
      break;
    case 'b':
    case 'B':
      tog.wires=false;
      break;
    case 'f':
    case 'F':
      friccion=true;
      break;
    case 'g':
    case 'G':
      friccion=false;
      break;
    case 'c':
    case 'C':
      tog.boxes=!tog.boxes;
      break;
    case 'i':
    case 'I':
    	cout<<"i="<<i<<endl;
    	di1=DepthImage(basepath,i);
        di1.bilateralDepthFilter();
    	tog.updateGrid(di1);
        cout <<"voxels="<< g.getVoxelsIdx().size() <<endl;
        tog.mesh.clear();
        tog.colors.clear();
        tog.buildMesh();
        cout <<"mesh="<< tog.mesh.size() <<endl;
        i++;
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
    for(;i<1;i+=1){
    	cout<<"i="<<i<<endl;
    	di1=DepthImage(basepath,i);
        di1.bilateralDepthFilter();
    	tog.updateGrid(di1);
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

