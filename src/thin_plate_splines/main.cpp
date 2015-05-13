/*
 *  Thin Plate Spline demo/example in C++
 *
 *  - a simple TPS editor, using the Boost uBlas library for large
 *    matrix operations and OpenGL + GLUT for 2D function visualization
 *    (curved plane) and user interface
 *
 *  Copyright (C) 2003,2005 by Jarno Elonen
 *
 *  TPSDemo is Free Software / Open Source with a very permissive
 *  license:
 *
 *  Permission to use, copy, modify, distribute and sell this software
 *  and its documentation for any purpose is hereby granted without fee,
 *  provided that the above copyright notice appear in all copies and
 *  that both that copyright notice and this permission notice appear
 *  in supporting documentation.  The authors make no representations
 *  about the suitability of this software for any purpose.
 *  It is provided "as is" without express or implied warranty.
 *
 *  TODO:
 *    - implement TPS approximation 3 as suggested in paper
 *      Gianluca Donato and Serge Belongie, 2002: "Approximation
 *      Methods for Thin Plate Spline Mappings and Principal Warps"
 */

#include <GL/glut.h>
#include <boost/numeric/ublas/matrix.hpp>

#include "linalg3d.h"
#include "ludecomposition.h"

#include <vector>
#include <cmath>
#include <cstdio>
#include <cstring>

using namespace boost::numeric::ublas;

// ========= BEGIN INTERESTING STUFF  =========

#define GRID_W 100
#define GRID_H 100
static float grid[GRID_W][GRID_H];

std::vector< Vec > control_points;
int selected_cp = -1;

double regularization = 0.0;
double bending_energy = 0.0;

static double tps_base_func(double r)
{
  if ( r == 0.0 )
    return 0.0;
  else
    return r*r * log(r);
}


/*
 *  Calculate Thin Plate Spline (TPS) weights from
 *  control points and build a new height grid by
 *  interpolating with them.
 */
static void calc_tps()
{
  // You We need at least 3 points to define a plane
  if ( control_points.size() < 3 )
    return;

  unsigned p = control_points.size();

  // Allocate the matrix and vector
  matrix<double> mtx_l(p+3, p+3);
  matrix<double> mtx_v(p+3, 1);
  matrix<double> mtx_orig_k(p, p);

  // Fill K (p x p, upper left of L) and calculate
  // mean edge length from control points
  //
  // K is symmetrical so we really have to
  // calculate only about half of the coefficients.
  double a = 0.0;
  for ( unsigned i=0; i<p; ++i )
  {
    for ( unsigned j=i+1; j<p; ++j )
    {
      Vec pt_i = control_points[i];
      Vec pt_j = control_points[j];
      pt_i.y = pt_j.y = 0;
      double elen = (pt_i - pt_j).len();
      mtx_l(i,j) = mtx_l(j,i) =
        mtx_orig_k(i,j) = mtx_orig_k(j,i) =
          tps_base_func(elen);
      a += elen * 2; // same for upper & lower tri
    }
  }
  a /= (double)(p*p);

  // Fill the rest of L
  for ( unsigned i=0; i<p; ++i )
  {
    // diagonal: reqularization parameters (lambda * a^2)
    mtx_l(i,i) = mtx_orig_k(i,i) =
      regularization * (a*a);

    // P (p x 3, upper right)
    mtx_l(i, p+0) = 1.0;
    mtx_l(i, p+1) = control_points[i].x;
    mtx_l(i, p+2) = control_points[i].z;

    // P transposed (3 x p, bottom left)
    mtx_l(p+0, i) = 1.0;
    mtx_l(p+1, i) = control_points[i].x;
    mtx_l(p+2, i) = control_points[i].z;
  }
  // O (3 x 3, lower right)
  for ( unsigned i=p; i<p+3; ++i )
    for ( unsigned j=p; j<p+3; ++j )
      mtx_l(i,j) = 0.0;


  // Fill the right hand vector V
  for ( unsigned i=0; i<p; ++i )
    mtx_v(i,0) = control_points[i].y;
  mtx_v(p+0, 0) = mtx_v(p+1, 0) = mtx_v(p+2, 0) = 0.0;

  // Solve the linear system "inplace"
  if (0 != LU_Solve(mtx_l, mtx_v))
  {
    puts( "Singular matrix! Aborting." );
    exit(1);
  }

  // Interpolate grid heights
  for ( int x=-GRID_W/2; x<GRID_W/2; ++x )
  {
    for ( int z=-GRID_H/2; z<GRID_H/2; ++z )
    {
      double h = mtx_v(p+0, 0) + mtx_v(p+1, 0)*x + mtx_v(p+2, 0)*z;
      Vec pt_i, pt_cur(x,0,z);
      for ( unsigned i=0; i<p; ++i )
      {
        pt_i = control_points[i];
        pt_i.y = 0;
        h += mtx_v(i,0) * tps_base_func( ( pt_i - pt_cur ).len());
      }
      grid[x+GRID_W/2][z+GRID_H/2] = h;
    }
  }

  // Calc bending energy
  matrix<double> w( p, 1 );
  for ( int i=0; i<p; ++i )
    w(i,0) = mtx_v(i,0);
  matrix<double> be = prod( prod<matrix<double> >( trans(w), mtx_orig_k ), w );
  bending_energy = be(0,0);
}

// ========= END INTERESTING STUFF  =========
// (The rest is essentially just visualization with OpenGL,
//  altough that can for sure be interesting, too.)


static int winW = 800, winH = 600;
static int mouseX = -999, mouseY = -999;
static bool mouseState[3] = {false};
static int modifiers = 0;
Vec cursor_loc;
static float camAlpha=30, camBeta=5, camZoom=0;
static bool screen_dirty=true;

static void clear_grid()
{
  for (int x=0; x<GRID_W; ++x)
    for (int z=0; z<GRID_H; ++z)
      grid[x][z] = 0;
}

#define SQUARE(x) ((x)*(x))

void draw_string (const char* str)
{
  for (unsigned i=0; i<strlen(str); i++)
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, str[i]);
};

// OGL: draw callback
static void display()
{
  unsigned i;

  // Uploads a Vec to OGL
  #define UPLOADVEC(v) glVertex3f( v.x, v.y, v.z )

  static GLfloat color1[] = {0.8, 0.8, 0.8, 1.0};
  static GLfloat color2[] = {1.0, 1.0, 1.0, 1.0};

  static GLfloat red[] = {1.0, 0.0, 0.0, 1.0};
  static GLfloat green[] = {0.0, 1.0, 0.0, 1.0};
  static GLfloat blue[] = {0.0, 0.0, 1.0, 1.0};

  // Make a rotation matrix out of mouse point

  const Mtx& rot = rotateY( camBeta ) * rotateX( camAlpha );

  // Rotate camera
  Vec cam_loc(0,0,-150), cam_up(0,1,0);
  cam_loc = cam_loc * rot;
  cam_up = cam_up * rot;

  // Clear the screen
  glClear( GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT );

  // Prepare zoom by changiqng FOV
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  float fov = 45 + camZoom;
  if ( fov < 5 ) fov = 5;
  if ( fov > 160 ) fov = 160;
  gluPerspective( fov, (float)winW/(float)winH, 1.0, 500.0 );

  gluLookAt( cam_loc.x, cam_loc.y, cam_loc.z,  // eye
    0, 0, 0, // target
    cam_up.x, cam_up.y, cam_up.z ); // up

  // Curve surface
  glEnable(GL_LIGHTING) ;
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glBegin( GL_QUADS );

  static GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
  static GLfloat mat_shininess[] = { 100.0 };

  for ( int x=-GRID_W/2; x<GRID_W/2-1; ++x )
  {
    for ( int z=-GRID_H/2; z<GRID_H/2-1; ++z )
    {
      if ( (x&8)^(z&8) )
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, color1);
      else
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, color2);

      glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
      glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);

      float a[] = { x+0, grid[x+0+GRID_W/2][z+0+GRID_H/2], z+0 };
      float b[] = { x+0, grid[x+0+GRID_W/2][z+1+GRID_H/2], z+1 };
      float c[] = { x+1, grid[x+1+GRID_W/2][z+1+GRID_H/2], z+1 };
      float d[] = { x+1, grid[x+1+GRID_W/2][z+0+GRID_H/2], z+0 };

      #define V_MINUS(A,B) {A[0]-B[0], A[1]-B[1], A[2]-B[2]}
      #define V_CROSS(A,B) \
        {A[1]*B[2]-A[2]*B[1], \
         A[2]*B[0]-A[0]*B[2], \
         A[0]*B[1]-A[1]*B[0]}
      float ab[] = V_MINUS(a,b);
      float cb[] = V_MINUS(c,b);
      float n[] = V_CROSS( cb, ab );

      glNormal3f( n[0],n[1],n[2] );
      glVertex3f( a[0],a[1],a[2] );
      glVertex3f( b[0],b[1],b[2] );
      glVertex3f( c[0],c[1],c[2] );
      glVertex3f( d[0],d[1],d[2] );
    }
  }
  glEnd();

#if 1 // visual helpers

  glDisable(GL_LIGHTING) ;

  // Flat grid
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  glBegin( GL_QUADS );
  glColor3ub( 128, 128, 128 );

  for ( int x=-GRID_W/2; x<GRID_W/2-4; x+=5 )
  {
    for ( int z=-GRID_H/2; z<GRID_H/2-4; z+=5 )
    {
      glVertex3f( x-0.5, -0.5f, z-0.5 );
      glVertex3f( x-0.5, -0.5f, z+4.5 );
      glVertex3f( x+4.5, -0.5f, z+4.5 );
      glVertex3f( x+4.5, -0.5f, z-0.5 );
    }
  }
  glEnd();

  // Axes: (x,y,z)=(r,g,b)
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glPushMatrix();
  glRotatef(90.0, 0.0, 1.0, 0.0);
  glColor3ub( 255, 0, 0 );
  glutSolidCone(0.5, 80.0, 5, 1);
  glPopMatrix();
  glPushMatrix();
  glRotatef(-90.0, 1.0, 0.0, 0.0);
  glColor3ub( 0, 255,  0 );
  glutSolidCone(0.5, 80.0, 5, 1);
  glPopMatrix();
  glPushMatrix();
  glRotatef(0.0, 1.0, 0.0, 0.0);
  glColor3ub( 0, 0, 255 );
  glutSolidCone(0.5, 80.0, 5, 1);
  glPopMatrix();

#endif

  // Control points
  int old_sel = selected_cp;
  if ( mouseState[0] == mouseState[1] == mouseState[2] == 0 )
    selected_cp = -1;

  for ( int i=0; i < control_points.size(); ++i )
  {
    const Vec& cp = control_points[i];
    if ( ( cp - cursor_loc ).len() < 2.0 )
    {
      selected_cp = i;
      glutSetCursor( GLUT_CURSOR_UP_DOWN );
    }

    glPushMatrix();
    glTranslatef(cp.x, cp.y, cp.z);
    if ( selected_cp == i )
      glColor3ub( 0, 255, 255 );
    else
      glColor3ub( 255, 255, 0 );
    glutSolidSphere(1.0,12,12);
    glPopMatrix();

    glBegin( GL_LINES );
    glVertex3f( cp.x, 0, cp.z );
    glVertex3f( cp.x, cp.y, cp.z );
    glEnd();
  }

  if ( selected_cp < 0 && old_sel != selected_cp )
    glutSetCursor( GLUT_CURSOR_CROSSHAIR );

	// Find out the world coordinates of mouse pointer
	// to locate the cursor
	if ( mouseState[0] == mouseState[1] == mouseState[2] == 0 )
	{
		GLdouble model[16], proj[16];
		GLint view[4];
		GLfloat z;
		GLdouble ox, oy, oz;

		glGetDoublev(GL_MODELVIEW_MATRIX, model);
		glGetDoublev(GL_PROJECTION_MATRIX, proj);
		glGetIntegerv(GL_VIEWPORT, view);

		glReadPixels(mouseX, view[3]-mouseY-1, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &z);
		gluUnProject(mouseX, view[3]-mouseY-1, z, model, proj, view, &ox, &oy, &oz);

    cursor_loc = Vec(ox, oy, oz);

    // Draw the cursor
    glPushMatrix();
    glDisable(GL_LIGHTING) ;
    glTranslatef(ox, oy, oz);
    glColor3ub( 255, 0, 255 );
    glutSolidSphere(0.5,12,12);
    glPopMatrix();
	}

  static char tmp_str[255];
  glDisable( GL_DEPTH_TEST );
  glLoadIdentity();
  gluOrtho2D(-1.0, 1.0, -1.0, 1.0);
  glColor3ub( 255, 255, 0 );
  glRasterPos2f (-0.95, -0.95);
  sprintf( tmp_str, "control points: %d, reqularization: %2.3f, bending energy: %4.3f",
    (int)control_points.size(), regularization, bending_energy );
  draw_string( tmp_str );
  glEnable( GL_DEPTH_TEST );

  glFlush();
	glutSwapBuffers();
  screen_dirty = false;
}

// OGL: window resize callback
static void reshape( int w, int h )
{
	winW = w;
	winH = h;
	glViewport( 0, 0, winW, winH );
	glEnable( GL_DEPTH_TEST );
	glDisable( GL_CULL_FACE );
  glEnable( GL_NORMALIZE );
  glDepthFunc( GL_LESS );
}

// OGL: key press callback
static void keyboard( unsigned char key, int, int )
{
  switch (key)
  {
    case 'a':
      control_points.push_back( cursor_loc );
      calc_tps();
      break;
    case 'd':
      if ( selected_cp >= 0 )
      {
        control_points.erase( control_points.begin() + selected_cp );
        selected_cp = -1;
        calc_tps();
      }
      break;
    case 'c':
      control_points.clear();
      clear_grid();
      break;
    case '+':
      regularization += 0.025;
      calc_tps();
      break;
    case '-':
      regularization -= 0.025;
      if (regularization < 0) regularization = 0;
      calc_tps();
      break;
    case '/': camZoom -= 1; break;
    case '*': camZoom += 1; break;
    case 'q': exit( 0 ); break;
  }
  screen_dirty=true;
}

// OGL: mouse button callback
static void mouse( int button, int state, int, int )
{
	mouseState[ button ] = (state==GLUT_DOWN);
	modifiers = glutGetModifiers();

  glutSetCursor( GLUT_CURSOR_CROSSHAIR );

  if ( button == 1 && state==GLUT_DOWN )
      glutSetCursor( GLUT_CURSOR_CYCLE );

	if ( button == 0 )
  {
    if ( state==GLUT_UP )
    {
		  calc_tps();
      screen_dirty=true;
    }
    else if ( state==GLUT_DOWN && selected_cp<0 )
      keyboard( 'a', 0,0 );
  }
}

// OGL: mouse movement callback
static void mouseMotion( int x, int y )
{
	if ( mouseState[0] && mouseX != -999 )
    if ( selected_cp >= 0 )
			 control_points[selected_cp].y += -(y - mouseY)/3;

	if ( mouseState[1] && mouseX != -999 )
	{
    camAlpha += -(y - mouseY);
    camBeta += (x - mouseX);

    screen_dirty=true;
	}

  if ( mouseX != x || mouseY != y )
  {
    mouseX = x;
    mouseY = y;
    screen_dirty=true;
  }
}

// OGL: keyboard press callback for special characters
static void keyboard_special( int key, int, int )
{
  switch (key)
  {
    case GLUT_KEY_UP: camAlpha += 5 ; break;
    case GLUT_KEY_DOWN: camAlpha += -5 ; break;
    case GLUT_KEY_RIGHT: camBeta += -5 ; break;
    case GLUT_KEY_LEFT: camBeta += 5 ; break;
  }
  screen_dirty=true;
}

// OGL: menu selection callback
static void menu_select(int mode)
{
  keyboard( (unsigned char)mode, 0,0 );
}

static void create_menu(void)
{
  glutCreateMenu(menu_select);
  glutAddMenuEntry(" d     Delete control point",'d');
  glutAddMenuEntry(" a     Add control point",'a');
  glutAddMenuEntry(" c     Clear all",'c');
  glutAddMenuEntry(" +     Relax more",'+');
  glutAddMenuEntry(" -     Rela less",'-');
  glutAddMenuEntry(" /     Zoom in",'/');
  glutAddMenuEntry(" *     Zoom out",'*');
  glutAddMenuEntry(" q     Exit",'q');
  glutAttachMenu(GLUT_RIGHT_BUTTON);
}

// Print out OpenGL errors, if any
static void glCheckErrors()
{
	GLenum errCode = glGetError();
	if ( errCode != GL_NO_ERROR )
	{
		const GLubyte *errString = gluErrorString( errCode );
		fprintf(stderr, "OpenGL error: %s\n", errString);
	}
}

// OGL: idle callback, continuosly keep redrawing
static void idlefunc()
{
  glCheckErrors();
  if (screen_dirty)
    glutPostRedisplay();
}

// Startup
int main( int argc, char *argv[] )
{
  glutInit( &argc, argv );
  glutInitDisplayMode( GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH );
  glutInitWindowSize( winW, winH );
  glutInitWindowPosition( 0, 0 );
  if ( !glutCreateWindow( "tpldemo" ) )
  {
    printf( "Couldn't open window.\n" );
    return 1;
  }
  glutDisplayFunc( display );
  glutIdleFunc( idlefunc );
  glutMouseFunc( mouse );
  glutMotionFunc( mouseMotion );
  glutPassiveMotionFunc( mouseMotion );
  glutKeyboardFunc( keyboard );
  glutSpecialFunc( keyboard_special );
  glutReshapeFunc( reshape );
  glutSetCursor( GLUT_CURSOR_CROSSHAIR );

  glEnable(GL_LIGHTING);
  GLfloat lightAmbient[] = {0.5f, 0.5f, 0.5f, 1.0f} ;
  glLightModelfv(GL_LIGHT_MODEL_TWO_SIDE, lightAmbient);

  glEnable(GL_LIGHT0);
  GLfloat light0Position[] = {0.7f, 0.5f, 0.9f, 0.0f} ;
  GLfloat light0Ambient[]  = {0.0f, 0.5f, 1.0f, 0.8f} ;
  GLfloat light0Diffuse[]  = {0.0f, 0.5f, 1.0f, 0.8f} ;
  GLfloat light0Specular[] = {0.0f, 0.5f, 1.0f, 0.8f} ;
  glLightfv(GL_LIGHT0, GL_POSITION,light0Position) ;
  glLightfv(GL_LIGHT0, GL_AMBIENT,light0Ambient) ;
  glLightfv(GL_LIGHT0, GL_DIFFUSE,light0Diffuse) ;
  glLightfv(GL_LIGHT0, GL_SPECULAR,light0Specular) ;

  glEnable(GL_LIGHT1);
  GLfloat light1Position[] = {0.5f, 0.7f, 0.2f, 0.0f} ;
  GLfloat light1Ambient[]  = {1.0f, 0.5f, 0.0f, 0.8f} ;
  GLfloat light1Diffuse[]  = {1.0f, 0.5f, 0.0f, 0.8f} ;
  GLfloat light1Specular[] = {1.0f, 0.5f, 0.0f, 0.8f} ;
  glLightfv(GL_LIGHT1, GL_POSITION,light1Position) ;
  glLightfv(GL_LIGHT1, GL_AMBIENT,light1Ambient) ;
  glLightfv(GL_LIGHT1, GL_DIFFUSE,light1Diffuse) ;
  glLightfv(GL_LIGHT1, GL_SPECULAR,light1Specular) ;

  glCheckErrors();
  clear_grid();

  create_menu() ;
  glutMainLoop();

  return 0;
}
