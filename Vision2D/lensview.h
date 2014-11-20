#ifndef _LENVIEW_H_
#define _LENVIEW_H_
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

//#pragma comment( lib, "opengl32.lib" )
//#pragma comment( lib, "glu32.lib" )
#pragma comment( lib, "glut32.lib" )

#define BGR 1.0f
#define BGG 1.0f
#define BGB 1.0f

#define LVW 770
#define LVH 480
#define LVSD 20

#define BARC 0.75f
#define BART 0.0f
#define BARH 25
//#define SBFONT GLUT_BITMAP_9_BY_15
#define SBFONT GLUT_BITMAP_HELVETICA_18
#define X 0
#define Y 1
#define maxv(a,b) (((a) > (b))?(a):(b))

static inline void drawLine(float x1, float y1, 
							float x2, float y2) 
{ 
	glBegin(GL_LINES); 
	glVertex2f(x1, y1); glVertex2f(x2, y2); 
	glEnd(); 
}

static inline void drawPoint(float size, float x, float y) 
{ 
	glPointSize(size); 
	glBegin(GL_POINTS); 
	glVertex2f(x, y); 
	glEnd(); 
	glPointSize(1.0f); 
}

static inline void drawSquare(float r, float g, float b, 
							  GLint xl, GLint yl, 
							  GLint xr, GLint yr) 
{ 
	glBegin(GL_QUADS); 
	glColor3f(r, g, b); 
	glVertex2i(xl, yl); 
	glVertex2i(xl, yr); 
	glVertex2i(xr, yr); 
	glVertex2i(xr, yl); 
	glEnd(); 
}

static inline void drawText(float x, float y, void *font,
							char* text)
{
	char *p;
	glRasterPos2f(x, y);
	for (p = text; *p; ++p) {
		glutBitmapCharacter(font, *p);
	}
}

#endif