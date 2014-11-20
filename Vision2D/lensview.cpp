#include <cstdarg>
//#include "eyemodel.h"
#include "eye.h"
#include "lensview.h"

#pragma comment( lib, "opengl32.lib" )
#pragma comment( lib, "glu32.lib" )
#pragma comment( lib, "glut32.lib" )

// Global lenSystem object
static HumanEye eye;
// Global interface variables
static int lvWidth    = LVW;
static int lvHeight   = LVH;
static float lvZoom   = 1.0f;
static float lvTransX = 0.0f;
static float lvTransY = 0.0f;
static bool lvPressed[3] = {false,false,false};
static int lastx = 0, lasty = 0;

// Printf style formatting for rendering text in OpenGL
static void lvText(int x, int y, void *font, int height, 
				   char *format, ...) {
	va_list args; char buffer[256];
	va_start(args, format);
	vsprintf(buffer, format, args);
	va_end(args);
	int w = lvWidth, h = lvHeight;

	glMatrixMode(GL_PROJECTION); 
	glPushMatrix(); glLoadIdentity();
	gluOrtho2D(0, w - 1, 0, h - 1);
	glMatrixMode(GL_MODELVIEW); 
	glPushMatrix(); glLoadIdentity();

	if ((x >= 0) && (y >= 0)) 
		glRasterPos2i(x,y);
	for (char *p = buffer; *p; ++p) 
		glutBitmapCharacter(font, *p);

	glMatrixMode(GL_PROJECTION); glPopMatrix();
	glMatrixMode(GL_MODELVIEW);  glPopMatrix();
}

static void applyTransforms(void) {
	glMatrixMode(GL_MODELVIEW);	 glLoadIdentity();
	glMatrixMode(GL_PROJECTION); glLoadIdentity();
	float orthoX = lvWidth / lvZoom;
	float orthoY = lvHeight / lvZoom;

	gluOrtho2D(orthoX, -orthoX, -orthoY, orthoY);
	glTranslated(lvTransX, lvTransY, 0.0f);
}

static void lvDisplay(void) {
	glClear(GL_COLOR_BUFFER_BIT);
	applyTransforms();
	glColor3f(0.75f, 0.75f, 0.75f);
	glLineWidth(2);
	drawLine(-lvWidth-lvTransX, 0, lvWidth-lvTransX, 0);
	drawLine(0, -lvHeight-lvTransY, 0, lvHeight-lvTransY);
	glLineWidth(1);

	eye.Draw();

	glMatrixMode(GL_PROJECTION); glPushMatrix();
	glLoadIdentity();
	gluOrtho2D(0, lvWidth, lvHeight, 0);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	drawSquare(0.25f, 0.25f, 0.25f, 0, 0, lvWidth, BARH);
	drawSquare(0.25f, 0.25f, 0.25f, 0, lvHeight, lvWidth, lvHeight-BARH);
	glLoadIdentity();

	gluOrtho2D(-lvWidth, lvWidth, lvHeight, -lvHeight);
	glColor3f(0, 0, 0);
	drawText(-lvWidth + 15, -lvTransY*lvZoom - 30, SBFONT, "+z");
	drawText(lvWidth - 50, -lvTransY*lvZoom - 30, SBFONT, "-z");
	drawText(-lvTransX*lvZoom + 10, -lvHeight + 85, SBFONT, "+y");
	drawText(-lvTransX*lvZoom + 10, lvHeight - 75, SBFONT, "-y");

	//drawText(-lvTransX*lvZoom + 105, lvHeight - 170, SBFONT, "eye model");
	//drawText(-lvTransX*lvZoom + 800, lvHeight - 170, SBFONT, "retina plane");

	glPopMatrix();
	glutSwapBuffers();
}

static void lvClick(int button, int state, int x, int y) {
	if (state == GLUT_DOWN)	
		lvPressed[button] = true;
	else if (state == GLUT_UP) 
		lvPressed[button] = false;

	lastx = x; lasty = y; 
	glutPostRedisplay();
}

static void lvMotion(int x, int y) {
	// Here we handle translation
	float norm = 1.0f / lvZoom;
	if (lvPressed[GLUT_LEFT_BUTTON]) {
		lvTransX += float(lastx - x)*2.0f*norm; 
		lvTransY += float(lasty - y)*2.0f*norm;
	}
	//if (lvPressed[GLUT_MIDDLE_BUTTON]) {
	//  // Use this to change camera focus
	//  lsystem.setImagePoint(lsystem.iPoint + Vector(0,lasty-y,lastx-x));
	//  //lsystem.setOPlane(min(lsystem.oPlane, 3000.f) + lastx-x);

	//}
	// Here we handle zooming
	if (lvPressed[GLUT_RIGHT_BUTTON]) {
		float zoom = 15.0f * float(x - lastx) / lvWidth;
		if ((lvZoom + zoom) >= 1.0f)
			lvZoom += zoom;
	}

	lastx = x; lasty = y; 
	glutPostRedisplay();
}

static void lvResize(int w, int h) {
	lvWidth = w; lvHeight = h;
	glViewport(0, 0, lvWidth, lvHeight);
}

static void lvKey(unsigned char key, int x, int y) {
	switch (key) {
	case 27: /* Escape key */ exit(0);
	case ',':
	case 'w':
		eye.MoveUp();
		break;
	case 'o':
	case 's':
		eye.MoveDown();
		break;
	case 'a':
		eye.MoveLeft();
		break;
	case 'e':
	case 'd':
		eye.MoveRight();
		break;
	case 'h':
	case 'j':
		eye.LessDP();
		break;
	case 't':
	case 'k':
		eye.MoreDP();
		break;
	default: break;
	}
	glutPostRedisplay();
}

int main(int argc, char** argv) {

	lvZoom = lvHeight / eye.MaxAper();
	lvZoom = maxv(1, lvZoom);
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowSize(lvWidth, lvHeight);
	glutCreateWindow("Navarro Human Eye Model 2D");
	glutDisplayFunc(lvDisplay);
	glutMotionFunc(lvMotion);
	glutMouseFunc(lvClick);
	glutKeyboardFunc(lvKey);
	glutReshapeFunc(lvResize);
	glutIdleFunc(NULL);

	glMatrixMode(GL_MODELVIEW);  glLoadIdentity();
	glMatrixMode(GL_PROJECTION); glLoadIdentity();
	glClearColor(BGR, BGG, BGB, 1);

	glutMainLoop();
	return (0);
}