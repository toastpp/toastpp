/**
 * Procedural Bump Mapping demo.  Uses the example shaders from
 * chapter 11 of the OpenGL Shading Language "orange" book.
 * 16 Jan 2007
 */

#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <GL/glut.h>
#include <GL/glu.h>
#include <GL/glext.h>
#include "extfuncs.h"
#include "shaderutil.h"


static char *FragProgFile = "CH11-bumpmap.frag";
static char *VertProgFile = "CH11-bumpmap.vert";

/* program/shader objects */
static GLuint fragShader;
static GLuint vertShader;
static GLuint program;


static struct uniform_info Uniforms[] = {
   { "LightPosition",    3, GL_FLOAT, { 0.57737, 0.57735, 0.57735, 0.0 }, -1 },
   { "SurfaceColor",     3, GL_FLOAT, { 0.8, 0.8, 0.2, 0 }, -1 },
   { "BumpDensity",      1, GL_FLOAT, { 10.0, 0, 0, 0 }, -1 },
   { "BumpSize",         1, GL_FLOAT, { 0.125, 0, 0, 0 }, -1 },
   { "SpecularFactor",   1, GL_FLOAT, { 0.5, 0, 0, 0 }, -1 },
   END_OF_UNIFORMS
};

static GLint win = 0;

static GLfloat xRot = 20.0f, yRot = 0.0f, zRot = 0.0f;

static GLuint tangentAttrib;

static GLboolean Anim = GL_FALSE;


static void
CheckError(int line)
{
   GLenum err = glGetError();
   if (err) {
      printf("GL Error %s (0x%x) at line %d\n",
             gluErrorString(err), (int) err, line);
   }
}

/*
 * Draw a square, specifying normal and tangent vectors.
 */
static void
Square(GLfloat size)
{
   glNormal3f(0, 0, 1);
   glVertexAttrib3f_func(tangentAttrib, 1, 0, 0);
   glBegin(GL_POLYGON);
   glTexCoord2f(0, 0);  glVertex2f(-size, -size);
   glTexCoord2f(1, 0);  glVertex2f( size, -size);
   glTexCoord2f(1, 1);  glVertex2f( size,  size);
   glTexCoord2f(0, 1);  glVertex2f(-size,  size);
   glEnd();
}


static void
Cube(GLfloat size)
{
   /* +X */
   glPushMatrix();
   glRotatef(90, 0, 1, 0);
   glTranslatef(0, 0, size);
   Square(size);
   glPopMatrix();

   /* -X */
   glPushMatrix();
   glRotatef(-90, 0, 1, 0);
   glTranslatef(0, 0, size);
   Square(size);
   glPopMatrix();

   /* +Y */
   glPushMatrix();
   glRotatef(90, 1, 0, 0);
   glTranslatef(0, 0, size);
   Square(size);
   glPopMatrix();

   /* -Y */
   glPushMatrix();
   glRotatef(-90, 1, 0, 0);
   glTranslatef(0, 0, size);
   Square(size);
   glPopMatrix();


   /* +Z */
   glPushMatrix();
   glTranslatef(0, 0, size);
   Square(size);
   glPopMatrix();

   /* -Z */
   glPushMatrix();
   glRotatef(180, 0, 1, 0);
   glTranslatef(0, 0, size);
   Square(size);
   glPopMatrix();

}


static void
Idle(void)
{
   GLint t = glutGet(GLUT_ELAPSED_TIME);
   yRot  = t * 0.05;
   glutPostRedisplay();
}


static void
Redisplay(void)
{
   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
   
   glPushMatrix();
   glRotatef(xRot, 1.0f, 0.0f, 0.0f);
   glRotatef(yRot, 0.0f, 1.0f, 0.0f);
   glRotatef(zRot, 0.0f, 0.0f, 1.0f);

   Cube(1.5);

   glPopMatrix();

   CheckError(__LINE__);

   glutSwapBuffers();
}


static void
Reshape(int width, int height)
{
   glViewport(0, 0, width, height);
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   glFrustum(-1.0, 1.0, -1.0, 1.0, 5.0, 25.0);
   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();
   glTranslatef(0.0f, 0.0f, -15.0f);
}


static void
CleanUp(void)
{
   glDeleteShader_func(fragShader);
   glDeleteShader_func(vertShader);
   glDeleteProgram_func(program);
   glutDestroyWindow(win);
}


static void
Key(unsigned char key, int x, int y)
{
   const GLfloat step = 2.0;
  (void) x;
  (void) y;

   switch(key) {
   case 'a':
      Anim = !Anim;
      glutIdleFunc(Anim ? Idle : NULL);
      break;
   case 'z':
      zRot += step;
      break;
   case 'Z':
      zRot -= step;
      break;
   case 27:
      CleanUp();
      exit(0);
      break;
   }
   glutPostRedisplay();
}


static void
SpecialKey(int key, int x, int y)
{
   const GLfloat step = 2.0;

  (void) x;
  (void) y;

   switch(key) {
   case GLUT_KEY_UP:
      xRot += step;
      break;
   case GLUT_KEY_DOWN:
      xRot -= step;
      break;
   case GLUT_KEY_LEFT:
      yRot -= step;
      break;
   case GLUT_KEY_RIGHT:
      yRot += step;
      break;
   }
   glutPostRedisplay();
}


static void
Init(void)
{
   if (!ShadersSupported())
      exit(1);

   GetExtensionFuncs();

   vertShader = CompileShaderFile(GL_VERTEX_SHADER, VertProgFile);
   fragShader = CompileShaderFile(GL_FRAGMENT_SHADER, FragProgFile);
   program = LinkShaders(vertShader, fragShader);

   glUseProgram_func(program);

   assert(glIsProgram_func(program));
   assert(glIsShader_func(fragShader));
   assert(glIsShader_func(vertShader));

   assert(glGetError() == 0);

   CheckError(__LINE__);

   InitUniforms(program, Uniforms);

   CheckError(__LINE__);

   tangentAttrib = glGetAttribLocation_func(program, "Tangent");
   printf("Tangent Attrib: %d\n", tangentAttrib);

   assert(tangentAttrib >= 0);

   CheckError(__LINE__);

   glClearColor(0.4f, 0.4f, 0.8f, 0.0f);

   glEnable(GL_DEPTH_TEST);

   glColor3f(1, 0, 0);
}


static void
ParseOptions(int argc, char *argv[])
{
   int i;
   for (i = 1; i < argc; i++) {
      if (strcmp(argv[i], "-fs") == 0) {
         FragProgFile = argv[i+1];
      }
      else if (strcmp(argv[i], "-vs") == 0) {
         VertProgFile = argv[i+1];
      }
   }
}


int
main(int argc, char *argv[])
{
   glutInit(&argc, argv);
   glutInitWindowPosition( 0, 0);
   glutInitWindowSize(400, 400);
   glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
   win = glutCreateWindow(argv[0]);
   glutReshapeFunc(Reshape);
   glutKeyboardFunc(Key);
   glutSpecialFunc(SpecialKey);
   glutDisplayFunc(Redisplay);
   ParseOptions(argc, argv);
   Init();
   glutMainLoop();
   return 0;
}

