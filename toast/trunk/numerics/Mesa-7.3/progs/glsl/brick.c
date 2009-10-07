/**
 * "Brick" shader demo.  Uses the example shaders from chapter 6 of
 * the OpenGL Shading Language "orange" book.
 * 10 Jan 2007
 */

#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <GL/gl.h>
#include <GL/glut.h>
#include <GL/glext.h>
#include "extfuncs.h"
#include "shaderutil.h"


static char *FragProgFile = "CH06-brick.frag";
static char *VertProgFile = "CH06-brick.vert";

/* program/shader objects */
static GLuint fragShader;
static GLuint vertShader;
static GLuint program;

static struct uniform_info Uniforms[] = {
   /* vert */
   { "LightPosition",     3, GL_FLOAT, { 0.1, 0.1, 9.0, 0}, -1 },
   /* frag */
   { "BrickColor",        3, GL_FLOAT, { 0.8, 0.2, 0.2, 0 }, -1 },
   { "MortarColor",       3, GL_FLOAT, { 0.6, 0.6, 0.6, 0 }, -1 },
   { "BrickSize",         2, GL_FLOAT, { 1.0, 0.3, 0, 0 }, -1 },
   { "BrickPct",          2, GL_FLOAT, { 0.9, 0.8, 0, 0 }, -1 },
   END_OF_UNIFORMS
};

static GLint win = 0;


static GLfloat xRot = 0.0f, yRot = 0.0f, zRot = 0.0f;




static void
Redisplay(void)
{
   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
   
   glPushMatrix();
   glRotatef(xRot, 1.0f, 0.0f, 0.0f);
   glRotatef(yRot, 0.0f, 1.0f, 0.0f);
   glRotatef(zRot, 0.0f, 0.0f, 1.0f);

   glBegin(GL_POLYGON);
   glTexCoord2f(0, 0);   glVertex2f(-2, -2);
   glTexCoord2f(1, 0);   glVertex2f( 2, -2);
   glTexCoord2f(1, 1);   glVertex2f( 2,  2);
   glTexCoord2f(0, 1);   glVertex2f(-2,  2);
   glEnd();

   glPopMatrix();

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
  (void) x;
  (void) y;

   switch(key) {
   case 'z':
      zRot -= 1.0;
      break;
   case 'Z':
      zRot += 1.0;
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
   const GLfloat step = 3.0f;

  (void) x;
  (void) y;

   switch(key) {
   case GLUT_KEY_UP:
      xRot -= step;
      break;
   case GLUT_KEY_DOWN:
      xRot += step;
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

   InitUniforms(program, Uniforms);

   assert(glGetError() == 0);

   glClearColor(0.4f, 0.4f, 0.8f, 0.0f);

   printf("GL_RENDERER = %s\n",(const char *) glGetString(GL_RENDERER));

   assert(glIsProgram_func(program));
   assert(glIsShader_func(fragShader));
   assert(glIsShader_func(vertShader));

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

