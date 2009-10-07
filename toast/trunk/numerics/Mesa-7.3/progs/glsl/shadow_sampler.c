/**
 * Test shadow2D() function.
 * Brian Paul
 * 11 April 2007
 */

#define GL_GLEXT_PROTOTYPES
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <GL/gl.h>
#include <GL/glut.h>
#include <GL/glext.h>
#include "extfuncs.h"


static char *FragProgFile = NULL;
static char *VertProgFile = NULL;

static GLuint fragShader;
static GLuint vertShader;
static GLuint program;

static GLint uTexture2D;
static GLint uTextureRect;

static GLint win = 0;

static GLenum Filter = GL_LINEAR;

static void
CheckError(int line)
{
   GLenum err = glGetError();
   if (err) {
      printf("GL Error %s (0x%x) at line %d\n",
             gluErrorString(err), (int) err, line);
   }
}


static void
PrintString(const char *s)
{
   while (*s) {
      glutBitmapCharacter(GLUT_BITMAP_8_BY_13, (int) *s);
      s++;
   }
}


static void
Redisplay(void)
{
   GLfloat q = 2.0;
   CheckError(__LINE__);
   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

   glPushMatrix();

   CheckError(__LINE__);
   glUseProgram_func(program);
   CheckError(__LINE__);

   glBegin(GL_POLYGON);
#define RECT 1
#if RECT
   glTexCoord4f( 0,  0, 0, q);  glVertex2f(-1, -1);
   glTexCoord4f(16,  0, 1, q);  glVertex2f( 1, -1);
   glTexCoord4f(16,  1, 1, q);  glVertex2f( 1,  1);
   glTexCoord4f( 0, 10, 0, q);  glVertex2f(-1,  1);
#else
   glTexCoord3f(0, 0, 0);  glVertex2f(-1, -1);
   glTexCoord3f(1, 0, 1);  glVertex2f( 1, -1);
   glTexCoord3f(1, 1, 1);  glVertex2f( 1,  1);
   glTexCoord3f(0, 1, 0);  glVertex2f(-1,  1);
#endif
   glEnd();

   glPopMatrix();

   glUseProgram_func(0);
   glWindowPos2iARB(80, 20);
   PrintString("white   black   white   black");

   glFlush();
   /*
   glutSwapBuffers();
   */
   CheckError(__LINE__);
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
   glTranslatef(0.0f, 0.0f, -8.0f);
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
   case 27:
      CleanUp();
      exit(0);
      break;
   }
   glutPostRedisplay();
}


static void
MakeTexture(void)
{
#define SZ 16
   GLfloat image[SZ][SZ];
   GLuint i, j;

   for (i = 0; i < SZ; i++) {
      for (j = 0; j < SZ; j++) {
         if (j < (SZ / 2)) {
            image[i][j] = 0.25;
         }
         else {
            image[i][j] = 0.75;
         }
      }
   }

   glActiveTexture(GL_TEXTURE0); /* unit 0 */
   glBindTexture(GL_TEXTURE_2D, 42);
   glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, SZ, SZ, 0,
                GL_DEPTH_COMPONENT, GL_FLOAT, image);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 1);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, Filter);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, Filter);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE_ARB,
                   GL_COMPARE_R_TO_TEXTURE_ARB);
   CheckError(__LINE__);

   glActiveTexture(GL_TEXTURE1); /* unit 1 */
   glBindTexture(GL_TEXTURE_RECTANGLE_ARB, 43);
   glTexImage2D(GL_TEXTURE_RECTANGLE_ARB, 0, GL_DEPTH_COMPONENT,
                SZ, 10, 0,/*16x10*/
                GL_DEPTH_COMPONENT, GL_FLOAT, image);
   glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MIN_FILTER, Filter);
   glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MAG_FILTER, Filter);
   glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_COMPARE_MODE_ARB,
                   GL_COMPARE_R_TO_TEXTURE_ARB);
   CheckError(__LINE__);
#undef SZ
}


static void
LoadAndCompileShader(GLuint shader, const char *text)
{
   GLint stat;
   glShaderSource_func(shader, 1, (const GLchar **) &text, NULL);
   glCompileShader_func(shader);
   glGetShaderiv_func(shader, GL_COMPILE_STATUS, &stat);
   if (!stat) {
      GLchar log[1000];
      GLsizei len;
      glGetShaderInfoLog_func(shader, 1000, &len, log);
      fprintf(stderr, "fslight: problem compiling shader:\n%s\n", log);
      exit(1);
   }
}


/**
 * Read a shader from a file.
 */
static void
ReadShader(GLuint shader, const char *filename)
{
   const int max = 100*1000;
   int n;
   char *buffer = (char*) malloc(max);
   FILE *f = fopen(filename, "r");
   if (!f) {
      fprintf(stderr, "fslight: Unable to open shader file %s\n", filename);
      exit(1);
   }

   n = fread(buffer, 1, max, f);
   printf("fslight: read %d bytes from shader file %s\n", n, filename);
   if (n > 0) {
      buffer[n] = 0;
      LoadAndCompileShader(shader, buffer);
   }

   fclose(f);
   free(buffer);
}


static void
CheckLink(GLuint prog)
{
   GLint stat;
   glGetProgramiv_func(prog, GL_LINK_STATUS, &stat);
   if (!stat) {
      GLchar log[1000];
      GLsizei len;
      glGetProgramInfoLog_func(prog, 1000, &len, log);
      fprintf(stderr, "Linker error:\n%s\n", log);
   }
}


static void
Init(void)
{
   static const char *fragShaderText =
      "uniform sampler2DShadow shadowTex2D; \n"
      "uniform sampler2DRectShadow shadowTexRect; \n"
      "void main() {\n"
#if RECT
      "   gl_FragColor = shadow2DRectProj(shadowTexRect, gl_TexCoord[0]); \n"
#else
      "   gl_FragColor = shadow2D(shadowTex2D, gl_TexCoord[0].xyz); \n"
#endif
      "}\n";
   static const char *vertShaderText =
      "void main() {\n"
      "   gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;\n"
      "   gl_TexCoord[0] = gl_MultiTexCoord0; \n"
      "}\n";
   const char *version;

   version = (const char *) glGetString(GL_VERSION);
   if (version[0] != '2' || version[1] != '.') {
      printf("This program requires OpenGL 2.x, found %s\n", version);
      exit(1);
   }
   printf("GL_RENDERER = %s\n",(const char *) glGetString(GL_RENDERER));

   GetExtensionFuncs();

   fragShader = glCreateShader_func(GL_FRAGMENT_SHADER);
   if (FragProgFile)
      ReadShader(fragShader, FragProgFile);
   else
      LoadAndCompileShader(fragShader, fragShaderText);

   vertShader = glCreateShader_func(GL_VERTEX_SHADER);
   if (VertProgFile)
      ReadShader(vertShader, VertProgFile);
   else
      LoadAndCompileShader(vertShader, vertShaderText);

   program = glCreateProgram_func();
   glAttachShader_func(program, fragShader);
   glAttachShader_func(program, vertShader);
   glLinkProgram_func(program);
   CheckLink(program);
   glUseProgram_func(program);

   uTexture2D = glGetUniformLocation_func(program, "shadowTex2D");
   uTextureRect = glGetUniformLocation_func(program, "shadowTexRect");
   printf("uTexture2D %d  uTextureRect %d\n", uTexture2D, uTextureRect);
   if (uTexture2D >= 0) {
      glUniform1i_func(uTexture2D, 0);  /* use texture unit 0 */
   }
   if (uTextureRect >= 0) {
      glUniform1i_func(uTextureRect, 1);  /* use texture unit 0 */
   }
   CheckError(__LINE__);

   glClearColor(0.3f, 0.3f, 0.3f, 0.0f);
   glColor3f(1, 1, 1);

   MakeTexture();
   CheckError(__LINE__);
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
   glutInitWindowSize(400, 300);
   glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE | GLUT_DEPTH);
   win = glutCreateWindow(argv[0]);
   glutReshapeFunc(Reshape);
   glutKeyboardFunc(Key);
   glutDisplayFunc(Redisplay);
   ParseOptions(argc, argv);
   Init();
   glutMainLoop();
   return 0;
}


