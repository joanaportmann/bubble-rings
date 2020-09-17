#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <GL/glut.h>
#include <vector>
struct point
{
        int x;
        int y;
        int z;
};

std::vector<point> points = {
                             {1, 5, 1},
                             {100, 100, 100},
                             {300, 300, 0},
                             {100, 500, 1},
                        //      {0, 3, 0},
                        //      {10, 10, 0},
                        //      {10, 50, 0},
                        //      {0, 0, 0}
                             };


void lineSegment(void)
{
        glClear(GL_COLOR_BUFFER_BIT);
        glColor3f(1.0, 1.0, 0.0);
        glBegin(GL_LINE_LOOP);
        for (int i = 0; i < points.size(); i++)
        {
                glColor3f(1.0, 1.0 * i / points.size(), 0.0);
                printf("%i\n", i);
                printf("%i, %i, %i\n", points[i].x, points[i].y, points[i].z);
                glVertex3d(points[i].x, points[i].y, points[i].z);
        }
        glEnd();
        glFlush();
}


void resize(int width, int height)
{
        glViewport(0, 0, (GLint)width, (GLint)height);
        glMatrixMode(GL_PROJECTION);
        gluOrtho2D(0.0, width, height, 0.0);
        glMatrixMode(GL_MODELVIEW);
}

int main(int argc, char **argv)
{
        glutInit(&argc, argv);
        glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
        glutInitWindowPosition(50, 100);
        glutInitWindowSize(512, 512);
        glutCreateWindow("Ceci n'est pas un bubble ring");
        glutDisplayFunc(lineSegment);
        glutReshapeFunc(resize);
        glutMainLoop();
}
