#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>

#include "extra.h"
#include "camera.h"

#include "TimeStepper.h"
#include "simpleSystem.h"
#include "pendulumSystem.h"
#include "ClothSystem.h"
#include "GridTestSystem.h"
#include "sphfluidsystem.h"

using namespace std;

// Globals here.
int numParticles = 103;

// For movement of the cloth along the z axis
float clothSystemZBoundaryPos = 6.0;
float clothSystemZBoundaryNeg = -6.0;
int pendulumSize = 5;
bool clothMoving;

namespace
{
    ParticleSystem *system;
    TimeStepper *timeStepper;
    float stepSize;

    enum IntegratorType
    {
        FORWARD_EULER,
        TRAPEZOIDAL,
        RUNGE_KUTTA
    };

    IntegratorType integratorType;

  void initSystem(int argc, char * argv[])
  {
    // Seed the random number generator with the current time
    srand( time( NULL ) );

    system = new SPHFluidSystem(numParticles);

	timeStepper = new RK4();

	integratorType = IntegratorType::RUNGE_KUTTA;

    stepSize = 0.02;
  }

  void stepSystem()
  {
      if(timeStepper!=0)
      {
          timeStepper->takeStep(system,stepSize);
      }

      // Collision detection
      /*
      vector<Vector3f> state = system->getState();
      for (int i = 0; i < numParticles; ++i)
      {
          Vector3f pos = state[2 * i];
          Vector3f vel = state[2 * i + 1];
      }
      */
  }


  // Draw the current particle positions
  void drawSystem()
  {
    // Base material colors (they don't change)
    GLfloat particleColor[] = {0.4f, 0.7f, 1.0f, 1.0f};
    GLfloat tankColor[] = {0.0f, 0.0f, 1.0f, 0.3f};
    
    glTranslatef(-1.0f, -1.0f, 0.0f); // Origin at (-1.0, -1.0, 0.0). Everything drawn relative to this point.

    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, particleColor);    

    system->draw();

    glEnable(GL_BLEND); //Enable alpha blending
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); //Set the blend function
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, tankColor);

    // Draw the tank

    // Back
    glPushMatrix();
    glTranslatef(1.0f, 1.0f, 0.0f);
    glScaled(2.0f, 2.0f, 0.1f);
    glutSolidCube(1);
    glPopMatrix();

    // Front
    glPushMatrix();
    glTranslatef(1.0f, 1.0f, 2.0f);
    glScaled(2.0f, 2.0f, 0.1f);
    glutSolidCube(1);
    glPopMatrix();

    // Bottom
    glPushMatrix();
    glTranslatef(1.0f, 0.0f, 1.0f);
    glScaled(2.0f, 0.1f, 2.0f);
    glutSolidCube(1);
    glPopMatrix();

    // Left
    glPushMatrix();
    glTranslatef(0.0f, 1.0f, 1.0f);
    glScaled(0.1f, 2.0f, 2.0f);
    glutSolidCube(1);
    glPopMatrix();

    // Right
    glPushMatrix();
    glTranslatef(2.0f, 1.0f, 1.0f);
    glScaled(0.1f, 2.0f, 2.0f);
    glutSolidCube(1);
    glPopMatrix();
}
        

    //-------------------------------------------------------------------
    
        
    // This is the camera
    Camera camera;

    // These are state variables for the UI
    bool g_mousePressed = false;

    // Declarations of functions whose implementations occur later.
    void arcballRotation(int endX, int endY);
    void keyboardFunc( unsigned char key, int x, int y);
    void specialFunc( int key, int x, int y );
    void mouseFunc(int button, int state, int x, int y);
    void motionFunc(int x, int y);
    void reshapeFunc(int w, int h);
    void drawScene(void);
    void initRendering();

    // This function is called whenever a "Normal" key press is
    // received.
    void keyboardFunc( unsigned char key, int x, int y )
    {
        switch ( key )
        {
            case 27: // Escape key
            {
                exit(0);
                break;
            }

            case ' ':
            {
                Matrix4f eye = Matrix4f::identity();
                camera.SetRotation( eye );
                camera.SetCenter( Vector3f::ZERO );
                break;
            }

            case 'r':
            {
                system->reinitializeSystem();
                break;
            }

            case '8':
            {
                if (integratorType != FORWARD_EULER)
                {
                    integratorType = FORWARD_EULER;
                    timeStepper = new ForwardEuler();
                    system->reinitializeSystem();
                }

                break;
            }

            case '9':
            {
                if (integratorType != TRAPEZOIDAL)
                {
                    integratorType = TRAPEZOIDAL;
                    timeStepper = new Trapezoidal();
                    system->reinitializeSystem();
                }

                break;
            }

            case '0':
            {
                if (integratorType != RUNGE_KUTTA)
                {
                    integratorType = RUNGE_KUTTA;
                    timeStepper = new RK4();
                    system->reinitializeSystem();
                }
                break;
            }


            default:
                cout << "Unhandled key press " << key << "." << endl;
                break;
            }

        glutPostRedisplay();
    }

    // This function is called whenever a "Special" key press is
    // received.  Right now, it's handling the arrow keys.
    void specialFunc( int key, int x, int y )
    {
        switch ( key )
        {

        }
        //glutPostRedisplay();
    }

    //  Called when mouse button is pressed.
    void mouseFunc(int button, int state, int x, int y)
    {
        if (state == GLUT_DOWN)
        {
            g_mousePressed = true;

            //No rotating the camera for now
            switch (button)
            {
            case GLUT_LEFT_BUTTON:
                camera.MouseClick(Camera::LEFT, x, y);
                break;
            case GLUT_MIDDLE_BUTTON:
                camera.MouseClick(Camera::MIDDLE, x, y);
                break;
            case GLUT_RIGHT_BUTTON:
                camera.MouseClick(Camera::RIGHT, x,y);
                break;
            default:
                break;
            }

        }
        else
        {
            camera.MouseRelease(x,y);
            g_mousePressed = false;
        }
        glutPostRedisplay();
    }

    // Called when mouse is moved while button pressed.
    void motionFunc(int x, int y)
    {
        camera.MouseDrag(x,y);        
    
        glutPostRedisplay();
    }

    // Called when the window is resized
    // w, h - width and height of the window in pixels.
    void reshapeFunc(int w, int h)
    {
        camera.SetDimensions(w,h);

        camera.SetViewport(0,0,w,h);
        camera.ApplyViewport();

        // Set up a perspective view, with square aspect ratio
        glMatrixMode(GL_PROJECTION);

        camera.SetPerspective(50);
        glLoadMatrixf( camera.projectionMatrix() );
    }

    // Initialize OpenGL's rendering modes
    void initRendering()
    {
        glEnable(GL_DEPTH_TEST);   // Depth testing must be turned on
        glEnable(GL_LIGHTING);     // Enable lighting calculations
        glEnable(GL_LIGHT0);       // Turn on light #0.

        glEnable(GL_NORMALIZE);

        // Setup polygon drawing
        glShadeModel(GL_SMOOTH);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        
        glEnable(GL_CULL_FACE);
        glCullFace(GL_BACK);

        // Clear to black
        glClearColor(0,0,0,1);
    }

    // This function is responsible for displaying the object.
    void drawScene(void)
    {
        // Clear the rendering window
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glMatrixMode( GL_MODELVIEW );  
        glLoadIdentity();              

        // Light color (RGBA)
        GLfloat Lt0diff[] = {1.0,1.0,1.0,1.0};
        GLfloat Lt0pos[] = {3.0,3.0,5.0,1.0};
        glLightfv(GL_LIGHT0, GL_DIFFUSE, Lt0diff);
        glLightfv(GL_LIGHT0, GL_POSITION, Lt0pos);

        glLoadMatrixf( camera.viewMatrix() );

        drawSystem();

        // This draws the coordinate axes when you're rotating, to
        // keep yourself oriented.
        if( g_mousePressed )
        {
            glPushMatrix();
            Vector3f eye = camera.GetCenter();
            glTranslatef( eye[0], eye[1], eye[2] );

            // Save current state of OpenGL
            glPushAttrib(GL_ALL_ATTRIB_BITS);

            // This is to draw the axes when the mouse button is down
            glDisable(GL_LIGHTING);
            glLineWidth(3);
            glPushMatrix();
            glScaled(5.0,5.0,5.0);
            glBegin(GL_LINES);
            glColor4f(1,0.5,0.5,1); glVertex3f(0,0,0); glVertex3f(1,0,0);
            glColor4f(0.5,1,0.5,1); glVertex3f(0,0,0); glVertex3f(0,1,0);
            glColor4f(0.5,0.5,1,1); glVertex3f(0,0,0); glVertex3f(0,0,1);

            glColor4f(0.5,0.5,0.5,1);
            glVertex3f(0,0,0); glVertex3f(-1,0,0);
            glVertex3f(0,0,0); glVertex3f(0,-1,0);
            glVertex3f(0,0,0); glVertex3f(0,0,-1);

            glEnd();
            glPopMatrix();

            glPopAttrib();
            glPopMatrix();
        }
                 
        // Dump the image to the screen.
        glutSwapBuffers();
    }

    void timerFunc(int t)
    {
        stepSystem();

        glutPostRedisplay();

        glutTimerFunc(t, &timerFunc, t);
    }
}

// Main routine.
// Set up OpenGL, define the callbacks and start the main loop
int main( int argc, char* argv[] )
{
    glutInit( &argc, argv );

    // We're going to animate it, so double buffer 
    glutInitDisplayMode( GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH );

    // Initial parameters for window position and size
    glutInitWindowPosition( 60, 60 );
    glutInitWindowSize( 1000, 1000 );
    
    camera.SetDimensions( 1000, 1000 );

    camera.SetDistance( 10 );
    camera.SetCenter(Vector3f::ZERO);
    
    glutCreateWindow("SPH Fluid Simulation");

    // Initialize OpenGL parameters.
    initRendering();

    // Setup particle system
    initSystem(argc,argv);

    // Set up callback functions for key presses
    glutKeyboardFunc(keyboardFunc); // Handles "normal" ascii symbols
    glutSpecialFunc(specialFunc);   // Handles "special" keyboard keys

    // Set up callback functions for mouse
    glutMouseFunc(mouseFunc);
    glutMotionFunc(motionFunc);

    // Set up the callback function for resizing windows
    glutReshapeFunc( reshapeFunc );

    // Call this whenever window needs redrawing
    glutDisplayFunc( drawScene );

    // Trigger timerFunc every 20 msec
    glutTimerFunc(20, timerFunc, 20);
        
    // Start the main loop.  glutMainLoop never returns.
    glutMainLoop();

    return 0;	// This line is never reached.
}
