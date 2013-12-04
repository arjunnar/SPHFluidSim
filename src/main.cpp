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

// Globals
float boxSizeX;
float boxSizeY;
float boxSizeZ;
float clothScale = 4.0;
int numIters = 0;
int numItersPerEmit = 1;

namespace
{
    SPHFluidSystem *system;
    TimeStepper *timeStepper;
    float stepSize;

    enum IntegratorType
    {
        FORWARD_EULER,
        LEAP_FROG
    };

    FluidSystemType fluidSystemType;

    IntegratorType integratorType;

  void initSystem(int argc, char * argv[])
  {
    // Seed the random number generator with the current time
    srand( time( NULL ) );

    fluidSystemType = FluidSystemType::System2DEmitter;

    switch (fluidSystemType)
    {
         case TwoDensitySystem2D:
            boxSizeX = 1.0;
            boxSizeY = 1.5;
            boxSizeZ = 1.0;
            break;

         case SystemSimple2D:
         case SystemLarge2D:
            boxSizeX = 0.8;
            boxSizeY = 0.9;
            boxSizeZ = 0.8;
            break;

         case SystemSimple3D:
            boxSizeX = 0.4;
            boxSizeY = 0.4;
            boxSizeZ = 0.4;
            break;

         case SystemLarge3D:
            boxSizeX = 0.45;
            boxSizeY = 1.1;
            boxSizeZ = 0.45;
            break;

         case System2DEmitter:
            boxSizeX = 1.1;
            boxSizeY = 1.1;
            boxSizeZ = 0.3;
            break;
glColor3f(0.0, 0.0, 1.0);
    }

    system = new SPHFluidSystem(boxSizeX, boxSizeY, boxSizeZ, fluidSystemType);

    timeStepper = new LeapFrog();

    integratorType = IntegratorType::LEAP_FROG;

    stepSize = 0.005;
  }

  bool fixCoord(float &coord, float &vel, float boxSize)
  {
      float collisionEpsilon = 0.02;
	  bool fixed = false;

      if (coord < 0.0f + collisionEpsilon)
	  {
		  coord = collisionEpsilon;
          vel = -0.5 * vel;
		  fixed = true;
	  }

      else if (coord > boxSize - collisionEpsilon)
	  {
          coord = boxSize - collisionEpsilon;
          vel = -0.5 * vel;
		  fixed = true;
	  }

	  return fixed;
   }


  void stepSystem()
  {
     if(timeStepper!=0)
     {
         if (numIters % numItersPerEmit == 0)
         {
            system->emitParticle();
            numIters = 0;
         }

         timeStepper->takeStep(system,stepSize);
         numIters++;
     }

      //system->advanceState();

      vector<Vector3f> state = system->getState();
      for (vector<Vector3f>::iterator iter = state.begin(); iter != state.end(); iter += 2)
      {
          Vector3f pos = *iter;
          Vector3f vel = *(iter + 1);

          float x = pos.x();
          float y = pos.y();
          float z = pos.z();

          float velX = vel.x();
          float velY = vel.y();
          float velZ = vel.z();

          bool fixedX = fixCoord(x, velX, boxSizeX);
          bool fixedY = fixCoord(y, velY, boxSizeY);
          bool fixedZ = fixCoord(z, velZ, boxSizeZ);

          /*
          if (fixedZ)
          {
              velY = 1.0;
              velX = 1.0;
          }
          */

          if (fixedX || fixedY || fixedZ)
          {
        	  pos = Vector3f(x, y, z);
        	  vel = Vector3f(velX, velY, velZ);
        	  *iter = pos;
        	  *(iter + 1) = vel;
          }
      }

      system->setState(state);
  }

  // Draw the current particle positions
  void drawSystem()
  {
    // Base material colors (they don't change)
    GLfloat tankColor[] = {160/255.0, 160/255.0, 160/255.0, 0.1f};
    GLfloat floorColor[] = {205.0/255.0, 133.0/255.0, 63/255.0, 0.80};

    //glPushMatrix();

    //glTranslatef(-boxSizeX/2, 0, boxSizeZ/2);

    system->draw();

    glEnable(GL_BLEND); //Enable alpha blending
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); //Set the blend function
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, tankColor);

    // Draw the tank

    // Back
    glPushMatrix();
    glTranslatef(boxSizeX/2.0f, boxSizeY/2.0f, 0.0f);
    glScaled(boxSizeX, boxSizeY, 0.05f);
    glutSolidCube(1);
    glPopMatrix();

    // Front
    glPushMatrix();
    glTranslatef(boxSizeX/2.0f, boxSizeY/2.0f, boxSizeZ);
    glScaled(boxSizeX, boxSizeY, 0.02f);
    glutSolidCube(1);
    glPopMatrix();

    // Bottom
    glPushMatrix();
    glTranslatef(boxSizeX/2.0f, 0.0f, boxSizeZ/2.0f);
    glScaled(boxSizeX, 0.02f, boxSizeZ);
    glutSolidCube(1);
    glPopMatrix();

    // Left
    glPushMatrix();
    glTranslatef(0.0f, boxSizeY/2.0f, boxSizeZ/2.0f);
    glScaled(0.02f, boxSizeY, boxSizeZ);
    glutSolidCube(1);
    glPopMatrix();

    // Right
    glPushMatrix();

    glTranslatef(boxSizeX,  boxSizeY/2.0f,  boxSizeZ/2.0f);
    glScaled(0.02f, boxSizeY, boxSizeZ);
    glutSolidCube(1);
    glPopMatrix();

    // Draw the floor
    glPushMatrix();
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, floorColor);
    glTranslatef(boxSizeX/2.0f, -0.1f, boxSizeZ/2.0f);
    glScaled(clothScale, 0.02f, clothScale);
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

        glPushMatrix();
        glTranslatef(-boxSizeX/2, -boxSizeY/2, boxSizeZ/2);

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

        glPopMatrix();

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

    camera.SetDistance( 3.0 );
    camera.SetCenter(Vector3f::ZERO);
    
    glutCreateWindow("Smoothed Particle Hydrodynamics Fluid Simulation");

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
