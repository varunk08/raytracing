
//-------------------------------------------------------------------------------
///
/// \file       viewport.cpp 
/// \author     Cem Yuksel (www.cemyuksel.com)
/// \version    7.0
/// \date       October 7, 2013
///
/// \brief Example source for CS 6620 - University of Utah.
///
//-------------------------------------------------------------------------------
 
#include "scene.h"
#include <GL/glut.h>
#include <time.h>
 
//-------------------------------------------------------------------------------
//void Sphere::ViewportDisplay() const
//{
//  static GLUquadric *q = NULL;
//  if ( q == NULL ) {
//      q = gluNewQuadric();
//      gluQuadricTexture(q,true);
//  }
//  gluSphere(q,1,50,50);
//}
//void Plane::ViewportDisplay() const
//{
//  const int resolution = 32;
//  float xyInc = 2.0f / resolution;
//  float uvInc = 1.0f / resolution;
//  glPushMatrix();
//  glNormal3f(0,0,1);
//  glBegin(GL_QUADS);
//  float y1=-1, y2=xyInc-1, v1=0, v2=uvInc;
//  for ( int y=0; y<resolution; y++ ) {
//      float x1=-1, x2=xyInc-1, u1=0, u2=uvInc;
//      for ( int x=0; x<resolution; x++ ) {
//          glTexCoord2f(u1, v1);
//          glVertex3f ( x1, y1, 0 );
//          glTexCoord2f(u2, v1);
//          glVertex3f ( x2, y1, 0 );
//          glTexCoord2f(u2, v2);
//          glVertex3f ( x2, y2, 0 );
//          glTexCoord2f(u1, v2);
//          glVertex3f ( x1, y2, 0 );
//          x1=x2; x2+=xyInc; u1=u2; u2+=uvInc;
//      }
//      y1=y2; y2+=xyInc; v1=v2; v2+=uvInc;
//  }
//  glEnd();
//  glPopMatrix();
//}
//void TriObj::ViewportDisplay() const
//{
//  glBegin(GL_TRIANGLES);
//  for ( unsigned int i=0; i<NF(); i++ ) {
//      for ( int j=0; j<3; j++ ) {
//          if ( HasTextureVertices() ) glTexCoord3fv( &VT( FT(i).v[j] ).x );
//          if ( HasNormals() ) glNormal3fv( &VN( FN(i).v[j] ).x );
//          glVertex3fv( &V( F(i).v[j] ).x );
//      }
//  }
//  glEnd();
//}
//void MtlBlinn::SetViewportMaterial() const
//{
//  cyColorA c;
//  c = diffuse.GetColor();
//  glMaterialfv( GL_FRONT, GL_AMBIENT_AND_DIFFUSE, &c.r );
//  c = specular.GetColor();
//  glMaterialfv( GL_FRONT, GL_SPECULAR, &c.r );
//  glMaterialf( GL_FRONT, GL_SHININESS, glossiness*1.5f );
//  const TextureMap *dm = diffuse.GetTexture();
//  if ( dm && dm->SetViewportTexture() ) {
//      glEnable( GL_TEXTURE_2D );
//      glMatrixMode( GL_TEXTURE );
//      Matrix3 tm = dm->GetInverseTransform();
//      Point3 p = tm * dm->GetPosition();
//      float m[16] = { tm[0],tm[1],tm[2],0, tm[3],tm[4],tm[5],0, tm[6],tm[7],tm[8],0, -p.x,-p.y,-p.z,1 };
//      glLoadMatrixf( m );
//      glMatrixMode( GL_MODELVIEW );
//  } else {
//      glDisable( GL_TEXTURE_2D );
//  }
//}
//void GenLight::SetViewportParam(int lightID, ColorA ambient, ColorA intensity, Point4 pos ) const
//{
//  glEnable ( GL_LIGHT0 + lightID );
//  glLightfv( GL_LIGHT0 + lightID, GL_AMBIENT,  &ambient.r );
//  glLightfv( GL_LIGHT0 + lightID, GL_DIFFUSE,  &intensity.r );
//  glLightfv( GL_LIGHT0 + lightID, GL_SPECULAR, &intensity.r );
//  glLightfv( GL_LIGHT0 + lightID, GL_POSITION, &pos.x );
//}
//-------------------------------------------------------------------------------
 
void BeginRender(); // Called to start rendering (renderer must run in a separate thread)
void StopRender();  // Called to end rendering (if it is not already finished)
 
extern Node rootNode;
extern Camera camera;
extern RenderImage renderImage;
extern LightList lights;
extern TexturedColor background;
 
//-------------------------------------------------------------------------------
 
enum Mode {
    MODE_READY,         // Ready to render
    MODE_RENDERING,     // Rendering the image
    MODE_RENDER_DONE    // Rendering is finished
};
 
enum ViewMode
{
    VIEWMODE_OPENGL,
    VIEWMODE_IMAGE,
    VIEWMODE_Z,
};
 
static Mode     mode        = MODE_READY;       // Rendering mode
static ViewMode viewMode    = VIEWMODE_OPENGL;  // Display mode
static int      startTime;                      // Start time of rendering
 
//-------------------------------------------------------------------------------
 
void GlutDisplay();
void GlutReshape(int w, int h);
void GlutIdle();
void GlutKeyboard(unsigned char key, int x, int y);
void GlutMouse(int button, int state, int x, int y);
void GlutMotion(int x, int y);
 
//-------------------------------------------------------------------------------
 
void ShowViewport()
{
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH );
    if (glutGet(GLUT_SCREEN_WIDTH) > 0 && glutGet(GLUT_SCREEN_HEIGHT) > 0){
        glutInitWindowPosition( (glutGet(GLUT_SCREEN_WIDTH) - camera.imgWidth)/2, (glutGet(GLUT_SCREEN_HEIGHT) - camera.imgHeight)/2 );
    }
    else glutInitWindowPosition( 50, 50 );
    glutInitWindowSize(camera.imgWidth, camera.imgHeight);
 
    glutCreateWindow("Ray Tracer - CS 6620");
    glutDisplayFunc(GlutDisplay);
    glutReshapeFunc(GlutReshape);
    glutIdleFunc(GlutIdle);
    glutKeyboardFunc(GlutKeyboard);
    glutMouseFunc(GlutMouse);
    glutMotionFunc(GlutMotion);
 
    Color bg = background.GetColor();
    glClearColor(bg.r,bg.g,bg.b,0);
 
    glPointSize(3.0);
    glEnable( GL_CULL_FACE );
 
    float zero[] = {0,0,0,0};
    glLightModelfv( GL_LIGHT_MODEL_AMBIENT, zero );
 
    glEnable(GL_NORMALIZE);
 
    glLineWidth(2);
 
    glutMainLoop();
}
 
//-------------------------------------------------------------------------------
 
void GlutReshape(int w, int h)
{
    if( w != camera.imgWidth || h != camera.imgHeight ) {
        glutReshapeWindow( camera.imgWidth, camera.imgHeight);
    } else {
        glViewport( 0, 0, w, h );
 
        glMatrixMode( GL_PROJECTION );
        glLoadIdentity();
        float r = (float) w / float (h);
        gluPerspective( camera.fov, r, 0.02, 1000.0);
 
        glMatrixMode( GL_MODELVIEW );
        glLoadIdentity();
    }
}
 
//-------------------------------------------------------------------------------
 
void DrawNode( Node *node )
{
    glPushMatrix();
 
    const Material *mtl = node->GetMaterial();
    if ( mtl ) mtl->SetViewportMaterial();
 
    Matrix3 tm = node->GetTransform();
    Point3 p = node->GetPosition();
    float m[16] = { tm[0],tm[1],tm[2],0, tm[3],tm[4],tm[5],0, tm[6],tm[7],tm[8],0, p.x,p.y,p.z,1 };
    glMultMatrixf( m );
 
    Object *obj = node->GetObject();
    if ( obj ) obj->ViewportDisplay();
 
    for ( int i=0; i<node->GetNumChild(); i++ ) {
        DrawNode( node->GetChild(i) );
    }
 
    glPopMatrix();
}
 
//-------------------------------------------------------------------------------
 
void DrawScene(bool capture=false)
{
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT );
 
    const TextureMap *bgMap = background.GetTexture();
    if ( bgMap ) {
        glDepthMask(GL_FALSE);
        glMatrixMode(GL_PROJECTION);
        glPushMatrix();
        glLoadIdentity();
        glMatrixMode(GL_MODELVIEW);
        Color c = background.GetColor();
        glColor3f(c.r,c.g,c.b);
        if ( bgMap->SetViewportTexture() ) {
            glEnable( GL_TEXTURE_2D );
            glMatrixMode( GL_TEXTURE );
            Matrix3 tm = bgMap->GetInverseTransform();
            Point3 p = tm * bgMap->GetPosition();
            float m[16] = { tm[0],tm[1],tm[2],0, tm[3],tm[4],tm[5],0, tm[6],tm[7],tm[8],0, -p.x,-p.y,-p.z,1 };
            glLoadMatrixf( m );
            glMatrixMode( GL_MODELVIEW );
        } else {
            glDisable( GL_TEXTURE_2D );
        }
        glBegin(GL_QUADS);
        glTexCoord2f(0,1);
        glVertex2f(-1,-1);
        glTexCoord2f(1,1);
        glVertex2f( 1,-1);
        glTexCoord2f(1,0);
        glVertex2f( 1, 1);
        glTexCoord2f(0,0);
        glVertex2f(-1, 1);
        glEnd();
        glMatrixMode(GL_PROJECTION);
        glPopMatrix();
        glMatrixMode(GL_MODELVIEW);
        glDepthMask(GL_TRUE);
 
        glDisable( GL_TEXTURE_2D );
    }
 
    glEnable( GL_LIGHTING );
    glEnable( GL_DEPTH_TEST );
 
    glPushMatrix();
    Point3 p = camera.pos;
    Point3 t = camera.pos + camera.dir;
    Point3 u = camera.up;
    gluLookAt( p.x, p.y, p.z,  t.x, t.y, t.z,  u.x, u.y, u.z );
 
    for ( unsigned int i=0; i<lights.size(); i++ ) {
        lights[i]->SetViewportLight(i);
    }
 
    DrawNode(&rootNode);
 
    glPopMatrix();
 
    glDisable( GL_DEPTH_TEST );
    glDisable( GL_LIGHTING );
    glDisable( GL_TEXTURE_2D );
 
    if ( capture ) {
        glReadPixels( 0, 0, camera.imgWidth, camera.imgHeight, GL_RGB, GL_UNSIGNED_BYTE, renderImage.GetPixels() );
    }
}
 
//-------------------------------------------------------------------------------
 
void DrawProgressBar()
{
    int rp = renderImage.GetNumRenderedPixels();
    int np = renderImage.GetWidth() * renderImage.GetHeight();
    if ( rp >= np ) return;
 
    float done = (float) rp / (float) np;
 
    glMatrixMode( GL_PROJECTION );
    glPushMatrix();
    glLoadIdentity();
 
    glBegin(GL_LINES);
    glColor3f(1,1,1);
    glVertex2f(-1,-1);
    glVertex2f(done*2-1,-1);
    glColor3f(0,0,0);
    glVertex2f(done*2-1,-1);
    glVertex2f(1,-1);
    glEnd();
 
    glPopMatrix();
    glMatrixMode( GL_MODELVIEW );
}
 
//-------------------------------------------------------------------------------
 
void GlutDisplay()
{
    switch ( viewMode ) {
    case VIEWMODE_OPENGL:
        DrawScene();
        break;
    case VIEWMODE_IMAGE:
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT );
        glDrawPixels( renderImage.GetWidth(), renderImage.GetHeight(), GL_RGB, GL_UNSIGNED_BYTE, renderImage.GetPixels() );
        DrawProgressBar();
        break;
    case VIEWMODE_Z:
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT );
        if ( ! renderImage.GetZBufferImage() ) renderImage.ComputeZBufferImage();
        glDrawPixels( renderImage.GetWidth(), renderImage.GetHeight(), GL_LUMINANCE, GL_UNSIGNED_BYTE, renderImage.GetZBufferImage() );
        break;
    }
 
    glutSwapBuffers();
}
 
//-------------------------------------------------------------------------------
 
void GlutIdle()
{
    static int lastRenderedPixels = 0;
    if ( mode == MODE_RENDERING ) {
        int nrp = renderImage.GetNumRenderedPixels();
        if ( lastRenderedPixels != nrp ) {
            lastRenderedPixels = nrp;
            if ( renderImage.IsRenderDone() ) {
                mode = MODE_RENDER_DONE;
                int endTime = (int) time(NULL);
                int t = endTime - startTime;
                int h = t / 3600;
                int m = (t % 3600) / 60;
                int s = t % 60;
                printf("\nRender time is %d:%02d:%02d.\n",h,m,s);
            }
            glutPostRedisplay();
        }
    }
}
 
//-------------------------------------------------------------------------------
 
void GlutKeyboard(unsigned char key, int x, int y)
{
    switch ( key ) {
    case 27:    // ESC
        exit(0);
        break;
    case ' ':
        switch ( mode ) {
        case MODE_READY: 
            mode = MODE_RENDERING;
            viewMode = VIEWMODE_IMAGE;
            DrawScene(true);
            startTime = time(NULL);
            BeginRender();
            break;
        case MODE_RENDERING:
            mode = MODE_READY;
            StopRender();
            glutPostRedisplay();
            break;
        case MODE_RENDER_DONE: 
            mode = MODE_READY;
            viewMode = VIEWMODE_OPENGL;
            glutPostRedisplay();
            break;
        }
        break;
    case '1':
        viewMode = VIEWMODE_OPENGL;
        glutPostRedisplay();
        break;
    case '2':
        viewMode = VIEWMODE_IMAGE;
        glutPostRedisplay();
        break;
    case '3':
        viewMode = VIEWMODE_Z;
        glutPostRedisplay();
        break;  }
}
 
//-------------------------------------------------------------------------------
 
void PrintPixelData(int x, int y)
{
    if ( x < renderImage.GetWidth() && y < renderImage.GetHeight() ) {
        Color24 *colors = renderImage.GetPixels();
        float *zbuffer = renderImage.GetZBuffer();
        int i = y*renderImage.GetWidth() + x;
        printf("Pixel [ %d, %d ] Color24: %d, %d, %d   Z: %f\n", x, y, colors[i].r, colors[i].g, colors[i].b, zbuffer[i] );
    } else {
        printf("-- Invalid pixel (%d,%d) --\n",x,y);
    }
}
 
//-------------------------------------------------------------------------------
 
void GlutMouse(int button, int state, int x, int y)
{
    if ( state == GLUT_DOWN ) {
        PrintPixelData(x,y);
    }
}
 
//-------------------------------------------------------------------------------
 
void GlutMotion(int x, int y)
{
    PrintPixelData(x,y);
}
 
//-------------------------------------------------------------------------------