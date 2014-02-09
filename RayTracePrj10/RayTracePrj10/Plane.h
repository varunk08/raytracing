//
//  Plane.h
//  RayTracePrj5
//
//  Created by Varun kumar Karuppannan on 28/09/13.
//  Copyright (c) 2013 Varun kumar Karuppannan. All rights reserved.
//
#include "scene.h"
#include <iostream>
#include <string.h>
#include <strings.h>
#include <cmath>
#include <GLUT/GLUT.h>

using namespace std;

class Plane : public Object
{
    
    
private:
    char* name;
public:
    
    void setName(const char* newName){
        cout<<newName<<endl;
        if ( newName ) {
            int n = strlen(newName);
            name = new char[n+1];
            for ( int i=0; i<n; i++ ) name[i] = newName[i];
            name[n] = '\0';
        } else { name = NULL; }
    };
    
    Box GetBoundBox() const { return Box(-1,-1,0,1,1,0); }
    
    bool IntersectRay( const Cone &ray, HitInfo &hInfo, int hitSide=HIT_FRONT )const
    {
        

        Point3 P(ray.p.x,ray.p.y,ray.p.z);
        Point3 D(ray.dir.x,ray.dir.y,ray.dir.z);
        Point3 N(0,0,1);
        
        float t = -(P.z/D.z);
        if(t >= 0 && t < BIGFLOAT && t<hInfo.z){
            Point3 pHit(P.x + t * D.x, P.y + t * D.y, P.z + t * D.z);
            if(pHit.x >= -1 && pHit.x <= 1 && pHit.y >=-1 && pHit.y <=1){
                hInfo.z = t;
//                cout<<"pHit: "<<pHit.x<<" pHit: "<<pHit.y<<" pHit: "<<pHit.z<<endl;
                hInfo.p.Set(pHit.x, pHit.y, pHit.z);
                hInfo.N = N;
                hInfo.uvw = Point3((pHit.x + 1)/2, (pHit.y +1)/2,0);
                
                Point3 majorAxis(0,0,0), minorAxis(0,0,0);
                ray.EllipseAt(hInfo.z, hInfo.N, majorAxis, minorAxis);
//                Point3 x = pHit + minorAxis * 2;
//                Point3 y = pHit + majorAxis * 2;
                hInfo.duvw[0]= minorAxis * 2;//(pHit - x) ;
                hInfo.duvw[1]= majorAxis * 2;//(pHit - y) ;
                
                if(N.Dot(D) < 0){
                    hInfo.front = false;
                }
                else hInfo.front = true;
                return true;
            }
        }
        
//        
//        Point3 Q = Point3(1,1,0);
//        Point3 P = ray.p;
//        Point3 D = ray.dir;
//        Point3 N = Point3(0,0,1);
//        
//        float t = N.Dot(P - Q) / -(N.Dot(D));
//        
//        if(t >= bias &&t < BIGFLOAT){
//                    hInfo.z = t;
//                    hInfo.p = Point3(P.x + t * D.x, P.y + t * D.y, 0);
//                    hInfo.N = Point3(hInfo.p.x,hInfo.p.y,1);
//                    hInfo.front = true;
//                    return true;
//        }
        return false;
        
    }
    
    
void ViewportDisplay() const
{
  const int resolution = 32;
  float xyInc = 2.0f / resolution;
  float uvInc = 1.0f / resolution;
  glPushMatrix();
  glNormal3f(0,0,1);
  glBegin(GL_QUADS);
  float y1=-1, y2=xyInc-1, v1=0, v2=uvInc;
  for ( int y=0; y<resolution; y++ ) {
      float x1=-1, x2=xyInc-1, u1=0, u2=uvInc;
      for ( int x=0; x<resolution; x++ ) {
          glTexCoord2f(u1, v1);
          glVertex3f ( x1, y1, 0 );
          glTexCoord2f(u2, v1);
          glVertex3f ( x2, y1, 0 );
          glTexCoord2f(u2, v2);
          glVertex3f ( x2, y2, 0 );
          glTexCoord2f(u1, v2);
          glVertex3f ( x1, y2, 0 );
          x1=x2; x2+=xyInc; u1=u2; u2+=uvInc;   
      }
      y1=y2; y2+=xyInc; v1=v2; v2+=uvInc;
  }
  glEnd();
  glPopMatrix();
}

    
    
};
extern Plane thePlane;