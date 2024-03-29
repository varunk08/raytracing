#include "scene.h"
#include <iostream>
#include <string.h>
#include <strings.h>
#include <cmath>
#include <GLUT/GLUT.h>

#define USE_MATH_DEFINES
using namespace std;

class Sphere : public Object
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
    }
    Box GetBoundBox() const { return Box(-1,-1,-1,1,1,1); }    
    bool IntersectRay( const Cone &ray, HitInfo &hInfo, int hitSide=HIT_FRONT )const
    {
        
        
        float r=1.0;
        Point3 pvec(0,0,0),dvec(0,0,0),pOrigin(0,0,0);
        
        //std::cout<<"Before ("<<ray.p.x<<","<<ray.p.y<<","<<ray.p.z<<")"<<std::endl;
        
        Ray tRay = ray;
        pvec.Set(tRay.p.x,tRay.p.y,tRay.p.z);
        dvec.Set(tRay.dir.x,tRay.dir.y,tRay.dir.z);
        
        //std::cout<<"After: ("<<tRay.dir.x<<","<<tRay.dir.y<<","<<tRay.dir.z<<")"<<std::endl;
        
        float a = dvec.Dot(dvec);
        float b = (2.0*(dvec.Dot(pvec - pOrigin)));
        float c = (pvec - pOrigin).Dot(pvec - pOrigin) - r*r;
        
        float det = b*b - 4.0*(a * c);
        // cout<<name<<" "<<t<<" "<<a<<" "<<b<<" "<<c<< endl;
        //cout<<"T: "<<t<<endl;
        if (det>=0) {
            
            float t1 = (-b + sqrtf(det))/(2.0 * a);  //cout<<"T1: :"<<t1<<endl;
            float t2 = (-b - sqrtf(det))/(2.0 * a);  //cout<<"T2: :"<<t2<<endl;
            float temp;
            if ( t1 < t2) { temp = t1; }
            if (t2 < t1) { temp = t2; }
            if (t1 == t2){temp = t1; }
            switch (hitSide) {
                    /*case HIT_FRONT:
                     if (hInfo.z > temp && temp >0){
                     hInfo.z = temp;
                     Point3 P(pvec.x+hInfo.z*dvec.x,pvec.y+hInfo.z*dvec.y,pvec.z+hInfo.z*dvec.z);
                     Point3 N = P;
                     hInfo.p = P;
                     hInfo.N = N;
                     //hInfo.front = true;
                     return true;
                     }
                     break;*/
                case HIT_FRONT:
                case HIT_FRONT_AND_BACK:
                    if (hInfo.z > temp){
                        float newZ = 0.0;
                        if( t1 <0 && t2>0){
                            newZ=t2; /*selecting the least value viz. > 1*/
                            hInfo.front = false;
                        }
                        if(t2 < 0 && t1>0){
                            newZ=t1;
                            hInfo.front = false;
                        }
                        if(t1 <0 && t2 <0) return false;
                        if(t1 >0 && t2 >0){
                            if(t1 < t2) newZ = t1;
                            else newZ = t2;
                            hInfo.front = true;
                        }
                        temp = newZ;
                        hInfo.z = temp;
                        Point3 P(pvec.x+hInfo.z*dvec.x,pvec.y+hInfo.z*dvec.y,pvec.z+hInfo.z*dvec.z);
                        Point3 N = P;
                        hInfo.p = P;
                        hInfo.N = N;
                        
                        /* calculate uvw */
                        float u = 0.5 - atan2(P.y,P.x)/ (2 * M_PI);
                        float v = 0.5 - asin(P.z)/M_PI;
                        hInfo.uvw = Point3(u,v,0);
                        Point3 minorAxis;
                        Point3 majorAxis;
                        ray.EllipseAt(hInfo.z, hInfo.N, majorAxis, minorAxis);
                        hInfo.duvw[0] = minorAxis * 2;
                        hInfo.duvw[1] = majorAxis * 2;
                        return true;
                    }
                    break;
                default:
                    break;
            }
            
        }
        
        return false;
        
        
    }
    
    
void ViewportDisplay() const
{
  static GLUquadric *q = NULL;
  if ( q == NULL ) {
      q = gluNewQuadric();
      gluQuadricTexture(q,true);
  }
  gluSphere(q,1,50,50);
}

};
extern Sphere theSphere;