#include "scene.h"
#include <iostream>
#include <string.h>
#include <strings.h>
#include <cmath>
using namespace std;
class Sphere : public Object
{
    
private:
    char* name;
public:
    
    void setName(const char* newName){
        if ( newName ) {
            int n = strlen(newName);
            name = new char[n+1];
            for ( int i=0; i<n; i++ ) name[i] = newName[i];
            name[n] = '\0';
        } else { name = NULL; }
    };
    
    
    bool IntersectRay( const Ray &ray, HitInfo& hInfo, int hitSide=HIT_FRONT )const
    {
        
      
        float r=1.0;
        Point3 pvec(0,0,0),dvec(0,0,0);
        
        //std::cout<<"Before ("<<ray.p.x<<","<<ray.p.y<<","<<ray.p.z<<")"<<std::endl;
        
        Ray tRay = ray;
        pvec.Set(tRay.p.x,tRay.p.y,tRay.p.z);
        dvec.Set(tRay.dir.x,tRay.dir.y,tRay.dir.z);
       
        //std::cout<<"After: ("<<tRay.dir.x<<","<<tRay.dir.y<<","<<tRay.dir.z<<")"<<std::endl;
        
        float a = dvec.Dot(dvec);
        float b = (2.0*(pvec.Dot(dvec)));
        float c = pvec.Dot(pvec) - r*r;
        
        float t = b*b - 4.0*(a * c);
        //std::cout<<this->name<<" "<<t<<" "<<a<<" "<<b<<" "<<c<<std::endl;
        //cout<<"T: "<<t<<endl;
        if (t>0) {
            
            float t1 = (-b + sqrtf(t))/(2.0 * a);
            float t2 = (-b - sqrtf(t))/(2.0 * a);
            float temp = hInfo.z;
            if (t1 < t2) { hInfo.z = t1; }
            if (t2 < t1) { hInfo.z = t2; }
            if (t1 == t2){ hInfo.z = t1; };
            if(hInfo.z < temp){
            //cout<<"Z: "<<hInfo.z<<" t1: "<<t1<<" t2: "<<t2<<endl;
            //p + td
            //Point3 zPt = pvec + (hInfo.z)*dvec;
            //Point3 zDist = zPt - pvec;
            //hInfo.z = zDist.Length();
            
            return true;
            }
            else hInfo.z=temp;
        }
        
        return false;

        
    }
    
    
    void ViewportDisplay() const
    {
        static GLUquadric *q = NULL;
        if ( q == NULL ) q = gluNewQuadric();
        gluSphere(q,1,50,50);
    }
    
};