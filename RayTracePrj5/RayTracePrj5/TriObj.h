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
#include "cyTriMesh.h"

using namespace std;

class TriObj : public Object, private cyTriMesh
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
   
    
    bool IntersectRay( const Ray &ray, HitInfo& hInfo, int hitSide=HIT_FRONT )const
    {
        bool hit = false;
        //get bound box
        Box aabb = GetBoundBox();
        
        //check interseciton with bound box
        bool boxHit = aabb.IntersectRay(ray, BIGFLOAT);
        //then check for interseciton with each triangle - already implemented
        if(boxHit){
            Point3 P(ray.p.x,ray.p.y,ray.p.z);
            Point3 D(ray.dir.x,ray.dir.y,ray.dir.z);
            Point3 N;
            
            for(unsigned int i = 0; i < NF(); i++){
                if(IntersectTriangle(ray, hInfo, hitSide, i)){
                    hit = true;
                    
                }
            }
        }
        return hit;
        
    }
    
    void ViewportDisplay() const
    {
        glBegin(GL_TRIANGLES);
        for ( unsigned int i=0; i< NF(); i++ ) {
            for ( int j=0; j<3; j++ ) {
                if ( HasTextureVertices() ) glTexCoord3fv( &VT( FT(i).v[j] ).x );
                if ( HasNormals() ) glNormal3fv( &VN( FN(i).v[j] ).x );
                glVertex3fv( &V( F(i).v[j] ).x );
            }
        }
        glEnd();
    }
    
    Box GetBoundBox() const { return Box(GetBoundMin(),GetBoundMax()); }
    
	bool Load(const char *filename)
	{
		if ( ! LoadFromFileObj( filename ) ) return false;
		if ( ! HasNormals() ) ComputeNormals();
		ComputeBoundingBox();
		return true;

    }
private:
    bool IntersectTriangle( const Ray &ray, HitInfo &hInfo, int hitSide, unsigned int faceID ) const
    {
        //get the vertices of each face with the faceID
        Point3 A, B, C, P, D, N, H;
        float t, a, a1, a2, alpha, beta;
        float bias = 1e-7f;
        Point2 A2d, B2d, C2d, H2d;
        A = V(F(faceID).v[0]);
        B = V(F(faceID).v[1]);
        C = V(F(faceID).v[2]);
        
        //Do Ray triangle intersection with the three vertices thus obtained
        P = ray.p;
        D = ray.dir;
        Point3 nTemp = (B - A).Cross(C - A);
        //N = nTemp/ nTemp.Length();
        nTemp.Normalize();
        N = nTemp;

        if(N.Dot(P-A) < bias) return false; /* we are in the back face */
        if(N.Dot(D) == 0) return false; /* ray is parallel */
        t = N.Dot(P-A)/-(N.Dot(D));
        
        //if the t-value is less than what is there in hit info then assign and signal hit
        if(t < hInfo.z && t >= bias && t< BIGFLOAT){
            H = P + t * D; /* point of intersection on the plane */
            int projAxis = 0; //0-xaxis, 1-yaxis, 2-zaxis
            //project H onto an axis aligned plane - find the maximum component
            if(abs(H.x) > abs(H.y)){
                if(abs(H.x) > abs(H.z)){
                    projAxis = 0;
                }
                else{
                    projAxis = 2;
                }
            }
            else if(abs(H.y) > abs(H.z)){
                projAxis = 1;
            }
            else{
                projAxis = 2;
            }
            //cout<<projAxis<<endl;
            switch (projAxis) {
                case 0: /* project onto the x-axis*/
                    A2d = Point2(A.y, A.z);
                    B2d = Point2(B.y, B.z);
                    C2d = Point2(C.y, C.z);
                    H2d = Point2(H.y, H.z);
                    break;
                case 1:/* project onto the y-axis*/
                    A2d = Point2(A.x, A.z);
                    B2d = Point2(B.x, B.z);
                    C2d = Point2(C.x, C.z);
                    H2d = Point2(H.x, H.z);
                    break;
                case 2:/* project onto the z-axis*/
                    A2d = Point2(A.x, A.y);
                    B2d = Point2(B.x, B.y);
                    C2d = Point2(C.x, C.y);
                    H2d = Point2(H.x, H.y);
                default:
                    break;
            }
            
            a = (B2d - A2d).Cross(C2d - A2d)/ 2.0;
            a1 = (B2d - A2d).Cross(H2d - A2d)/ 2.0;
            a2 = (H2d - A2d).Cross(C2d - A2d)/ 2.0;
            alpha = a1/a;
            beta = a2/a;
            if(alpha < -bias || beta < -bias || alpha + beta > 1+bias){
                return false; /* point H is outside triangle */
            }
//            if((a <0 && a1 < 0 && a2 <0) || (a >=0 && a1 >= 0 && a2 >=0)){ /* sign check */
//                if(a1 + a2 <= a){ /* H is inside */
//                    
//                    float b0 = 1 - beta - alpha;
//                    float b1 = beta;
//                    float b2 = alpha;
//                    H = b0 * A + b1 * B + b2 * C;
//                    
//                    hInfo.front = true;
//                    hInfo.p = H;
//                    hInfo.N = N;
//                    hInfo.z = t;
//                    return true;
//                }
//            }
            float b0 = 1 - beta - alpha;
            float b1 = beta;
            float b2 = alpha;
            H = b0 * A + b1 * B + b2 * C;
//            Point3 Na = (B - A).Cross(C - A)/((B - A).Cross(C - A)).Length();
//            Point3 Nb = (C-B).Cross(A - B)/((C-B).Cross(A - B)).Length();
//            Point3 Nc = (A-C).Cross(B - C)/((A-C).Cross(B - C)).Length();
//            Point3 Na = VN(F(faceID).v[0]);
//            Point3 Nb = VN(F(faceID).v[1]);
//            Point3 Nc = VN(F(faceID).v[2]);
            Point3 Normal = GetNormal(faceID, Point3(b0,b1,b2));
            N = Normal;
            //N = b0 * Na + b1 * Nb + b2 * Nc;
            hInfo.front = true;
            hInfo.p = H;
            hInfo.N = N;
            hInfo.z = t;
            return true;
            
        }
        
        return false;
    }
};