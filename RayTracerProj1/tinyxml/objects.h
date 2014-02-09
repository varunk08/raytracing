
//-------------------------------------------------------------------------------
///
/// \file       objects.h 
/// \author     Cem Yuksel (www.cemyuksel.com)
/// \version    7.0
/// \date       October 7, 2013
///
/// \brief Example source for CS 6620 - University of Utah.
///
//-------------------------------------------------------------------------------
 
#ifndef _OBJECTS_H_INCLUDED_
#define _OBJECTS_H_INCLUDED_
 
#include "scene.h"
#include "cyTriMesh.h"
#include "cyBVH.h"
 
//-------------------------------------------------------------------------------
 
class Sphere : public Object
{
public:
    virtual bool IntersectRay( const Cone &ray, HitInfo &hInfo, int hitSide=HIT_FRONT ) const;
    virtual void ViewportDisplay() const;
    virtual Box GetBoundBox() const { return Box(-1,-1,-1,1,1,1); }
};
 
extern Sphere theSphere;
 
//-------------------------------------------------------------------------------
 
class Plane : public Object
{
public:
    virtual bool IntersectRay( const Cone &ray, HitInfo &hInfo, int hitSide=HIT_FRONT ) const;
    virtual void ViewportDisplay() const;
    virtual Box GetBoundBox() const { return Box(-1,-1,0,1,1,0); }
};
 
extern Plane thePlane;
 
//-------------------------------------------------------------------------------
 
class TriObj : public Object, private cyTriMesh
{
public:
    virtual bool IntersectRay( const Cone &ray, HitInfo &hInfo, int hitSide=HIT_FRONT ) const;
    virtual void ViewportDisplay() const;
    virtual Box GetBoundBox() const { return Box(GetBoundMin(),GetBoundMax()); }
 
    bool Load(const char *filename)
    {
        bvh.Clear();
        if ( ! LoadFromFileObj( filename ) ) return false;
        if ( ! HasNormals() ) ComputeNormals();
        ComputeBoundingBox();
        bvh.SetMesh(this,4);
        return true;
    }
 
private:
    cyBVHTriMesh bvh;
    bool IntersectTriangle( const Cone &ray, HitInfo &hInfo, int hitSide, unsigned int faceID ) const;
    bool TraceBVHNode( const Cone &ray, HitInfo &hInfo, int hitSide, unsigned int nodeID ) const;
};
 
//-------------------------------------------------------------------------------
 
#endif