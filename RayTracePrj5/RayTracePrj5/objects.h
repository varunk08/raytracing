//-------------------------------------------------------------------------------
///
/// \file       objects.h
/// \author     Cem Yuksel (www.cemyuksel.com)
/// \version    5.0
/// \date       September 23, 2013
///
/// \brief Example source for CS 6620 - University of Utah.
///
//-------------------------------------------------------------------------------

#ifndef _OBJECTS_H_INCLUDED_
#define _OBJECTS_H_INCLUDED_

#include "scene.h"
#include "cyTriMesh.h"
//-------------------------------------------------------------------------------

//class Sphere : public Object
//{
//public:
//	virtual bool IntersectRay( const Ray &ray, HitInfo &hInfo, int hitSide=HIT_FRONT ) const;
//	virtual void ViewportDisplay() const;
//	virtual Box GetBoundBox() const { return Box(-1,-1,-1,1,1,1); }
//};

//extern Sphere theSphere;

//-------------------------------------------------------------------------------

//class Plane : public Object
//{
//public:
//	virtual bool IntersectRay( const Ray &ray, HitInfo &hInfo, int hitSide=HIT_FRONT ) const;
//	virtual void ViewportDisplay() const;
//	virtual Box GetBoundBox() const { return Box(-1,-1,0,1,1,0); }
//};
//
//extern Plane thePlane;

//-------------------------------------------------------------------------------
//
//class TriObj : public Object, private cyTriMesh
//{
//public:
//	virtual bool IntersectRay( const Ray &ray, HitInfo &hInfo, int hitSide=HIT_FRONT ) const;
//	virtual void ViewportDisplay() const;
//    virtual Box GetBoundBox() const { return Box(GetBoundMin(),GetBoundMax()); }
//bool Load(const char *filename)
//{
//    if ( ! LoadFromFileObj( filename ) ) return false;
//    if ( ! HasNormals() ) ComputeNormals();
//    ComputeBoundingBox();
//    return true;
//}
//private:
//	bool IntersectTriangle( const Ray &ray, HitInfo &hInfo, int hitSide, unsigned int faceID ) const;
//};

//-------------------------------------------------------------------------------

#endif