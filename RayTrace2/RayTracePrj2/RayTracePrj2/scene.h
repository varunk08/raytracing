//-------------------------------------------------------------------------------
///
/// \file       scene.h
/// \author     Cem Yuksel (www.cemyuksel.com)
/// \version    2.0
/// \date       September 2, 2013
///
/// \brief Example source for CS 6620 - University of Utah.
///
//-------------------------------------------------------------------------------

#ifndef _SCENE_H_INCLUDED_
#define _SCENE_H_INCLUDED_

//-------------------------------------------------------------------------------

#include <string.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <vector>

#include "cyPoint.h"
typedef cyPoint2f Point2;
typedef cyPoint3f Point3;
typedef cyPoint4f Point4;

#include "cyMatrix3.h"
typedef cyMatrix3f Matrix3;

#include "cyColor.h"
typedef cyColor Color;
typedef cyColorA ColorA;
typedef cyColor24 Color24;

typedef unsigned char uchar;

//-------------------------------------------------------------------------------

#ifndef min
# define min(a,b) ((a)<(b)?(a):(b))
#endif

#ifndef max
# define max(a,b) ((a)>(b)?(a):(b))
#endif

#define BIGFLOAT 1.0e30f

//-------------------------------------------------------------------------------

using namespace std;

class Ray
{
public:
	Point3 p, dir;
    
	Ray() {}
	Ray(const Point3 &_p, const Point3 &_dir) : p(_p), dir(_dir) {}
	void Normalize() { dir.Normalize(); }
};

//-------------------------------------------------------------------------------

class Node;

#define HIT_NONE			0
#define HIT_FRONT			1
#define HIT_BACK			2
#define HIT_FRONT_AND_BACK	(HIT_FRONT|HIT_BACK)

struct HitInfo
{
	float z;			// the distance from the ray center to the hit point
	Point3 p;			// position of the hit point
	Point3 N;			// surface normal at the hit point
	const Node *node;	// the object node that was hit
	bool front;			// true if the ray hits the front side, false if the ray hits the back side
    
	HitInfo() { Init(); }
	void Init() { z=BIGFLOAT; node=NULL; front=true; }
};

//-------------------------------------------------------------------------------

class ItemBase
{
private:
	char *name;					// The name of the item
    
public:
	ItemBase() : name(NULL) {}
	virtual ~ItemBase() { if ( name ) delete [] name; }
    
	const char* GetName() const { return name ? name : ""; }
	void SetName(const char *newName)
	{
		if ( name ) delete [] name;
		if ( newName ) {
			int n = strlen(newName);
			name = new char[n+1];
			for ( int i=0; i<n; i++ ) name[i] = newName[i];
			name[n] = '\0';
		} else { name = NULL; }
	}
};

template <class T> class ItemList : public std::vector<T*>
{
public:
	virtual ~ItemList() { DeleteAll(); }
	void DeleteAll() { int n=(int)this->size(); for ( int i=0; i<n; i++ ) if ( this->at(i) ) delete this->at(i); }
};


template <class T> class ItemFileList
{
public:
	void Clear() { list.DeleteAll(); }
	void Append( T* item, const char *name ) { list.push_back( new FileInfo(item,name) ); }
	T* Find( const char *name ) const { int n=list.size(); for ( int i=0; i<n; i++ ) if ( list[i] && strcmp(name,list[i]->GetName())==0 ) return list[i]->GetObj(); return NULL; }
    
private:
	class FileInfo : public ItemBase
	{
	private:
		T *item;
	public:
		FileInfo() : item(NULL) {}
		FileInfo(T *_item, const char *name) : item(_item) { SetName(name); }
		~FileInfo() { Delete(); }
		void Delete() { if (item) delete item; item=NULL; }
		void SetObj(T *_item) { Delete(); item=_item; }
		T* GetObj() { return item; }
	};
    
	ItemList<FileInfo> list;
};

//-------------------------------------------------------------------------------

class Transformation
{
private:
	Matrix3 tm;						// Transformation matrix to the local space
	Point3 pos;						// Translation part of the transformation matrix
	mutable Matrix3 itm;			// Inverse of the transformation matrix (cached)
public:
	Transformation() : pos(0,0,0) { tm.SetIdentity(); itm.SetIdentity(); }
	const Matrix3& GetTransform() const { return tm; }
	const Point3& GetPosition() const { return pos; }
	const Matrix3&	GetInverseTransform() const { return itm; }
    
	Point3 TransformTo( const Point3 &p ) const { return itm * (p - pos); }	// Transform to the local coordinate system
	Point3 TransformFrom( const Point3 &p ) const { return tm*p + pos; }	// Transform from the local coordinate system
    
	// Transforms a vector to the local coordinate system (same as multiplication with the inverse transpose of the transformation)
	Point3 VectorTransformTo( const Point3 &dir ) const { return TransposeMult(tm,dir); }
    
	// Transforms a vector from the local coordinate system (same as multiplication with the inverse transpose of the transformation)
	Point3 VectorTransformFrom( const Point3 &dir ) const { return TransposeMult(itm,dir); }
    
	void Translate(Point3 p) { pos+=p; }
	void Rotate(Point3 axis, float degree) { Matrix3 m; m.SetRotation(axis,degree*(float)M_PI/180.0f); Transform(m); }
	void Scale(float sx, float sy, float sz) { Matrix3 m; m.Zero(); m[0]=sx; m[4]=sy; m[8]=sz; Transform(m); }
	void Transform(const Matrix3 &m) { tm*=m; pos=m*pos; tm.GetInverse(itm); }
    
	void InitTransform() { pos.Zero(); tm.SetIdentity(); itm.SetIdentity(); }
    
private:
	// Multiplies the given vector with the transpose of the given matrix
	static Point3 TransposeMult( const Matrix3 &m, const Point3 &dir )
	{
		Point3 d;
		d.x = m.GetColumn(0) % dir;
		d.y = m.GetColumn(1) % dir;
		d.z = m.GetColumn(2) % dir;
		return d;
	}
};

//-------------------------------------------------------------------------------

// Base class for all object types
class Object
{
public:
	virtual bool IntersectRay( const Ray &ray, HitInfo &hInfo, int hitSide=HIT_FRONT ) const=0;
	virtual void ViewportDisplay() const {}	// used for OpenGL display
};

typedef ItemFileList<Object> ObjFileList;

//-------------------------------------------------------------------------------

class Light : public ItemBase
{
public:
	virtual Color	Illuminate(Point3 p) const=0;
	virtual Point3	Direction (Point3 p) const=0;
	virtual bool	IsAmbient () const { return false; }
	virtual void	SetViewportLight(int lightID) const {}	// used for OpenGL display
};

class LightList : public ItemList<Light> {};

//-------------------------------------------------------------------------------

class Material : public ItemBase
{
public:
	// The main method that handles the shading by calling all the lights in the list.
	// ray: incoming ray and the ray differentials,
	// hInfo: hit information for the point that is being shaded, lights: the light list,
	virtual Color Shade(const Ray &ray, const HitInfo &hInfo, const LightList &lights) const=0;
    
	virtual void SetViewportMaterial() const {}	// used for OpenGL display
   // virtual ~Material(){};
};

class MaterialList : public ItemList<Material>
{
public:
	Material* Find( const char *name ) { int n=size(); for ( int i=0; i<n; i++ ) if ( at(i) && strcmp(name,at(i)->GetName())==0 ) return at(i); return NULL; }
};

//-------------------------------------------------------------------------------

class Node : public ItemBase, public Transformation
{
private:
	Node **child;				// Child nodes
	int numChild;				// The number of child nodes
	Object *obj;				// Object reference (merely points to the object, but does not own the object, so it doesn't get deleted automatically)
	Material *mtl;				// Material used for shading the object
public:
	Node() : child(NULL), numChild(0), obj(NULL) {}
    
	void Init() { DeleteAllChildNodes(); obj=NULL; SetName(NULL); InitTransform(); } // Initialize the node deleting all child nodes
    
	// Hierarchy management
	int	 GetNumChild() const { return numChild; }
	void SetNumChild(int n, int keepOld=false)
	{
		if ( n < 0 ) n=0;	// just to be sure
		Node **nc = NULL;	// new child pointer
		if ( n > 0 ) nc = new Node*[n];
		for ( int i=0; i<n; i++ ) nc[i] = NULL;
		if ( keepOld ) {
			int sn = min(n,numChild);
			for ( int i=0; i<sn; i++ ) nc[i] = child[i];
		}
		if ( child ) delete [] child;
		child = nc;
		numChild = n;
	}
	const Node*	GetChild(int i) const		{ return child[i]; }
	Node*		GetChild(int i)				{ return child[i]; }
	void		SetChild(int i, Node *node)	{ child[i]=node; }
	void		AppendChild(Node *node)		{ SetNumChild(numChild+1,true); SetChild(numChild-1,node); }
	void		RemoveChild(int i)			{ for ( int j=i; j<numChild-1; j++) child[j]=child[j-1]; SetNumChild(numChild-1); }
	void		DeleteAllChildNodes()		{ for ( int i=0; i<numChild; i++ ) { child[i]->DeleteAllChildNodes(); delete child[i]; } SetNumChild(0); }
    
	// Object management
	const Object*	GetObject() const { return obj; }
	Object*			GetObject() { return obj; }
	void			SetObject(Object *object) { obj=object; }
    
	// Material management
	const Material* GetMaterial() const { return mtl; }
	void			SetMaterial(Material *material) { mtl=material; }
    
	// Transformations
	Ray ToNodeCoords( const Ray &ray ) const
	{
		Ray r;
		r.p   = TransformTo(ray.p);
		r.dir = TransformTo(ray.p + ray.dir) - r.p;
		return r;
	}
	void FromNodeCoords( HitInfo &hInfo ) const
	{
		hInfo.p = TransformFrom(hInfo.p);
		hInfo.N = VectorTransformFrom(hInfo.N).GetNormalized();
	}
    
};

//-------------------------------------------------------------------------------

class Camera
{
public:
	Point3 pos, dir, up;
	float fov;
	int imgWidth, imgHeight;
    
	void Init()
	{
		pos.Set(0,0,0);
		dir.Set(0,0,-1);
		up.Set(0,1,0);
		fov = 40;
		imgWidth = 200;
		imgHeight = 150;
	}
};

//-------------------------------------------------------------------------------

class RenderImage
{
private:
	Color24	*img;
	float	*zbuffer;
	uchar	*zbuffer8;
	int		width, height;
	int		numRenderedPixels;
public:
	RenderImage() : img(NULL), zbuffer(NULL), zbuffer8(NULL), width(0), height(0), numRenderedPixels(0) {}
	void Init(int w, int h)
	{
		width=w;
		height=h;
		if (img) delete [] img;
		img = new Color24[width*height];
		if (zbuffer) delete [] zbuffer;
		zbuffer = new float[width*height];
		if (zbuffer8) delete [] zbuffer8;
		zbuffer8 = NULL;
		ResetNumRenderedPixels();
	}
	int			GetWidth() const	{ return width; }
	int			GetHeight() const	{ return height; }
	Color24*	GetPixels()			{ return img; }
	float*		GetZBuffer()		{ return zbuffer; }
	uchar*		GetZBufferImage()	{ return zbuffer8; }
    
	void	ResetNumRenderedPixels()		{ numRenderedPixels=0; }
	int		GetNumRenderedPixels() const	{ return numRenderedPixels; }
	void	IncrementNumRenderPixel(int n)	{ numRenderedPixels+=n; }	// not thread-safe
	bool	IsRenderDone() const			{ return numRenderedPixels >= width*height; }
    //custom
    void PutPixel(int index, Color24 color,float zbuf)
    {
        img[index]=color;
        zbuffer[index] = zbuf;
        //cout<<"zbuf: "<<index<<endl;
        IncrementNumRenderPixel(1);
    }
	void	ComputeZBufferImage()
	{
		int size = width * height;
		if (zbuffer8) delete [] zbuffer8;
		zbuffer8 = new unsigned char[size];
        
		float zmin=BIGFLOAT, zmax=0;
		for ( int i=0; i<size; i++ ) {
            //cout<<"zbuf: "<<zbuffer[i]<<endl;
			if ( zbuffer[i] == BIGFLOAT ) continue;
			if ( zmin > zbuffer[i] ) zmin = zbuffer[i];
			if ( zmax < zbuffer[i] ) zmax = zbuffer[i];
		}
        cout<<"Zmin: "<<zmin<<" Zmax: "<<zmax<<endl;
		for ( int i=0; i<size; i++ ) {
			if ( zbuffer[i] == BIGFLOAT ) {
                //cout<<"Hurray"<<endl;
                zbuffer8[i] = 0;}
			else {
				float f = (zmax-zbuffer[i])/(zmax-zmin);
				int c = int(f * 255);
				if ( c < 0 ) f = 0;
				if ( c > 255 ) f = 255;
				zbuffer8[i] = c;
                
			}
		}
	}
    
	bool SaveImage (const char *filename, bool flipped=false) const { return SavePPM(filename,&img[0].r,3,flipped); }
	bool SaveZImage(const char *filename, bool flipped=false) const { return SavePPM(filename,zbuffer8, 1,flipped); }
    
private:
	bool SavePPM(const char *filename, uchar *data, int compCount, bool flipped) const
	{
		FILE *fp = fopen(filename,"wb");
		if ( !fp ) return false;
		fprintf(fp,"P6\n%d %d\n255\n", width, height);
		switch( compCount ) {
            case 1:
                if ( flipped ) {
                    for ( int y=height-1; y>=0; y-- ) {
                        int i = y*width;
                        for ( int x=0; x<width; x++, i++ ) {
                            uchar d[3] = { data[i], data[i], data[i] };
                            fwrite(d,3,1,fp);
                        }
                    }
                } else {
                    for ( int i=0; i<width*height; i++ ) {
                        uchar d[3] = { data[i], data[i], data[i] };
                        fwrite(d,3,1,fp);
                    }
                }
                break;
            case 3:
                if ( flipped ) {
                    for ( int y=height-1; y>=0; y-- ) {
                        fwrite(&data[y*width*3],3,width,fp);
                    }
                } else {
                    fwrite(data,3,width*height,fp);
                }
                break;
		}
		fclose(fp);
		return true;
	}
};

//-------------------------------------------------------------------------------

#endif