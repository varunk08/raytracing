//-------------------------------------------------------------------------------
///
/// \file       materials.h
/// \author     Cem Yuksel (www.cemyuksel.com)
/// \version    1.0
/// \date       September 2, 2013
///
/// \brief Example source for CS 6620 - University of Utah.
///
//-------------------------------------------------------------------------------

#ifndef _MATERIALS_H_INCLUDED_
#define _MATERIALS_H_INCLUDED_
//#include "sphere.h"
#include "scene.h"

//-------------------------------------------------------------------------------

class MtlBlinn : public Material
{
public:
	MtlBlinn() : diffuse(0.5f,0.5f,0.5f), specular(0.7f,0.7f,0.7f), glossiness(20.0f) {}
    virtual Color Shade(const Ray &ray, const HitInfo &hInfo, const LightList &lights) const;
    
	void SetDiffuse(Color dif) { diffuse = dif; }
	void SetSpecular(Color spec) { specular = spec; }
	void SetGlossiness(float gloss) { glossiness = gloss; }
    
    virtual void SetViewportMaterial()const;// used for OpenGL display
    
     ~MtlBlinn(){};
    
private:
	Color diffuse, specular;
	float glossiness;
};

//-------------------------------------------------------------------------------

#endif