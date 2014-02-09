
//-------------------------------------------------------------------------------
///
/// \file       materials.h 
/// \author     Cem Yuksel (www.cemyuksel.com)
/// \version    7.0
/// \date       October 7, 2013
///
/// \brief Example source for CS 6620 - University of Utah.
///
//-------------------------------------------------------------------------------
 
#ifndef _MATERIALS_H_INCLUDED_
#define _MATERIALS_H_INCLUDED_
 
#include "scene.h"
 
//-------------------------------------------------------------------------------
 
class MtlBlinn : public Material
{
public:
    MtlBlinn() : diffuse(0.5f,0.5f,0.5f), specular(0.7f,0.7f,0.7f), glossiness(20.0f), 
        reflection(0,0,0), refraction(0,0,0), absorption(0,0,0), ior(1), viewportTextureID(0) {}
    virtual Color Shade(const Cone &ray, const HitInfo &hInfo, const LightList &lights, int bounceCount) const;
 
    void SetDiffuse(Color dif) { diffuse.SetColor(dif); }
    void SetSpecular(Color spec) { specular.SetColor(spec); }
    void SetGlossiness(float gloss) { glossiness = gloss; }
    void SetReflection(Color reflect) { reflection.SetColor(reflect); }
    void SetRefraction(Color refract) { refraction.SetColor(refract); }
    void SetAbsorption(Color absorp ) { absorption = absorp; }
    void SetRefractionIndex(float _ior) { ior = _ior; }
 
    void SetDiffuseTexture(TextureMap *map) { diffuse.SetTexture(map); }
    void SetSpecularTexture(TextureMap *map) { specular.SetTexture(map); }
    void SetReflectionTexture(TextureMap *map) { reflection.SetTexture(map); }
    void SetRefractionTexture(TextureMap *map) { refraction.SetTexture(map); }
 
    virtual void SetViewportMaterial() const;   // used for OpenGL display
 
private:
    TexturedColor diffuse, specular, reflection, refraction;
    float glossiness;
    Color absorption;
    float ior;  // index of refraction
    unsigned int viewportTextureID;
};
 
//-------------------------------------------------------------------------------
 
#endif