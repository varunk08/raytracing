
//-------------------------------------------------------------------------------
///
/// \file       materials.h
/// \author     Cem Yuksel (www.cemyuksel.com)
/// \version    10.0
/// \date       November 4, 2013
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
    virtual Color Shade(const Cone &ray, const HitInfo &hInfo, const LightList &lights, int bounceCount, bool globalAmbient=false) const;
    
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
    void SetReflectionGlossiness(float gloss) { reflectionGlossiness=gloss; }
    void SetRefractionGlossiness(float gloss) { refractionGlossiness=gloss; }
    
    virtual void SetViewportMaterial() const;   // used for OpenGL display
    // Photon Extensions
    virtual bool IsPhotonSurface() const { return diffuse.GetColor().Grey() > 0; }   // if this method returns true, the photon will be stored
    virtual bool RandomPhotonBounce(Ray &r, Color &c, const HitInfo &hInfo) const;  // if this method returns true, a new photon with the given direction and color will be traced
    
private:
    TexturedColor diffuse, specular, reflection, refraction;
    float glossiness;
    Color absorption;
    float ior;  // index of refraction
    float reflectionGlossiness, refractionGlossiness;
    unsigned int viewportTextureID;
};



//-------------------------------------------------------------------------------

#endif