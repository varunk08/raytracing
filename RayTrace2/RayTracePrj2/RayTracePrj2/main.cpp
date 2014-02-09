#include "tinyxml/tinyxml.h"
#include <GLUT/GLUT.h>
#include "viewport.cpp"
#include "scene.h"
#include "xmlload.cpp"
#include <iostream>
#include <math.h>
#include <cmath>
#include "cyColor.h"

using namespace std;

#define _USE_MATH_DEFINES

Camera camera;
Node rootNode;
RenderImage renderImage;
Sphere theSphere;
MaterialList materials;
LightList lights;

Color24 white(cyColor({255,255,255}));
Color24 black(cyColor({0,0,0}));

bool RayTrace(HitInfo &hitInfo, Node* node,Ray ray,int PixIndex);


void BeginRender()
{
   
    cout<<"\nBeginning Render...";
    
    float alpha = camera.fov;
    float l = 1.0;
    float h = l * tan(alpha/2.0 *(M_PI/180));
    float aspectRatio = (float)camera.imgWidth/camera.imgHeight;
    float s = aspectRatio * abs(h);
    float dx = (2 * abs(s))/camera.imgWidth;
    float dy = -(2 * abs(h))/camera.imgHeight;
    float dxx = dx/2,dyy=dy/2;
    Point3 K(-s,h,-l);
    
    K.x += dxx;
    K.y += dyy;
   // cout<<"K: "<<K.x<<" "<<K.y<<" "<<K.z<<endl;
    for(int i = 0; i< camera.imgHeight; i++){
        for(int j = 0; j<camera.imgWidth; j++){
            int PixIndex = i * camera.imgWidth + j;
            bool pixelHit=false;
            K.x += dx;
            Matrix3 RotMat;
            cyPoint3f f = camera.dir;
            f.Normalize();
            cyPoint3f s = f.Cross(camera.up);
            s.Normalize();
            cyPoint3f u = s.Cross(f);
            const float pts[9]={s.x,u.x,-f.x,s.y,u.y,-f.y,s.z,u.z,-f.z};
            RotMat.Set(pts);
           // RotMat.SetView(camera.dir, camera.up);
            Ray r(camera.pos, K);
            r.dir=r.dir*RotMat;
            r.dir.Normalize();
            //cout<<"RAY: "<<r.p.x<<" "<<r.p.y<<" "<<r.p.z<<" "<<r.dir.x<<" "<<r.dir.y<<" "<<r.dir.z<<endl;
            HitInfo hitInfo;
            hitInfo.Init(); 

            Color shade(255,255,255);
            if(rootNode.GetNumChild()>0){
                for(int k=0; k < rootNode.GetNumChild(); ++k){
                    Node* node = rootNode.GetChild(k);
                    
                    if(RayTrace(hitInfo, node,r,PixIndex)){
                        pixelHit=true;
                        if(hitInfo.node){
                            if(hitInfo.node->GetMaterial()){
                                //hitInfo.node->FromNodeCoords(hitInfo);
                                shade = hitInfo.node->GetMaterial()->Shade(r, hitInfo, lights);
                            }
                    
                        renderImage.PutPixel(PixIndex, shade, hitInfo.z);
                        }
                        
                    }
                    
                }
            } 
            
            if(!pixelHit){renderImage.PutPixel(PixIndex, black, BIGFLOAT);}
        }
        K.x = -s;
        K.x += dxx;
        K.y += dy;
    }
    cout<<"Render Complete"<<endl;
    renderImage.ComputeZBufferImage();
    renderImage.SaveZImage("/Users/varunk/Desktop/RayTracerProj1/RayTrace2/RayTracePrj2/RayTracePrj2/zbuffer.ppm");
    renderImage.SaveImage("/Users/varunk/Desktop/RayTracerProj1/RayTrace2/RayTracePrj2/RayTracePrj2/renderimage.ppm");
    
}
void MtlBlinn::SetViewportMaterial() const{
    ColorA c;
    c = diffuse;
    glMaterialfv( GL_FRONT, GL_AMBIENT_AND_DIFFUSE, &c.r );
    c = specular;
    glMaterialfv( GL_FRONT, GL_SPECULAR, &c.r );
    glMaterialf( GL_FRONT, GL_SHININESS, glossiness );
}
Color MtlBlinn::Shade(const Ray &ray, const HitInfo &hInfo, const LightList &lights) const{
    //cout<<"shading..";
    Color shade;
    const Material *mat;
    mat = hInfo.node->GetMaterial();
    const MtlBlinn* mb =static_cast<const MtlBlinn*>(mat);
    Color ambInt = mb->diffuse;
    Color allOther(0,0,0);
    Color diffuse = mb->diffuse;;
    Color ambComponent(0,0,0);
    
    for ( unsigned int i=0; i<lights.size(); i++ ) {
        
        if(lights[i]->IsAmbient()){
            //cout<<"ambient";
            Color intensity = lights[i]->Illuminate(hInfo.p);
            ambComponent += (ambInt * intensity);
           // cout<<"ambComponene: "<<ambComponent.r<<" "<<ambComponent.g<<endl;
            continue;
        }
        else{
            
            Point3 L = -lights[i]->Direction(hInfo.p);
            //cout<<"L: "<<L.x<<" "<<L.y<<" "<<L.z<<endl;
            L.Normalize();
            Point3 V = camera.pos - hInfo.p;
            V.Normalize();
            
            Point3 LplusV = L + V;
            
            Point3 H = (L+V)/LplusV.Length();
            H.Normalize();
            
            float alpha = mb->glossiness;
            Point3 N = hInfo.N;
            //N.Normalize();
            float S = H.Dot(N);
            S = pow(S,alpha);
            
            float costheta = L.Dot(N)/(L.Length() * N.Length());
            
            Color intensity = lights[i]->Illuminate(hInfo.p);
            
            allOther += intensity * (costheta>0?costheta:0) * (diffuse + S * (mb->specular)) ;
            
            
            
        }
        // finally add inta*cola + intall*costheta*(cold + s* colS)
        shade = ambComponent  + allOther;
    }
    return shade;
};


bool RayTrace(HitInfo &hitInfo, Node* curnode, Ray ray, int PixIndex)
{
    Node* node = curnode;
    bool hitTest = false;
    const Object *obj = node->GetObject();
    
    ray = node->ToNodeCoords(ray);
    HitInfo tempHitInfo;
    tempHitInfo.Init();
    tempHitInfo.node = node;
    hitTest = obj->IntersectRay(ray, tempHitInfo);
    node->FromNodeCoords(tempHitInfo);
    if(hitTest && tempHitInfo.z < hitInfo.z){
        hitInfo = tempHitInfo;
        //cout<<hitInfo.z<<endl;
    }
    else hitTest=false;
        if(node->GetNumChild()>0)
        {
            
            for(int i=0;i<node->GetNumChild();i++)
            {
                node = node->GetChild(i);
                
                if(RayTrace(hitInfo, node, ray, PixIndex)){
                    curnode->FromNodeCoords(hitInfo);
                    hitTest = true;
                }
            }
        }
   
    if(hitTest) return true;
    else return false;
    
}

void StopRender()
{
    
}
int main(int argc, char* argv[])
{
    const char* filename = "/Users/varunk/Desktop/RayTracerProj1/RayTrace2/RayTracePrj2/RayTracePrj2/scene.xml";
    LoadScene(filename);
    
    glutInit(&argc,argv);
    ShowViewport();
    
    
	return 0;
}
