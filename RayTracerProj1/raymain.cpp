#include "tinyxml/tinyxml.h"
#include <GLUT/GLUT.h>
#include "viewport.cpp"
#include "scene.h"
#include "xmlload.cpp"
#include <iostream>
#include <math.h>
#include <cmath>

using namespace std;
#define _USE_MATH_DEFINES
Camera camera;
Node rootNode;
RenderImage renderImage;
Color24 white = {255,255,255};
Color24 black = {0,0,0};
void RayTrace(Node* node,Ray ray,int PixIndex);
bool TraceNode(const Ray &r, HitInfo &hInfo,const Node &node);

bool RayTrace_2(const Ray &ray, HitInfo &hitInfo)
{
    return TraceNode(ray, hitInfo, rootNode);
    
}

bool TraceNode(const Ray &r, HitInfo &hInfo,const Node &node)
{
    Ray ray;
    ray = node.ToNodeCoords(r);
    const Object* obj = node.GetObject();
    bool objHitTest = false;
    
    if (node.GetNumChild() > 0)
    {
        
        for (int i = 0; i < node.GetNumChild(); i++)
        {
            bool childHit = false;
            const Node &childNode = *node.GetChild(i);
            childHit = TraceNode(ray, hInfo, childNode);
            if(childHit)
            {
                objHitTest = true;
            }
        }
    }
    if(obj)
    {
        if(obj->IntersectRay(ray, hInfo)){
            objHitTest=true;
            hInfo.node = &node;
        }
        if(objHitTest){
            //hInfo.node = &node;
            //node.FromNodeCoords(hInfo);
        }
    }
    
    return objHitTest;
}

void BeginRender()
{

    cout<<"\nBeginning Render...";
    
    float alpha = camera.fov;
    float l = 1.0;
    float h = l * tan(alpha/2.0 *(M_PI/180));
    float aspectRatio = (float)camera.imgWidth/camera.imgHeight;
    float s = aspectRatio * h;
    float dx = (2 * s)/camera.imgWidth;
    float dy = -(2 * h)/camera.imgHeight;
    float dxx = dx/2,dyy=dy/2;
    Point3 K(-s,h,-l);
    K.x += dxx;
    K.y += dyy;
    for(int i = 0; i< camera.imgHeight; i++){
        
        for(int j = 0; j<camera.imgWidth; j++){
            
            K.x += dx;
            Matrix3 RotMat;
            
            Point3 dvec = camera.dir - camera.pos;
            Point3 svec = camera.up.Cross(dvec);
            dvec.Normalize();
            svec.Normalize();
            camera.up.Normalize();
            RotMat.Set(svec,camera.up, dvec);
            Ray r(camera.pos, K);
            
            r.dir=r.dir*RotMat;
            
            r.dir.Normalize();
            
            HitInfo hInfo;
            hInfo.Init();
            if(rootNode.GetNumChild()>0){
//                for(int k=0; k < rootNode.GetNumChild(); ++k){
//                    RayTrace(rootNode.GetChild(k),r,i * camera.imgWidth + j);
//                }
                if(RayTrace_2(r, hInfo))
                {
                    renderImage.PutPixel(i *camera.imgWidth+j, white, hInfo.z);
                }
                else renderImage.PutPixel(i *camera.imgWidth+j, black, BIGFLOAT);
            }
            
            
        }
        K.x = -s;
        K.x += dxx;
        K.y += dy;
    }
    cout<<"Render Complete"<<endl;
    renderImage.ComputeZBufferImage();
    renderImage.SaveZImage("/Users/varunk/Desktop/RayTracerProj1/RayTracerProj1/zbuffer.ppm");
    renderImage.SaveImage("/Users/varunk/Desktop/RayTracerProj1/RayTracerProj1/renderimage.ppm");
    
}




void RayTrace(Node* node, Ray ray, int PixIndex)
{
    HitInfo hitInfo;
    hitInfo.node=node;
    hitInfo.z = 0;
    bool hitTest = false;
    Object *obj = node->GetObject();
    
    ray = node->ToNodeCoords(ray);
    hitTest = obj->IntersectRay(ray, hitInfo);
    //cout<<"Z: "<<hitInfo.z<<endl;
   
    if(hitTest){
        
        renderImage.PutPixel(PixIndex, white, hitInfo.z);
        return;
        
    }
    
    else{
        
        if(node->GetNumChild()>0)
        {
            
            for(int i=0;i<node->GetNumChild();i++)
            {
                RayTrace(node->GetChild(i),ray,PixIndex);
            }
        }
        
    }
    
}

void StopRender()
{
    
}
int main(int argc, char* argv[])
{
    const char* filename = "/Users/varunk/Desktop/RayTracerProj1/RayTracerProj1/originalscene.xml";
    LoadScene(filename);
    
    glutInit(&argc,argv);
    ShowViewport();
    
    
	return 0;
}
