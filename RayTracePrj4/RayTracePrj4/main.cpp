
#include <GLUT/GLUT.h>
#include <iostream>
#include <math.h>
#include <cmath>
#include <pthread.h>

#include "viewport.cpp"
#include "scene.h"
#include "xmlload.cpp"
#include "cyColor.h"
#include "tinyxml/tinyxml.h"

using namespace std;

#define _USE_MATH_DEFINES
struct RenderParams{
    int startX;
    int startY;
    int endX;
    int endY;
};
Camera camera;
Node rootNode;
RenderImage renderImage;
Sphere theSphere;
MaterialList materials;
LightList lights;
const int NUM_THREADS = 2;

pthread_t threadID[NUM_THREADS];
RenderParams args[NUM_THREADS];

Color24 white(cyColor({255,255,255}));
Color24 black(cyColor({0,0,0}));

bool TraceNode(const Ray &r,HitInfo &hitInfo,const Node &curNode);
bool RayTrace(HitInfo &hitInfo, Node* node,Ray ray,int PixIndex);
bool RayTrace_2(const Ray &ray, HitInfo &hitInfo);

void *doRender(void* arg){
    
    RenderParams rarg = *((RenderParams *)arg);
    cout<<"\nBeginning Render...";
     cout<<"startX: "<<rarg.startX<<endl;
         cout<<"startY: "<<rarg.startY<<endl;
         cout<<"endX: "<<rarg.endX<<endl;
         cout<<"endY: "<<rarg.endY<<endl;
    float alpha = camera.fov;
    float l = 1.0;
    float h = l * tan(alpha/2.0 *(M_PI/180));
    
    float aspectRatio = (float)camera.imgWidth/camera.imgHeight;
    float s = aspectRatio * abs(h);
//    cout<<"S: "<<s<<endl;
    float dx = (2 * abs(s))/camera.imgWidth;
    float dy = -(2 * abs(h))/camera.imgHeight;
    float dxx = dx/2,dyy=dy/2;
    if(rarg.startX == 400){
        s = 0;
    }
    Point3 K(-s,h,-l);
//    cout<<"K.x: "<<K.x<<"K.y: "<<K.y<<"K.z: "<<K.z<<endl;
    K.x += (dxx );
    K.y += (dyy );
    
    for(int i = rarg.startY; i<= rarg.endY ; i++){
        for(int j = rarg.startX; j<= rarg.endX; j++){
            int PixIndex = i * camera.imgWidth + j;
            bool pixelHit=false;
            //K.x += rarg.startX * dx;
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
            
            //   cout<<rootNode.GetNumChild()<<endl;
            Color shade(255,255,255);
            if(rootNode.GetNumChild()>0){
                //                for(int k=0; k < rootNode.GetNumChild(); ++k){
                //                    Node* node = rootNode.GetChild(k);
                ////                    cout<<"K: "<<k<<endl;
                //                    if(RayTrace(hitInfo, node,r,PixIndex)){
                //
                //                        pixelHit=true;
                //                        if(hitInfo.node){
                //                            if(hitInfo.node->GetMaterial()){
                //
                ////                                cout<<"Shading "<<PixIndex<<" "<<hitInfo.z<< endl;
                //                                shade = hitInfo.node->GetMaterial()->Shade(r, hitInfo, lights, 1);
                //                            }
                //
                //                            renderImage.PutPixel(PixIndex, shade, hitInfo.z);
                //                        }
                //
                //                    }
                ////                    cout<<"not hit it"<<endl;
                //
                //                }
                if(RayTrace_2(r, hitInfo)) {
                    pixelHit=true;
                    shade = hitInfo.node->GetMaterial()->Shade(r, hitInfo, lights, 16);
                }
                
                renderImage.PutPixel(PixIndex, shade, hitInfo.z);
                
            }
            
            if(!pixelHit){renderImage.PutPixel(PixIndex, black, BIGFLOAT);}
        }
        K.x = -s;
        K.x += dxx;
        K.y += dy;
    }
    cout<<"Render Complete"<<endl;
    renderImage.ComputeZBufferImage();
    renderImage.SaveZImage("/Users/varunk/Desktop/RayTracerProj1/RayTracePrj4/RayTracePrj4/zbuffer.ppm");
    renderImage.SaveImage("/Users/varunk/Desktop/RayTracerProj1/RayTracePrj4/RayTracePrj4/renderimage.ppm");
    
    return NULL;
}

void BeginRender()
{
    

    int threadError;
    args[0].startX = 0;
    args[0].startY = 0;
    args[0].endX = 399;
    args[0].endY = 599;
    
    args[1].startX = 400;
    args[1].startY = 0;
    args[1].endX = 799;
    args[1].endY = 599;
    
    for(int i=0; i<NUM_THREADS; i++){
       /* args[i].startX = i * camera.imgWidth/2;
        if(args[i].startX >=camera.imgWidth) args[i].startX -= camera.imgWidth;
        if(i < NUM_THREADS/2) args[i].startY = 0;
        else args[i].startY = camera.imgHeight/2;
        
        args[i].endX = args[i].startX + camera.imgWidth/2 - 1;
        args[i].endY = args[i].startY + camera.imgHeight/2 - 1; */
        //cout<<args[i].startX<<","<<args[i].startY<<endl;
        //cout<<args[i].endX<<","<<args[i].endY<<endl;
        threadError = pthread_create(&threadID[i],NULL,&doRender,&args[i]);
        if(threadError != 0){
            cout<<"THREAD ERROR..."<<i<<endl;
        }
    }
    


    
    
    
    
}

float GenLight::Shadow(Ray ray, float t_max)
{
   // cout<<"Calculating shadow"<<endl;
    //add bias
    float eps = 1e-3f;
    ray.p = Point3(ray.p.x+eps*ray.dir.x, ray.p.y+eps*ray.dir.y,ray.p.z+eps*ray.dir.z);
    //cout<<"RAY: "<<ray.p.x<<" "<<ray.p.y<<" "<<ray.p.z<<" "<<endl;
    //ray trace - if hit return shadow - else intensity
    int PixIndex=0;
    //ray.dir.Normalize();
    HitInfo hitInfo;
    
    //cout<<rootNode.GetNumChild()<<endl;
//    for(int i=0; i < rootNode.GetNumChild(); i++){
//        Node* curnode = rootNode.GetChild(i);
//        
//        hitInfo.Init();
//        hitInfo.node = curnode;
//       //cout<<"Ray tracing "<<i<<endl;
//        if(RayTrace( hitInfo, curnode,ray, PixIndex))
//        {
//            //cout<<"HIT "<<hitInfo.z<<endl;
//
//            if(hitInfo.z>0 && hitInfo.z < t_max){
//               
//                return 0.0; //occluded
//            }
//            
//        }
//    }
    if(RayTrace_2(ray, hitInfo))
    {
        //cout<<"HIT "<<hitInfo.z<<endl;
        
        if(hitInfo.z>0 && hitInfo.z < t_max){
            
            return 0.0; //occluded
        }
        
    }
    return 1.0; //direct
}

Color MtlBlinn::Shade(const Ray &ray, const HitInfo &hInfo, const LightList &lights, int bounceCount) const{
    float bias = 0.001f;
    Color shade;
    Color rShade = Color(0,0,0);
    Color tShade = Color(0,0,0);
    const Material *mat;
    mat = hInfo.node->GetMaterial();
    const MtlBlinn* mb =static_cast<const MtlBlinn*>(mat);
//    cout<<"HInfo front: "<<hInfo.front<<endl;
    /* local copy */
    Point3 P = hInfo.p;
    Ray iRay = ray;
    
    Color ambInt = mb->diffuse;
    Color allOther = Color(0,0,0);
    Color diffuse = mb->diffuse;;
    Color ambComponent = Color(0,0,0);
    
    for ( unsigned int i=0; i<lights.size(); i++ ) {
        if(lights[i]->IsAmbient()){
//            cout<<"ambient "<<endl;
            Color intensity = lights[i]->Illuminate(hInfo.p);
            ambComponent += (ambInt * intensity);
            continue;
        }
        else{
//            cout<<"other lighting  "<<endl;
            Point3 L = -lights[i]->Direction(P);
            L.Normalize();
            
            Point3 V = ray.p - P;
            V.Normalize();
            
            Point3 LplusV = L + V;
            Point3 H = (L+V)/LplusV.Length();
            H.Normalize();
            
            float alpha = mb->glossiness;
            Point3 N = hInfo.N;
            float S = H.Dot(N);
            S = pow(S,alpha);
            float costheta = L.Dot(N)/(L.Length() * N.Length());
            Color intensity = lights[i]->Illuminate(P);
//            cout<<"costheta "<<endl;
            allOther += intensity * (costheta>0?costheta:0) * (diffuse + S * (mb->specular)) ;
        }
        /* finally add inta*cola + intall*costheta*(cold + s* colS)*/
        shade = ambComponent  + allOther;
    }
    
    /* Calculate refraction */
    if(refraction.Grey()>0 && bounceCount>0){
        Color reflShade = Color(0,0,0);
        float R0, Refl = 0.0f, Trans = 0.0f;
        HitInfo temp;
        temp.Init();
        
        Point3 N = hInfo.N;
//        Point3 V = Point3(iRay.p.x -  hInfo.p.x, iRay.p.y - hInfo.p.y, iRay.p.z - hInfo.p.z);
        Point3 V = Point3(hInfo.p.x - iRay.p.x, hInfo.p.y - iRay.p.y, hInfo.p.z - iRay.p.z);
        V.Normalize();
        float n1 = 1, n2 = 1;
        if(hInfo.front){ /* Hitting from outside */
            //temp.front = false;
            n2 = ior;
//            cout<<"outside "<<endl;
        }
        else if(!hInfo.front){ /* Transmission from the inside */
            //temp.front = true;
            n1 = ior;
//            cout<<"intside... "<<endl;
            N = -hInfo.N;
        }
        float ratio_n = n1 / n2;
        
        float costheta_v = -V.Dot(N);        /* refer: http://graphics.stanford.edu/courses/cs148-10-summer/docs/2006--degreve--reflection_refraction.pdf */

        float sin2theta_t = ratio_n * ratio_n * (1 - costheta_v * costheta_v);
        Point3 T =   ratio_n * V + (ratio_n * costheta_v - sqrtf(1 - sin2theta_t)) * N ;
//        cout<<ratio_n<<" "<<"cos_v "<<costheta_v<<" sin2theta_t "<<sin2theta_t<<endl;
        Ray tRay = Ray(hInfo.p,T);
        
        //tRay.dir.Normalize();
        tRay.p.x = tRay.p.x + bias *tRay.dir.x; /* add bias */
        tRay.p.y = tRay.p.y + bias *tRay.dir.y;
        tRay.p.z = tRay.p.z + bias *tRay.dir.z;
//        cout<<"B temp front: "<< temp.front<<endl;
        if(sin2theta_t <= 1){
            if(RayTrace_2(tRay, temp)){
//                bounceCount--;
//                cout<<"A temp front: "<< temp.front<<endl;
                tShade  =  temp.node->GetMaterial()->Shade(tRay,temp,lights,bounceCount);
                tShade.r *= exp(-absorption.r * temp.z);
                tShade.g *= exp(-absorption.g * temp.z);
                tShade.b *= exp(-absorption.b * temp.z);
//                shade = tShade; /* remove later */
//                return shade;
               
                
                /* Calculate Schlick's approximation */
                
                R0 = (n1 - n2)/(n1 + n2);
                R0 *= R0;
                double  X = 0.0;
//                if(n1 > n2){
//                    X = 1.0 - sqrtf(1.0 - sin2theta_t);
//                }
//                else{ X = 1.0 - costheta_v; }
                X = 1.0 - costheta_v;
                Refl = R0 + (1.0 - R0) *  X * X * X * X * X;
                Trans = 1.0 - Refl;
                
            }
        }
        else {/* Total internal reflection */
            Refl = 1.0f;
        }
        
        /* Calculate reflection due to reflectance */
        if(bounceCount >0){
            N = hInfo.N;
            Point3 V = Point3(iRay.p.x -  P.x, iRay.p.y - P.y, iRay.p.z - P.z);
            //V.Normalize();
            Point3 VR = 2 * V.Dot(N) * N - V;
            //VR.Normalize();
            Ray rRay = Ray(P, VR);
            //rRay.dir.Normalize();
            rRay.p.x = rRay.p.x + bias *rRay.dir.x;
            rRay.p.y = rRay.p.y + bias *rRay.dir.y;
            rRay.p.z = rRay.p.z + bias *rRay.dir.z;
            HitInfo temp1;
            temp1.Init();
            if(rootNode.GetNumChild()>0){
                if(RayTrace_2(rRay, temp1)){
                    bounceCount --;
                    reflShade =   temp1.node->GetMaterial()->Shade(rRay, temp1, lights, bounceCount);
                }
            }
        }
        
//        cout<<"Refl: "<<Refl<<"Trans "<<Trans<<endl;
        tShade = refraction * (Trans * tShade + Refl * reflShade);
        
        
    }
    





    /* calculate reflection*/
    if(reflection.Grey()>0 && bounceCount > 0){

        Point3 N = hInfo.N;
        Point3 V = Point3(iRay.p.x -  P.x, iRay.p.y - P.y, iRay.p.z - P.z);
       // V.Normalize();
        Point3 VR = 2 * V.Dot(N) * N - V;
        Ray rRay = Ray(hInfo.p, VR);
        //rRay.dir.Normalize();
        rRay.p.x = rRay.p.x + bias *rRay.dir.x;
        rRay.p.y = rRay.p.y + bias *rRay.dir.y;
        rRay.p.z = rRay.p.z + bias *rRay.dir.z;
        HitInfo temp;
        temp.Init();
        if(rootNode.GetNumChild()>0){
            if(RayTrace_2(rRay, temp)){
                bounceCount--;
                rShade = reflection * temp.node->GetMaterial()->Shade(rRay, temp, lights, bounceCount);
            }
        }
    }
    
  
    
    /* Add shade with reflected and refracted colors */
    shade += (rShade + tShade);
    return shade;
};


bool RayTrace(HitInfo &hitInfo, Node* curnode, Ray ray, int PixIndex)
{
    Node* node = curnode;
    bool hitTest = false;
    
    const Object *obj = node->GetObject();
    ray = curnode->ToNodeCoords(ray);
    if(obj){
//       cout<<"Transforming to..."<<endl;
        
        HitInfo tempHitInfo;
        tempHitInfo.Init();
        tempHitInfo.node = node;
        tempHitInfo.z = hitInfo.z;
        hitTest = obj->IntersectRay(ray, tempHitInfo);
        node->FromNodeCoords(tempHitInfo);
        if(hitTest && tempHitInfo.z < hitInfo.z){
            hitInfo = tempHitInfo;
            //cout<<hitInfo.z<<endl;
        }
        //else hitTest=false;
    }
    if(node->GetNumChild()>0)
    {
        //cout<<"Children "<<node->GetNumChild()<<endl;
        for(int i=0;i<curnode->GetNumChild();++i)
        {
//            cout<<"Child "<<i<<endl;
            node = curnode->GetChild(i);
            HitInfo temp;
            temp.Init();
            temp = hitInfo;
            if(RayTrace(hitInfo, node, ray, PixIndex)){
                curnode->FromNodeCoords(hitInfo);
                
//                cout<<"Transforming from "<<curnode->GetNumChild()<<endl;
                if(temp.z > hitInfo.z) hitTest = true;
                else{
                   // hitInfo = temp;
                    hitTest = false;
                    continue;
                }
            }
            
        }
    }
    
    
    if(hitTest) return true;
    else return false;
    
}

bool RayTrace_2(const Ray &ray, HitInfo &hitInfo)
{
//    cout<<"Ray trace "<<endl;
    return TraceNode(ray, hitInfo, rootNode);
  
}

bool TraceNode(const Ray &r, HitInfo &hInfo,const Node &node)
{
//    cout<<"Tracenode"<<endl;
    Ray ray = r;
    ray = node.ToNodeCoords(ray);
//    cout<<"...To node coords"<<endl;
    const Object* obj = node.GetObject();
    bool objHitTest = false;
    bool childHit = false;
//    int hitSide = HIT_FRONT;
    if(obj)
    {
//        if(!hInfo.front){ hitSide = HIT_FRONT_AND_BACK; } /* Back face hit for refraction */
        if(obj->IntersectRay(ray, hInfo)){
            objHitTest=true;
            hInfo.node = &node;
            //node.FromNodeCoords(hInfo);
//            cout<<" - Hit -"<<endl;
//            cout<<"..From node coords"<<endl;
        }
    }
    if (node.GetNumChild() > 0)
    {

        for (int i = 0; i < node.GetNumChild(); i++)
        {
            
//            cout<<"Child "<<i<<endl;
            
            const Node &childNode = *node.GetChild(i);
            if(TraceNode(ray, hInfo, childNode)){
                    childHit = true;
            }
            
//            cout<<"Child "<<i<<" out"<<endl;
        }
        if(childHit)
        {
            //node.FromNodeCoords(hInfo);
            objHitTest = true;
        }
    }
    
    
    if(objHitTest){
        //hInfo.node = &node;
//        cout<<"..From node coords - child hit not parent"<<endl;
       node.FromNodeCoords(hInfo);
    }
    return objHitTest;
}

void StopRender()
{
    
}

void MtlBlinn::SetViewportMaterial() const{
    ColorA c;
    c = diffuse;
    glMaterialfv( GL_FRONT, GL_AMBIENT_AND_DIFFUSE, &c.r );
    c = specular;
    glMaterialfv( GL_FRONT, GL_SPECULAR, &c.r );
    glMaterialf( GL_FRONT, GL_SHININESS, glossiness );
}

int main(int argc, char* argv[])
{
    const char* filename = "/Users/varunk/Desktop/RayTracerProj1/RayTracePrj4/RayTracePrj4/scene.xml";
    LoadScene(filename);
    
    glutInit(&argc,argv);
    ShowViewport();
    
    
	return 0;
}
