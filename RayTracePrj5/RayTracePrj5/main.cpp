#include "tinyxml/tinyxml.h"
#include <GLUT/GLUT.h>
#include "viewport.cpp"
#include "scene.h"
#include "xmlload.cpp"
#include <iostream>
#include <math.h>
#include <cmath>
#include "cyColor.h"
#include <pthread.h>
#include <unistd.h>
#include <vector>
#include "threadpool.h"
using namespace std;
#define BIAS_SHADOW 1e-4f
#define BIAS_SHADING 0.001f
#define _USE_MATH_DEFINES
#define NUM_THREADS  2
struct RenderParams{
    Point3 K;
    int pixIndex;
    Point2 pixLocation;
    bool renderComplete = false;
    Ray ray;
};


struct ImageParams{
    unsigned int NUM_PIXELS;
    vector<int> PixIndex;
    vector<bool> rendered;
    vector<Point2> PixLocation;
    vector<Point3> K;
    vector<Ray> Ray;
}imageParams;
ThreadPool tp(NUM_THREADS);
Camera camera;
Node rootNode;
RenderImage renderImage;
Sphere theSphere;
Plane thePlane;
MaterialList materials;
LightList lights;
ObjFileList objList;


pthread_t threadID[NUM_THREADS];
pthread_mutex_t getPix_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t setPix_mutex = PTHREAD_MUTEX_INITIALIZER;
RenderParams args[NUM_THREADS];

Color24 white(cyColor({255,255,255}));
Color24 black(cyColor({0,0,0}));
RenderParams giveMeAPixelToRender();
bool TraceNode(const Ray &r,HitInfo &hitInfo,const Node &curNode);
bool RayTrace(HitInfo &hitInfo, Node* node,Ray ray,int PixIndex);
bool RayTrace_2(const Ray &ray, HitInfo &hitInfo);
void PopulateImageParams();
void doRender(void* arg);
void SetPixelAsRendered(int pixelIndex);

void PopulateImageParams()
{
    cout<<"Populating ImageParams..."<<endl;
    float alpha = camera.fov;
    float l = 1.0;
    float h = l * tan(alpha/2.0 *(M_PI/180));
    
    float aspectRatio = (float)camera.imgWidth/camera.imgHeight;
    float s = aspectRatio * abs(h);
    float dx = (2 * abs(s))/camera.imgWidth;
    float dy = -(2 * abs(h))/camera.imgHeight;
    float dxx = dx/2,dyy=dy/2;
    Point3 K(-s,h,-l);
    K.x += (dxx );
    K.y += (dyy );
    
    for(int i = 0; i< camera.imgHeight ; i++){
        for(int j = 0; j< camera.imgWidth; j++){
            K.x += dx;
            Matrix3 RotMat;
            cyPoint3f f = camera.dir;
            f.Normalize();
            cyPoint3f s = f.Cross(camera.up);
            s.Normalize();
            cyPoint3f u = s.Cross(f);
            const float pts[9]={s.x,u.x,-f.x,s.y,u.y,-f.y,s.z,u.z,-f.z};
            RotMat.Set(pts);
            Ray r(camera.pos, K);
            r.dir=r.dir*RotMat;
            r.dir.Normalize();
            
            /* Populating the Struct */
            Point2 pixLoc = Point2(j,i);
            //imageParams.K.push_back(K);
            imageParams.rendered.push_back(false);
            imageParams.PixLocation.push_back(pixLoc);
            imageParams.PixIndex.push_back( i * camera.imgWidth + j);
            imageParams.Ray.push_back(r);
                        
        }
        K.x = -s;
        K.x += dxx;
        K.y += dy;
    }
    
   

}
void BeginRender()
{
    //Load ImageParams
    imageParams.NUM_PIXELS = camera.imgHeight * camera.imgWidth;
    PopulateImageParams();
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  /*  for(int i = 0; i< camera.imgHeight ; i++){
        for(int j = 0; j< camera.imgWidth; j++){
            HitInfo hitInfo;
            hitInfo.Init();
            int PixIndex = i * camera.imgWidth + j;
            Point2 pixLoc = imageParams.PixLocation.at(PixIndex);
            Ray r = imageParams.Ray.at(PixIndex);
            bool pixelHit = false;
            Color shade(255,255,255);

            if(rootNode.GetNumChild()>0){
                
                if(RayTrace_2(r, hitInfo)) {
                    pixelHit=true;
                    shade = hitInfo.node->GetMaterial()->Shade(r, hitInfo, lights, 5);
                }
                
                renderImage.PutPixel(PixIndex, shade, hitInfo.z);
                
                
            }
            
            if(!pixelHit){
                
                renderImage.PutPixel(PixIndex, black, BIGFLOAT);
                
            }

            
        }
    }
   
    int threadError;
    cout<<"Creating threads"<<endl;
    for(int i=0; i<NUM_THREADS; i++){
        args[i] = giveMeAPixelToRender();
        if(args[i].renderComplete != 1){
            threadError = pthread_create(&threadID[i],&attr,&doRender,&args[i]);
            cout<<"Thread ID: "<<threadID[i]<<endl;
            if(threadError != 0){
                cout<<"THREAD ERROR..."<<i<<endl;
            }
            
        }
    }
    for(int i =0;i<NUM_THREADS;i++){
        pthread_join(threadID[i], NULL);
    }
    
*/
    
    int ret = tp.initialize_threadpool();
    if (ret == -1) {
        cerr << "Failed to initialize thread pool!" << endl;
        //return 0;
    }
    else{
    for (int i = 0; i < imageParams.NUM_PIXELS; i++) {
        RenderParams* x = new RenderParams();
        x->pixIndex = imageParams.PixIndex.at(i);
        x->pixLocation = imageParams.PixLocation.at(i);
        x->ray = imageParams.Ray.at(i);

        Task* t = new Task(&doRender, (void*) x);
        //    cout << "Adding to pool, task " << i+1 << endl;
        tp.add_task(t);
//            cout << "Added to pool, task " << i+1 << endl;
    }
    
    //usleep(20);
    
    
    }
    cout<<"Render Complete"<<endl;
    renderImage.ComputeZBufferImage();
    renderImage.SaveZImage("/Users/varunk/Desktop/RayTracerProj1/RayTracePrj5/RayTracePrj5/zbuffer.ppm");
    renderImage.SaveImage("/Users/varunk/Desktop/RayTracerProj1/RayTracePrj5/RayTracePrj5/renderimage.ppm");
    
    
}
RenderParams giveMeAPixelToRender()
{
    //cout<<"Getting Pixel"<<endl;
    pthread_mutex_lock(&getPix_mutex);
    RenderParams rParams;
    int i =0;
    for ( i = 0; i<imageParams.NUM_PIXELS; i++) {
        if(imageParams.rendered[i] != 1){
            //rParams.K = imageParams.K.at(i);
            rParams.pixIndex = imageParams.PixIndex.at(i);
            rParams.pixLocation = imageParams.PixLocation.at(i);
            rParams.ray = imageParams.Ray.at(i);
            break;
        }
        
    }
    if(i == imageParams.NUM_PIXELS){
        rParams.renderComplete = true;

    }
    pthread_mutex_unlock(&getPix_mutex);
    return rParams;
}

void doRender(void* arg){
    
    RenderParams rarg = *((RenderParams *)arg);
    //cout<<"Do render...."<<endl;
            bool pixelHit=false;

    
            HitInfo hitInfo;
            hitInfo.Init();
            
            Point2 pixLoc = rarg.pixLocation;
            Ray r = rarg.ray;
            int PixIndex = rarg.pixIndex;
            Color shade(255,255,255);
            if(rootNode.GetNumChild()>0){
                
                if(RayTrace_2(r, hitInfo)) {
                    pixelHit=true;
//                    cout<<"Shading...."<<endl;
                    shade = hitInfo.node->GetMaterial()->Shade(r, hitInfo, lights, 5);
                }
                
                renderImage.PutPixel(PixIndex, shade, hitInfo.z);
                
                
            }
            
            if(!pixelHit){
                
                renderImage.PutPixel(PixIndex, black, BIGFLOAT);
            
            }
    
   
    
//    RenderParams renderArg = giveMeAPixelToRender();
//    if(renderArg.renderComplete != 1){
//            doRender(&renderArg);
//    }
//    
    

}

void SetPixelAsRendered(int pixIndex)
{
    pthread_mutex_lock(&setPix_mutex);
     imageParams.rendered.at(pixIndex) = true;
    pthread_mutex_unlock(&setPix_mutex);
}
float GenLight::Shadow(Ray ray, float t_max)
{
   // cout<<"Calculating shadow"<<endl;
    //add bias
    float eps = BIAS_SHADOW;
    ray.p = Point3(ray.p.x+eps*ray.dir.x, ray.p.y+eps*ray.dir.y,ray.p.z+eps*ray.dir.z);
    HitInfo hitInfo;
    
    if(RayTrace_2(ray, hitInfo))
    {
        
        if(hitInfo.z>0 && hitInfo.z < t_max){
            
            return 0.0; //occluded
        }
        
    }
    return 1.0; //direct
}

bool Box::IntersectRay(const Ray &ray, float t_max) const
{
    Ray r = ray;
    float tmin  =  -t_max;
    float tmax = t_max;
    //if ray is inside box - return true
    if (IsInside(r.p)) return true;
    //get pairs of planes - x , y, z
    // 0:(x_min,y_min,z_min), 1:(x_max,y_min,z_min)
    // 2:(x_min,y_max,z_min), 3:(x_max,y_max,z_min)
    // 4:(x_min,y_min,z_max), 5:(x_max,y_min,z_max)
    // 6:(x_min,y_max,z_max), 7:(x_max,y_max,z_max)
    float xl = Corner(0).x;
    float xh = Corner(3).x;
    float yl = Corner(0).y;
    float yh = Corner(2).y;
    float zl = Corner(0).z;
    float zh = Corner(5).z;
    
    //Check intersection for X planes
    if(r.p.x != 0.0){
        float tx1 = (xl - ray.p.x)/ ray.dir.x;
        float tx2 = (xh - ray.p.x)/ ray.dir.x;
        
        tmin = max(tmin, min(tx1,tx2));
        tmax = min(tmax, max(tx1, tx2));
    }
    //Check intersection for Y planes
    if(r.p.y != 0.0){
        float tx1 = (yl - ray.p.y)/ ray.dir.y;
        float tx2 = (yh - ray.p.y)/ ray.dir.y;
        
        tmin = max(tmin, min(tx1,tx2));
        tmax = min(tmax, max(tx1, tx2));
    }
    
    //Chec intersection for Z planes
    if(r.p.z != 0.0){
        float tx1 = (xl - ray.p.z)/ ray.dir.z;
        float tx2 = (xh - ray.p.z)/ ray.dir.z;
        
        tmin = max(tmin, min(tx1,tx2));
        tmax = min(tmax, max(tx1, tx2));
    }
    
    return tmax>=tmin;
    
}

Color MtlBlinn::Shade(const Ray &ray, const HitInfo &hInfo, const LightList &lights, int bounceCount) const{
    float bias = BIAS_SHADING;
    Color shade;
    Color rShade = Color(0,0,0);
    Color tShade = Color(0,0,0);
    const Material *mat;
    mat = hInfo.node->GetMaterial();
    const MtlBlinn* mb =static_cast<const MtlBlinn*>(mat);
//    cout<<"HInfo front: "<<hInfo.front<<endl;
    /* local copy */
    Point3 P;
    P.Set(hInfo.p.x,hInfo.p.y,hInfo.p.z);
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
//            temp.front = false;
            n2 = ior;
//            cout<<"outside "<<endl;
        }
        else if(!hInfo.front){ /* Transmission from the inside */
//            temp.front = true;
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


bool RayTrace_2(const Ray &ray, HitInfo &hitInfo)
{
//    cout<<"Ray trace "<<endl;
    return TraceNode(ray, hitInfo, rootNode);
  
}

bool TraceNode(const Ray &r, HitInfo &hInfo,const Node &node)
{
//    cout<<"Tracenode"<<endl;
    Ray ray = r;
    ray = node.ToNodeCoords(r);
//    cout<<"...To node coords"<<endl;
    const Object* obj = node.GetObject();
    bool objHitTest = false;
    bool childHit = false;
    
    if(obj)
    {
//        if(!hInfo.front){ hitSide = HIT_FRONT_AND_BACK; } /* Back face hit for refraction */
        
        if(obj->IntersectRay(ray, hInfo)){
            
            objHitTest=true;
            hInfo.node = &node;
        }
    }
    if (node.GetNumChild() > 0)
    {

        for (int i = 0; i < node.GetNumChild(); i++)
        {
            
//            cout<<"Child "<<i<<endl;
            Ray r = ray;
            const Node &childNode = *node.GetChild(i);
            if(TraceNode(r, hInfo, childNode)){
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

int main(int argc, char* argv[])
{
    const char* filename = "/Users/varunk/Desktop/RayTracerProj1/RayTracePrj5/RayTracePrj5/scene.xml";
    LoadScene(filename);
    
    glutInit(&argc,argv);
    ShowViewport();
   tp.destroy_threadpool();
	return 0;
}

