#include "tinyxml/tinyxml.h"
#include <GLUT/GLUT.h>
#include "viewport.cpp"
#include "scene.h"
#include "xmlload.cpp"
#include <iostream>
#include <math.h>
#include <cmath>
#include "cyColor.h"
#include "cyIrradianceMap.h"
#include <pthread.h>
#include <unistd.h>
#include <vector>
#include <functional>
#include "threadpool.h"
#include "texture.cpp"
#include "PhotonMap.h"
using namespace std;

#define _USE_MATH_DEFINES
#define NUM_THREADS  2

#define BIAS_SHADOW 1.e-3f
#define BIAS_SHADING 1.e-3f
#define BIAS_GLOBAL_AMB 0.01
#define H_BASE_1 2
#define H_BASE_2 3
#define SHADE_THRESHOLD 1e-4f
#define MIN_N_SAMPLES 8
#define MAX_N_SAMPLES 32
#define MIN_SHADOW_SAMPLES 2
#define MAX_SHADOW_SAMPLES 3

#define USE_PHOTONMAP_DIRECT false
#define USE_PHOTONMAP_INDIRECT true
#define MAX_PHOTON_BOUNCE_COUNT 10
#define NUM_EMITTED_PHOTONS 1000000
#define MAX_PHOTONS_STORED 1000000
#define PHOTON_MAX_DIST 5.0f
#define MAX_PHOTONS_CONSIDERED 100
#define PHOTON_BOUNCE_COUNT 5

struct RenderParams{
    Point3 K;
    int pixIndex;
    Point2 pixLocation;
    bool renderComplete = false;
    Cone ray;
    Point2 PixParams;
    vector<Point3> ConfCirclePts;
};

struct ImageParams{
    unsigned int NUM_PIXELS;
    vector<int> PixIndex;
    vector<bool> rendered;
    vector<Point2> PixLocation;
    vector<Point3> K;
    vector<Cone> Ray;
    vector<Point3> SampleLoc;
    vector<Point2> PixParams;
    vector<vector<Point3 >> ConfusionCirclePoints;
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
TexturedColor background;
TexturedColor environment;
TextureList textureList;
cyPoint3f _f;
cyPoint3f _s;
cyPoint3f _u;
BalancedPhotonMap *_bmap;
int _numberOfPhotonsEmitted = 0;
pthread_t threadID[NUM_THREADS];
pthread_mutex_t getPix_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t setPix_mutex = PTHREAD_MUTEX_INITIALIZER;

RenderParams args[NUM_THREADS];
Color24 white(cyColor({255,255,255}));
Color24 black(cyColor({0,0,0}));
PhotonMap *map;



RenderParams giveMeAPixelToRender();
bool TraceNode(const Cone &r,HitInfo &hitInfo,const Node &curNode);
bool RayTrace(HitInfo &hitInfo, Node* node,Ray ray,int PixIndex);
bool RayTrace_2(const Cone &ray, HitInfo &hitInfo);
void PopulateImageParams();
void doRender(void* arg);

void PopulateImageParams()
{
    _f = camera.dir;
    _f.Normalize();
    _s = _f.Cross(camera.up);
    _s.Normalize();
    _u = _s.Cross(_f);
    
    cout<<"Populating ImageParams..."<<endl;
    float alpha = camera.fov;
    float l = camera.focaldist;
    float h = l * tan(alpha/2.0 *(M_PI/180.0));
    
    float aspectRatio = (float)camera.imgWidth/camera.imgHeight;
    float s = aspectRatio * abs(h);
    float dx = (2 * abs(s))/camera.imgWidth;
    float dy = -(2 * abs(h))/camera.imgHeight;
    float dxx = dx/2.0 , dyy=dy/2.0;
    Point3 K(-s,h,-l);
    K.x += (dxx );
    K.y += (dyy );
    
    for(int i = 0; i< camera.imgHeight ; i++){
        for(int j = 0; j< camera.imgWidth; j++){
            K.x += dx;
            Matrix3 RotMat;
            const float pts[9]={_s.x,_u.x,-_f.x,_s.y,_u.y,-_f.y,_s.z,_u.z,-_f.z};
            RotMat.Set(pts);
//            K = RotMat*K;
            Cone r = Cone(camera.pos, K);
            r.dir = r.dir * RotMat;
            r.dir.Normalize();
            r.radius = 0.0;
            r.tanAngle = tan(abs(dyy)/(float)l);
            /* Populating the Struct */
            Point2 pixLoc = Point2(j,i);
            
            imageParams.K.push_back(K);
            imageParams.rendered.push_back(false);
            imageParams.PixLocation.push_back(pixLoc);
            imageParams.PixIndex.push_back( i * camera.imgWidth + j);
            imageParams.Ray.push_back(r);
            Point2 pixDimensions = Point2(dx,dy);
            imageParams.PixParams.push_back(pixDimensions);
            vector<Point3> ConfCirclePts;
           
            float randAng = rand()/ (float) RAND_MAX;
            randAng *= M_PI * 2.0;
            for(int i = 1; i<=MAX_N_SAMPLES; i++){
                float hx = camera.dof * Halton(i, H_BASE_1);
                float hy = Halton(i, H_BASE_2);
                
                float r =  sqrtf(hx);
                float theta = hy * M_PI * 2.0 + randAng;
                
                float x =  r * cosf(theta);
                float y =  r * sinf(theta);
                Point3 newCamPos(x, y, 0);
                newCamPos = newCamPos * RotMat;
                newCamPos += camera.pos;
                ConfCirclePts.push_back(newCamPos);
            }
            
            imageParams.ConfusionCirclePoints.push_back(ConfCirclePts);
        }
        K.x = -s;
        K.x += dxx;
        K.y += dy;
    }
    
   

}
void ComputeIrradianceMap()
{
    // irradMap.initialize();
    
}

Ray PointLight::RandomPhoton() const
{
    Ray ray;
    float r = 1.0f;
    float theta = 0.0f; float phi = 0.0f;
    float f_U = (float)rand()/(float)RAND_MAX;
    float f_V = (float)rand()/(float)RAND_MAX;
    theta	=	2.0 * M_PI * f_U;
    phi	=	acos(2 * f_V -1);
//    theta = ((float)rand()/(float)RAND_MAX) * 2.0 * M_PI;
//    phi   = ((float)rand()/(float)RAND_MAX) * 2.0 * M_PI;
    
    float x = r * sin(phi) * cos(theta);
    float y = r * sin(phi) * sin(theta);
    float z = r * cos(phi);
    
    Point3 randomPointOnSphere = Point3(x,y,z) + position; //translate
    ray.p = position;
    ray.dir = randomPointOnSphere - position; //target - souce
    ray.dir.Normalize();
    
    
    return ray;
}
bool MtlBlinn::RandomPhotonBounce(Ray &r, Color &c, const HitInfo &hitInfo) const
{
    float f_absorptionProb, f_diffuseRefProb, f_specularRefProb, f_transProb, f_randNum, f_totalProb;
    bool newRaySet = false;
    int outCome = -1; //0:abs, 1: diffuse reflection, 2: specular, 3: trans
    Color chosenColor(0, 0, 0);
    f_diffuseRefProb = diffuse.GetColor().Grey();
    f_specularRefProb = specular.GetColor().Grey();
    f_transProb = refraction.GetColor().Grey();
    f_absorptionProb = absorption.Grey();

   // f_absorptionProb = 1 - (f_diffuseRefProb + f_specularRefProb + f_transProb);
//    f_absorptionProb  = f_absorptionProb < 0? 0.1:f_absorptionProb;
//    f_absorptionProb  = f_absorptionProb > 1? 0.9:f_absorptionProb;
    f_totalProb = f_absorptionProb + f_diffuseRefProb + f_specularRefProb + f_transProb;
    
    f_diffuseRefProb /= f_totalProb;
    f_specularRefProb /= f_totalProb;
    f_transProb /= f_totalProb;
    f_absorptionProb /= f_totalProb;
//        cout<<"diff: "<<f_diffuseRefProb<<"spec: "<<f_specularRefProb<<"trans: "<<f_transProb<<"ABS: "<<f_absorptionProb<<endl;
    f_randNum = (float)rand()/(float)RAND_MAX;
    
    if(f_randNum < f_absorptionProb){
        outCome = 0;
    }
    else if(f_randNum <= (f_absorptionProb + f_transProb) && f_randNum > f_absorptionProb){
        outCome = 3;

    }
    else if(f_randNum <= (f_absorptionProb + f_transProb + f_specularRefProb) && f_randNum > (f_absorptionProb + f_transProb)){
        outCome = 2;
    }
    else if(f_randNum <= (f_absorptionProb + f_transProb + f_specularRefProb + f_diffuseRefProb ) && f_randNum > (f_absorptionProb + f_transProb + f_specularRefProb)){
        outCome = 1;
    }
   // if(f_transProb == 1) outCome = 3;

//         cout<<outCome<<endl;
    switch (outCome) {
        case 0:
            //kill ray
            newRaySet = false;
            break;
        case 1:
        {
            //diffuse reflection
            //create random point over hemisphere
            Point3 zaxis = Point3(hitInfo.N.x, hitInfo.N.y, hitInfo.N.z), axis, xaxis, yaxis, dir;
            
            //choose an axis with smallest component in normal
            if (zaxis.x < zaxis.y  && zaxis.x < zaxis.z)
            {
                axis = Point3(1.0, 0.0, 0.0);
            }
            else if (zaxis.y < zaxis.z)
            {
                axis = Point3(0.0, 1.0, 0.0);
            }
            else
                axis = Point3(0.0, 0.0, 1.0);
            
            xaxis = zaxis.Cross(axis).GetNormalized();
            yaxis = zaxis.Cross(xaxis).GetNormalized();

            float theta = ((float)rand()/(float)RAND_MAX) * (M_PI_2);
            float phi = (float)rand()/((float)RAND_MAX/ (2.0 * M_PI));
            
            dir = (xaxis * (cos(phi) * sin(theta))) + ( yaxis * (sin(phi) * sin(theta))) + (zaxis * (cos(theta)));
            dir.Normalize();
            
            //create ray
            Ray ray;
            ray.p = hitInfo.p + BIAS_GLOBAL_AMB * dir;
            ray.dir = dir;
            r = ray;

            newRaySet = true;
        }
            break;
        
        case 2: {//specular reflection
            Point3 N = hitInfo.N;
            Point3 V = Point3(r.p.x -  hitInfo.p.x, r.p.y - hitInfo.p.y, r.p.z - hitInfo.p.z);
            // V.Normalize();
            Point3 VR = 2 * V.Dot(N) * N - V;
            Ray ray = Ray(hitInfo.p + BIAS_SHADING * VR, VR);
            ray.dir.Normalize();
            
            r = ray;
            newRaySet = true;
        }
            break;
        
        case 3:{
            //transmission
            Point3 N = hitInfo.N;
            Point3 V = Point3(hitInfo.p.x - r.p.x, hitInfo.p.y - r.p.y, hitInfo.p.z - r.p.z);
            V.Normalize();
            float n1 = 1, n2 = 1;
            if(hitInfo.front){ /* Hitting from outside */
                //            temp.front = false;
                n2 = ior;
                //            cout<<"outside "<<endl;
            }
            else if(!hitInfo.front){ /* Transmission from the inside */
                //            temp.front = true;
                n1 = ior;
                //            cout<<"intside... "<<endl;
                //            N = -hInfo.N;
                N *= -1;
            }
            float ratio_n = n1 / n2;
            
            float costheta_v = -V.Dot(N);        /* refer: http://graphics.stanford.edu/courses/cs148-10-summer/docs/2006--degreve--reflection_refraction.pdf */
            
            float sin2theta_t = ratio_n * ratio_n * (1 - costheta_v * costheta_v);
            Point3 T =   ratio_n * V + (ratio_n * costheta_v - sqrtf(1 - sin2theta_t)) * N ;
            //        cout<<ratio_n<<" "<<"cos_v "<<costheta_v<<" sin2theta_t "<<sin2theta_t<<endl;
            Ray tRay = Ray(hitInfo.p,T);
            
            tRay.dir.Normalize();
            tRay.p = tRay.p + BIAS_SHADING * tRay.dir;
            tRay.dir.Normalize();
            r = tRay;
            newRaySet = true;
        }
            break;
        
        default:
            break;
    }
    //generate random ray
    Color newColor = c * diffuse.GetColor(); //modulate with color - changing the photon
    c = newColor;
    return newRaySet;
}
void DoRecursivePhotonBounce(HitInfo hInfo, Ray photonRay, Color color, int &numStoredPhotons, int bouncesAvailable)
{
    
    //spawn new ray based on prob of abs, refl_diffuse, refl_spec, trans.
    int bounceCount = bouncesAvailable;
    const MtlBlinn *material = static_cast<const MtlBlinn *>(hInfo.node->GetMaterial());
    if(bounceCount>0){
        if(material->RandomPhotonBounce(photonRay, color, hInfo)){ //material has bounced the ray: ray, color have been changed

           HitInfo hitInfo;
            hitInfo.Init();
            if(RayTrace_2(photonRay, hitInfo)){ //ray has hit something after second bounce

                float pos[3], dir[3], power[3];
                hitInfo.p.GetValue(pos); //hit point of current hit
                photonRay.dir.GetValue(dir); //direction of photon pre-bounce
                Color intensity = color; //color of the photon pre-bounce
                intensity.GetValue(power);
                const MtlBlinn *hitMaterial = static_cast<const MtlBlinn *>(hitInfo.node->GetMaterial());
                if(hitMaterial->IsPhotonSurface()){
                    storePhoton(map, power, pos, dir);
                    numStoredPhotons++;
                }
                bounceCount--;
                DoRecursivePhotonBounce(hitInfo, photonRay, color, numStoredPhotons, bounceCount);
            }
        }
    }
    //store subsequenct hits - decrement bounce count
}

void CreatePhotonMap()
{
    cout<<"Creating Photon map..."<<endl;
    map = createPhotonMap(NUM_EMITTED_PHOTONS);
    Color clr_totalIntensity(0, 0, 0);
    LightList list_PointLights{};
    list_PointLights.DeleteAll();
    for(int i=0; i < lights.size(); i++){
        if(lights.at(i)->IsAmbient()) continue;
        if(PointLight * curLight = static_cast<PointLight *>( lights.at(i)))
        {

            clr_totalIntensity += lights.at(i)->GetPhotonIntensity();
            list_PointLights.push_back(lights.at(i));
        }
    }
   // std::sort(list_PointLights.begin(),list_PointLights.end());
    /* Right now implemented for single pointlight source */
    int numStoredPhotons = 0, numEmittedPhotons = 0;
    int bounceCount = PHOTON_BOUNCE_COUNT;
    HitInfo hInfo;
    while(numStoredPhotons <= MAX_PHOTONS_STORED){
        PointLight *_curLight = (PointLight *)(list_PointLights.at(0));
        Ray photonRay = _curLight->RandomPhoton();
        numEmittedPhotons++;
        hInfo.Init();
        
        //ray trace
        
        //don't store first hit

        if(RayTrace_2((Cone)photonRay, hInfo)){
            float pos[3], dir[3], power[3];
            hInfo.p.GetValue(pos);
            photonRay.dir.GetValue(dir);
            Color intensity = _curLight->GetPhotonIntensity();
            intensity.GetValue(power);
            if(USE_PHOTONMAP_DIRECT){
                const MtlBlinn *mtl = static_cast<const MtlBlinn *>(hInfo.node->GetMaterial());
                if (mtl->IsPhotonSurface()) {
                    storePhoton(map, power, pos, dir);
                    numStoredPhotons++;
                }
            }
            if(USE_PHOTONMAP_INDIRECT){
                DoRecursivePhotonBounce(hInfo,photonRay,intensity,numStoredPhotons,bounceCount);
            }
        }
        _curLight = NULL;
    }
    
    
    list_PointLights.clear();
    BalancedPhotonMap *bmap = balancePhotonMap(map);
    char *filename = "/Users/varunk/Desktop/RayTracerProj1/RayTracePrj11/RayTracePrj11/photonmap.dat";
    savePhotonMap(bmap, filename);
    _bmap = loadPhotonMap(filename);
    _numberOfPhotonsEmitted = numEmittedPhotons;
}

void BeginRender()
{
    //Load ImageParams
    imageParams.NUM_PIXELS = camera.imgHeight * camera.imgWidth;
    PopulateImageParams();
    
    //Compute Irradiance Map
//    ComputeIrradianceMap();
    if(USE_PHOTONMAP_DIRECT || USE_PHOTONMAP_INDIRECT) {
        CreatePhotonMap();
    }
    cout<<"Beginning to render..."<<endl;
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    srand(time(NULL));
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
        x->K = imageParams.K.at(i);
        x->ConfCirclePts = imageParams.ConfusionCirclePoints.at(i);
        x->PixParams = imageParams.PixParams.at(i);
        Task* t = new Task(&doRender, (void*) x);
        //    cout << "Adding to pool, task " << i+1 << endl;
        tp.add_task(t);
//            cout << "Added to pool, task " << i+1 << endl;
    }
    
//    usleep(20);
    
    
    }

    
    
}
Color AverageShades(vector<Color> shades, int n)
{
    Color shade(0,0,0);
    for(int i=0; i<shades.size(); i++){
        shade += shades.at(i);
    }
    shade /= n;
    return shade;
}

bool VarianceOverThreshold(vector<Color> shades)
{
    float sumR=0.0, sumB=0.0, sumG=0.0, Ninv=0,avgR=0,avgG=0,avgB=0;
    Ninv = 1.0 / MIN_N_SAMPLES;
    for(int i=0; i<shades.size(); i++){
        
        sumR += shades.at(i).r;
        sumG += shades.at(i).g;
        sumB += shades.at(i).b;
        
        
    }
    
    avgR = Ninv * sumR;
    avgG = Ninv * sumG;
    avgB = Ninv * sumB;
    sumR=sumG=sumB=0.0f;
    for(int i = 0; i< shades.size(); i++){
        sumR += pow(shades.at(i).r - avgR,2.0);
        sumG += pow(shades.at(i).g - avgG,2.0);
        sumB += pow(shades.at(i).b - avgB,2.0);
    }
    
    
    float varR = sumR /(float)(MIN_N_SAMPLES-1);
    float varG = sumG /(float)(MIN_N_SAMPLES-1);
    float varB = sumB /(float)(MIN_N_SAMPLES-1);
    
    if(varR > SHADE_THRESHOLD || varG > SHADE_THRESHOLD || varB > SHADE_THRESHOLD){
        return true;
    }
    return false;
    
}
void doRender(void* arg){
    
    RenderParams rarg = *((RenderParams *)arg);
    //cout<<"Do render...."<<endl;
            bool pixelHit=false;

    
            HitInfo hitInfo;
            hitInfo.Init();
            
            Point2 pixLoc = rarg.pixLocation;
            Cone r = rarg.ray;
            int PixIndex = rarg.pixIndex;
            Color shade(0,0,0);
    
    vector<Point2> haltonXY;
    float dx = rarg.PixParams.x;
    float dy = rarg.PixParams.y;
    float x=0;
    float y=0;
    _f.Normalize();
    _s.Normalize();
    const float pts[9]={_s.x,_u.x,-_f.x,_s.y,_u.y,-_f.y,_s.z,_u.z,-_f.z};

    Matrix3 RotMat;
    RotMat.Set(pts);
    for(int i=0; i < MIN_N_SAMPLES; i++){
        
        x = dx * Halton(i+1, H_BASE_1);
        y = dy * Halton(i+1, H_BASE_2);
        
        if(x > dx * 0.5) { x -= dx; }
        if(y < dy * 0.5) { y -= dy; }
       
        x += rarg.K.x;
        y += rarg.K.y;
        
        Point2 sampleLoc = Point2(x,y);
        haltonXY.push_back(sampleLoc);
    }
    
    vector<Color> shades;
            if(rootNode.GetNumChild()>0){
                for(int i=0; i< MIN_N_SAMPLES;i++){
                    
                    Point3 sampleDir = Point3(haltonXY.at(i).x, haltonXY.at(i).y, rarg.K.z);
                    
                    int rindex = rand() % MAX_N_SAMPLES; //rindex = i;
                    Point3 randPos = rarg.ConfCirclePts.at(rindex);
                    Cone sampleRay = Cone(randPos, sampleDir);
                    sampleRay.dir = sampleRay.dir * RotMat;
                    sampleRay.dir -= randPos - camera.pos;
                    sampleRay.dir.Normalize();
                    sampleRay.radius = r.radius;
                    sampleRay.tanAngle = r.tanAngle;
                    r = sampleRay;
                    
                    if(RayTrace_2(r, hitInfo)) {
                        pixelHit=true;
                        shade = hitInfo.node->GetMaterial()->Shade(r, hitInfo, lights, 8);
                        shades.push_back(shade);
                    }
                    hitInfo.Init();
                }

                if(VarianceOverThreshold(shades)){
                    renderImage.SetSampleCountPixel(PixIndex, 255);
                    
                    hitInfo.Init();
                    
                    for(int i=MIN_N_SAMPLES; i < MAX_N_SAMPLES; i++){
                        
                        x = dx * Halton(i+1, H_BASE_1);
                        y = dy * Halton(i+1, H_BASE_2);
    
                        if(x > dx * 0.5) { x -= dx;}
                        if(y < dy * 0.5) { y -= dy;}
                        
                        x += rarg.K.x;
                        y += rarg.K.y;
                        
                        Point2 sampleLoc = Point2(x,y);
                        haltonXY.push_back(sampleLoc);
                        
                        Point3 sampleDir = Point3(haltonXY.at(i).x, haltonXY.at(i).y, rarg.K.z);
                        
                        int rindex = rand() %  MAX_N_SAMPLES; //rindex = i;
                        Point3 randPos = rarg.ConfCirclePts.at(rindex);
                        Cone sampleRay = Cone(randPos, sampleDir);
                        sampleRay.dir = sampleRay.dir * RotMat;
                        sampleRay.dir -= randPos - camera.pos;
                        
                        sampleRay.dir.Normalize();
                        sampleRay.radius = r.radius;
                        sampleRay.tanAngle = r.tanAngle;
                        r = sampleRay;
                        
                        if(RayTrace_2(r, hitInfo)) {
                            pixelHit=true;
                            shade = hitInfo.node->GetMaterial()->Shade(r, hitInfo, lights, 5);
                            shades.push_back(shade);
                        }
                        hitInfo.Init();
                    }
                    shade = AverageShades(shades, (int)shades.size());
                }
                else{
                    shade = AverageShades(shades, (int)shades.size());
                    renderImage.SetSampleCountPixel(PixIndex, 0);
                }
                shade.r = pow(shade.r, 1.0/2.2);
                shade.g = pow(shade.g, 1.0/2.2);
                shade.b = pow(shade.b, 1.0/2.2);
                renderImage.PutPixel(PixIndex,shade, hitInfo.z);
            }
            
            if(!pixelHit){
                Point3 uvw(pixLoc.x/ camera.imgWidth,pixLoc.y/camera.imgHeight,0);
                shade = background.Sample(uvw);
                renderImage.PutPixel(PixIndex, shade, BIGFLOAT);
            
            }
    pthread_mutex_lock(&setPix_mutex);
    renderImage.IncrementNumRenderPixel(1);
    pthread_mutex_unlock(&setPix_mutex);

}
Color PointLight::Illuminate(const Point3 &hitPt ,const Point3 &normal) const
{
    //create a ray towards light
    Point3 dir = position - hitPt;
    //dir.Normalize();
    
    Ray sRay = Ray(hitPt, dir);
//    return intensity * Shadow(sRay, 1);
    
    Point3 xAxis(1,0,0), yAxis(0,1,0), v1;
    
    if(dir.Dot(xAxis) > 0.8) {
        v1 = yAxis.Cross(dir);
    }
    else{
        v1 = xAxis.Cross(dir);
    }
    
    float shadow = 0.0;
    Point3 v2 = v1.Cross(dir);
    v2.Normalize();
    v1.Normalize();
    Point3 xv1 = v1;
    Point3 yv2 = v2;
//    srand(time(NULL));
    float random = 0.0;

    for(int i =0; i< MIN_SHADOW_SAMPLES; i++){
        random = rand() / (float) RAND_MAX;
        float rRadius = sqrtf(random) * size;
        random = rand() / (float) RAND_MAX;
        float rAngle = random * (2.0 * M_PI);
        float xv = rRadius * cos(rAngle);
        float yv = rRadius * sin(rAngle);
        xv1 *= xv;
        yv2 *= yv;
        
    
    sRay.dir = (position + xv1.Length()+ yv2.Length() ) - hitPt;
//    sRay.dir.Normalize();
    shadow += Shadow(sRay, 1);
        xv1 = v1; yv2 = v2;
    }
    shadow /= (float)MIN_SHADOW_SAMPLES;
    if(shadow != 0.0 && shadow != 1.0)
    {
//        cout<<"SHADOW: "<<shadow<<endl;
        shadow = 0.0;
        for(int i = 0; i< MAX_SHADOW_SAMPLES; i++){
            random = rand() / (float) RAND_MAX;
            float rRadius = sqrtf(random) * size;
            random = rand() / (float) RAND_MAX;
            float rAngle = random * (2.0 * M_PI);
            float xv = rRadius * cos(rAngle);
            float yv = rRadius * sin(rAngle);
            xv1 *= -xv;
            yv2 *= -yv;
            
            sRay.dir = (position +xv1.Length() + yv2.Length() ) - hitPt;
//            sRay.dir.Normalize();
            shadow += Shadow(sRay, 1);
            xv1 = v1; yv2 = v2;
        }
            shadow /= (float)(MAX_SHADOW_SAMPLES);
    }

    /*inverse square fall off*/
    Point3 d = position - hitPt;
    float distance = d.LengthSquared();
    
    return (1.0 / distance) * intensity * shadow;
    
}

float GenLight::Shadow(Ray ray, float t_max)
{
   // cout<<"Calculating shadow"<<endl;
    //add bias
//    float eps = BIAS_SHADOW;
//    ray.p = Point3(ray.p.x+eps*ray.dir.x, ray.p.y+eps*ray.dir.y,ray.p.z+eps*ray.dir.z);
    ray.p = ray.p + BIAS_SHADOW * ray.dir;
    HitInfo hitInfo;
    
    if(RayTrace_2(ray, hitInfo))
    {
        
        if(hitInfo.z>0 && hitInfo.z < t_max){
            
            return 0.0; //occluded
        }
        
    }
    return 1.0; //direct
}
Color AmbientLight::Illuminate(const Point3 &P, const Point3 &N)const
{
        Color color = Color(0.0, 0.0, 0.0);
    
    
        Point3 zaxis = Point3(N.x, N.y, N.z), axis, xaxis, yaxis, dir;

        //choose an axis with smallest component in normal
        if (zaxis.x < zaxis.y  && zaxis.x < zaxis.z)
        {
            axis = Point3(1.0, 0.0, 0.0);
        }
        else if (zaxis.y < zaxis.z)
        {
            axis = Point3(0.0, 1.0, 0.0);
        }
        else
            axis = Point3(0.0, 0.0, 1.0);
        
        xaxis = zaxis.Cross(axis).GetNormalized();
        yaxis = zaxis.Cross(xaxis).GetNormalized();
//        for (int i=0; i < 1 ; i++)
        {
            float X = (float)rand()/(float)RAND_MAX;
            float theta = 0.5 * acos(1.0 - (2.0 * X));
            float phi = (float)rand()/((float)RAND_MAX/ (2.0 * M_PI));
            
            dir = (xaxis * (cos(phi) * sin(theta))) + ( yaxis * (sin(phi) * sin(theta))) + (zaxis * (cos(theta)));
            dir.Normalize();
            //Trace the ray
            HitInfo hInfo;
            hInfo.Init();

            Cone ray;
            ray.p = P + BIAS_GLOBAL_AMB * dir;
            ray.dir = dir;
            //ray.dir.Normalize();
            if(RayTrace_2(ray, hInfo))
               {
                   if(USE_PHOTONMAP_INDIRECT){
                       const Material *material = hInfo.node->GetMaterial();
                       
                       float f_irrad[3], f_pos[3],f_normal[3], f_maxDist;
                       int i_nPhotons;
                       hInfo.p.GetValue(f_pos);
                       hInfo.N.GetValue(f_normal);
                       f_maxDist = PHOTON_MAX_DIST;
                       i_nPhotons = MAX_PHOTONS_CONSIDERED;
                       irradianceEstimate(_bmap, f_irrad, f_pos, f_normal, f_maxDist, i_nPhotons);
                       f_irrad[0] *= (4.0 * M_PI)/(_numberOfPhotonsEmitted);
                       f_irrad[1] *= (4.0 * M_PI)/(_numberOfPhotonsEmitted);
                       f_irrad[2] *= (4.0 * M_PI)/(_numberOfPhotonsEmitted);
                       color += Color(f_irrad[0],f_irrad[1],f_irrad[2]);

                       color *= static_cast<const MtlBlinn*>(material)->Shade (ray, hInfo, lights, -10, true);
                       
                   }
                   else{
                       const Material *material = hInfo.node->GetMaterial();
                       color += static_cast<const MtlBlinn*>(material)->Shade (ray, hInfo, lights, 1, true);
                   }
               }
               else
               {
                   color += environment.SampleEnvironment(dir);
               }
         }
    

               return color;
}





Color MtlBlinn::Shade(const Cone &ray, const HitInfo &hInfo, const LightList &lights, int bounceCount, bool globalAmbient) const
{
    float bias = BIAS_SHADING;
    Color shade;
    Color rShade = Color(0,0,0);
    Color tShade = Color(0,0,0);
    const Material *mat;
    mat = hInfo.node->GetMaterial();
    const MtlBlinn* mb =static_cast<const MtlBlinn*>(mat);
    
    /* local copy */
    Point3 P;
    P.Set(hInfo.p.x,hInfo.p.y,hInfo.p.z);
    Cone iRay = ray;
    
    Color ambInt = mb->diffuse.Sample(hInfo.uvw, hInfo.duvw);
    Color allOther = Color(0,0,0);
    Color diffuse = mb->diffuse.Sample(hInfo.uvw, hInfo.duvw);
    Color ambComponent = Color(0,0,0);
    Point3 newN = hInfo.N;
    if(bounceCount == -10){
//        cout<<bounceCount<<endl;
        return diffuse;
    }
    
    for ( unsigned int i=0; i<lights.size(); i++ ) {
        Light *curLight = lights.at(i);
        if(curLight->IsAmbient()){
            Color intensity;
            if(globalAmbient == false){
                intensity = lights.at(i)->Illuminate(hInfo.p, hInfo.N);
                ambComponent += (ambInt * intensity);
            }
            
//          continue;
        }
        else{
            if(USE_PHOTONMAP_DIRECT){
                float f_irrad[3], f_pos[3],f_normal[3], f_maxDist;
                int i_nPhotons;
                hInfo.p.GetValue(f_pos);
                hInfo.N.GetValue(f_normal);
                f_maxDist = PHOTON_MAX_DIST;
                i_nPhotons = MAX_PHOTONS_CONSIDERED;
                irradianceEstimate(_bmap, f_irrad, f_pos, f_normal, f_maxDist, i_nPhotons);
                f_irrad[0] *= (4.0 * M_PI)/(_numberOfPhotonsEmitted);
                f_irrad[1] *= (4.0 * M_PI)/(_numberOfPhotonsEmitted);
                f_irrad[2] *= (4.0 * M_PI)/(_numberOfPhotonsEmitted);
                allOther += Color(f_irrad[0],f_irrad[1],f_irrad[2]);
                allOther *= diffuse; //modulate the color with the irradiance
                
            }
            else{
                
            Point3 L = -lights[i]->Direction(P);
            L.Normalize();
            Point3 V = ray.p - P;
            V.Normalize();
            Point3 LplusV = L + V;
            Point3 H = (L+V)/LplusV.Length();
            H.Normalize();
            float alpha = mb->glossiness;
            Point3 N = newN;
            float S = H.Dot(N);
            S = pow(S,alpha);
            float costheta = L.Dot(N)/(L.Length() * N.Length());
            Color intensity = lights[i]->Illuminate(P, hInfo.N);
//            cout<<"costheta "<<endl;
            allOther += intensity * (costheta>0?costheta:0) * (diffuse + S * (mb->specular.Sample(hInfo.uvw, hInfo.duvw))) ;
            }
        }
        /* finally add inta*cola + intall*costheta*(cold + s* colS)*/
        shade = ambComponent  + allOther;
    }
    
    /* Calculate refraction */
    if(refraction.GetColor().r>0 && bounceCount>0){
        //compute new jittered normal
        float gloss = refractionGlossiness;
        if(gloss){
            float random = rand()/(float)RAND_MAX;
            float rRadius = sqrtf(random) * gloss;
            random = rand()/(float)RAND_MAX;
            float rAngle =  random * 2.0 * M_PI;
            float x = rRadius * cos(rAngle);
            float y = rRadius * sin(rAngle);
            Point3 xAxis(1,0,0), yAxis(0,1,0), v1, v2, normalDir;
            normalDir = hInfo.N;
            //    normalDir.Normalize();
            if(normalDir.Dot(xAxis) > 0.7)  v1 = normalDir.Cross(yAxis);
            else v1 = normalDir.Cross(xAxis);
            v2 = v1.Cross(normalDir);
            v1.Normalize(); v2.Normalize();
            v1 *= x;
            v2 *= y;
            
            newN = hInfo.N + v1.Length() + v2.Length();
            newN.Normalize();
        }
        else{
            newN = hInfo.N;
        }
        //-------------------------------------
        
        Color reflShade = Color(0,0,0);
        float R0, Refl = 0.0f, Trans = 0.0f;
        HitInfo temp;
        temp.Init();
        
//        Point3 N = hInfo.N;
        Point3 N = newN;
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
//            N = -hInfo.N;
            N *= -1;
        }
        float ratio_n = n1 / n2;
        
        float costheta_v = -V.Dot(N);        /* refer: http://graphics.stanford.edu/courses/cs148-10-summer/docs/2006--degreve--reflection_refraction.pdf */

        float sin2theta_t = ratio_n * ratio_n * (1 - costheta_v * costheta_v);
        Point3 T =   ratio_n * V + (ratio_n * costheta_v - sqrtf(1 - sin2theta_t)) * N ;
//        cout<<ratio_n<<" "<<"cos_v "<<costheta_v<<" sin2theta_t "<<sin2theta_t<<endl;
        Cone tRay = Cone(hInfo.p,T);
        
        //tRay.dir.Normalize();
        tRay.p.x = tRay.p.x + bias *tRay.dir.x; /* add bias */
        tRay.p.y = tRay.p.y + bias *tRay.dir.y;
        tRay.p.z = tRay.p.z + bias *tRay.dir.z;
//        cout<<"B temp front: "<< temp.front<<endl;
        if(sin2theta_t <= 1){
            if(RayTrace_2(tRay, temp)){ /* ray tracing after refraction */
//                bounceCount--;
                tShade  =  temp.node->GetMaterial()->Shade(tRay,temp,lights,bounceCount);
            }
            else{ /* no hit after refraction */
                tShade = environment.SampleEnvironment(tRay.dir);
                
            }
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
            
            tShade.r *= exp(-absorption.r * temp.z);
            tShade.g *= exp(-absorption.g * temp.z);
            tShade.b *= exp(-absorption.b * temp.z);
        }
        else {/* Total internal reflection */
            Refl = 1.0f;
        }
        
            
        
        
        /* Calculate reflection due to reflectance */
        if(bounceCount >0){
//            N = hInfo.N;
            N = newN;
            Point3 V = Point3(iRay.p.x -  P.x, iRay.p.y - P.y, iRay.p.z - P.z);
            //V.Normalize();
            Point3 VR = 2 * V.Dot(N) * N - V;
            //VR.Normalize();
            Cone rRay = Cone(P + BIAS_SHADING * VR, VR);
            rRay.dir.Normalize();
            HitInfo temp1;
            temp1.Init();
            if(rootNode.GetNumChild()>0){
                if(RayTrace_2(rRay, temp1)){
                    bounceCount --;
                    reflShade = temp1.node->GetMaterial()->Shade(rRay, temp1, lights, bounceCount);
                }
                else{
                    reflShade = environment.SampleEnvironment(rRay.dir);
//                    reflShade = Color(1,100,1);
                }
            }
        }
        
//        cout<<"Refl: "<<Refl<<"Trans "<<Trans<<endl;
        tShade = refraction.GetColor().r * (Trans * tShade + Refl * reflShade);
        
        
    }

    /* calculate reflection*/
    if(reflection.GetColor().r>0 && bounceCount > 0){
        //compute new jittered normal
        float gloss = reflectionGlossiness;
        if(gloss){
            float random = rand()/(float)RAND_MAX;
            float rRadius = sqrtf(random) * gloss;
            random = rand()/(float)RAND_MAX;
            float rAngle =  random * 2.0 * M_PI;
            float x = rRadius * cos(rAngle);
            float y = rRadius * sin(rAngle);
            Point3 xAxis(1,0,0), yAxis(0,1,0), v1, v2, normalDir;
            normalDir = hInfo.N;
            //    normalDir.Normalize();
            if(normalDir.Dot(xAxis) > 0.7)  v1 = normalDir.Cross(yAxis);
            else v1 = normalDir.Cross(xAxis);
            v2 = v1.Cross(normalDir);
            v1.Normalize(); v2.Normalize();
            v1 *= x;
            v2 *= y;
            
            newN = hInfo.N + v1+ v2;
            newN.Normalize();
        }
        else{
            newN = hInfo.N;
        }
        //-------------------------------------
        
        Point3 N = newN;//hInfo.N;
        Point3 V = Point3(iRay.p.x -  P.x, iRay.p.y - P.y, iRay.p.z - P.z);
       // V.Normalize();
        Point3 VR = 2 * V.Dot(N) * N - V;
        Cone rRay = Cone(P + BIAS_SHADING * VR, VR);
        rRay.dir.Normalize();
        HitInfo temp;
        temp.Init();
        if(rootNode.GetNumChild()>0){
            if(RayTrace_2(rRay, temp)){
                bounceCount--;
                rShade = reflection.GetColor().r * temp.node->GetMaterial()->Shade(rRay, temp, lights, bounceCount);
            }
            else{
                rShade =  reflection.GetColor().r *environment.SampleEnvironment(rRay.dir);
//                rShade = Color(1,111,1);
            }
        }
    }
    
  
    
    /* Add shade with reflected and refracted colors */
    shade += (rShade + tShade);
    return shade;
};


bool RayTrace_2(const Cone &ray, HitInfo &hitInfo)
{
//    cout<<"Ray trace "<<endl;
    return TraceNode(ray, hitInfo, rootNode);
  
}

bool TraceNode(const Cone &r, HitInfo &hInfo,const Node &node)
{
//    cout<<"Tracenode"<<endl;
    Cone ray = r;
    if(node.IsInMotion()){
        float timeIns = (rand()/(float)RAND_MAX) * camera.shutter;
        ray = node.ToNodeCoords(r, timeIns);
    }
    else{    ray = node.ToNodeCoords(r); }
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
            Cone r = ray;
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
    cout<<"Stopping Render..."<<endl;
}

int main(int argc, char* argv[])
{
    const char* filename = "/Users/varunk/Desktop/RayTracerProj1/RayTracePrj11/RayTracePrj11/scene.xml";
    LoadScene(filename);
    
    glutInit(&argc,argv);
    ShowViewport();
    
   tp.destroy_threadpool();
	return 0;
}

