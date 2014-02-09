//-------------------------------------------------------------------------------
///
/// \file       texture.cpp 
/// \author     Cem Yuksel (www.cemyuksel.com)
/// \version    1.0
/// \date       October 7, 2013
///
/// \brief Example source for CS 6620 - University of Utah.
///
//-------------------------------------------------------------------------------
 
#include "texture.h"
#include <GLUT/glut.h>
 
//-------------------------------------------------------------------------------
 
int ReadLine( FILE *fp, int size, char *buffer )
{
    int i;
    for ( i=0; i<size; i++ ) {
        buffer[i] = fgetc(fp);
        if ( feof(fp) || buffer[i] == '\n' || buffer[i] == '\r' ) {
            buffer[i] = '\0';
            return i+1;
        }
    }
    return i;
}
 
//-------------------------------------------------------------------------------
 
bool LoadPPM( FILE *fp, int &width, int &height, std::vector<Color24> &data )
{
    const int bufferSize = 1024;
    char buffer[bufferSize];
    ReadLine(fp,bufferSize,buffer);
    if ( buffer[0] != 'P' && buffer[1] != '6' ) return false;
     
    ReadLine(fp,bufferSize,buffer);
    while ( buffer[0] == '#' ) ReadLine(fp,bufferSize,buffer);  // skip comments
     
    sscanf(buffer,"%d %d",&width,&height);
     
    ReadLine(fp,bufferSize,buffer);
    while ( buffer[0] == '#' ) ReadLine(fp,bufferSize,buffer);  // skip comments
 
    // last read line should be "255\n"
 
    data.resize(width*height);
    fread( data.data(), sizeof(Color24), width*height, fp );
 
    return true;
}
 
//-------------------------------------------------------------------------------
 
bool TextureFile::Load()
{
    data.clear();
    width = 0;
    height = 0;
    FILE *fp = fopen( GetName(), "rb" );
    if ( ! fp ) return false;
 
    bool success = false;
    success = LoadPPM(fp,width,height,data);
 
    fclose(fp);
    return success;
}
 
//-------------------------------------------------------------------------------
 
Color TextureFile::Sample(const Point3 &uvw) const
{
    if ( width + height == 0 ) return Color(0,0,0);
 
    Point3 u = TileClamp(uvw);
    float x = width * u.x;
    float y = height * u.y;
    int ix = (int)x;
    int iy = (int)y;
    float fx = x - ix;
    float fy = y - iy;
 
    if ( ix < 0 ) ix -= (ix/width - 1)*width;
    if ( ix >= width ) ix -= (ix/width)*width;
    int ixp = ix+1;
    if ( ixp >= width ) ixp -= width;
 
    if ( iy < 0 ) iy -= (iy/height - 1)*height;
    if ( iy >= height ) iy -= (iy/height)*height;
    int iyp = iy+1;
    if ( iyp >= height ) iyp -= height;
 
    return  data[iy *width+ix ].ToColor() * ((1-fx)*(1-fy)) +
            data[iy *width+ixp].ToColor() * (   fx *(1-fy)) +
            data[iyp*width+ix ].ToColor() * ((1-fx)*   fy ) +
            data[iyp*width+ixp].ToColor() * (   fx *   fy );
}
 
//-------------------------------------------------------------------------------
 
bool TextureFile::SetViewportTexture() const
{
    if ( viewportTextureID == 0 ) {
        glGenTextures(1,&viewportTextureID);
        glBindTexture(GL_TEXTURE_2D,viewportTextureID);
        gluBuild2DMipmaps( GL_TEXTURE_2D, 3, width, height, GL_RGB, GL_UNSIGNED_BYTE, &data[0].r );
        glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR );
        glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
        glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT );
        glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT );
    }
    glBindTexture(GL_TEXTURE_2D,viewportTextureID);
    return true;
}
 
//-------------------------------------------------------------------------------
 
Color TextureChecker::Sample(const Point3 &uvw) const
{
    Point3 u = TileClamp(uvw);
    if ( u.x <= 0.5f ) {
        return u.y <= 0.5f ? color1 : color2;
    } else {
        return u.y <= 0.5f ? color2 : color1;
    }
}
 
//-------------------------------------------------------------------------------
 
bool TextureChecker::SetViewportTexture() const
{
    if ( viewportTextureID == 0 ) {
        const int texSize = 256;
        glGenTextures(1,&viewportTextureID);
        glBindTexture(GL_TEXTURE_2D,viewportTextureID);
        Color24 c[2] = { color1, color2 };
        Color24 *tex = new Color24[texSize*texSize];
        for ( int i=0; i<texSize*texSize; i++ ) {
            int ix = (i%texSize) < 128 ? 0 : 1;
            if ( i/256 >= 128 ) ix = 1 - ix;
            tex[i] = c[ix];
        }
        gluBuild2DMipmaps( GL_TEXTURE_2D, 3, texSize, texSize, GL_RGB, GL_UNSIGNED_BYTE, &tex[0].r );
        delete [] tex;
        glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR );
        glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
        glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT );
        glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT );
    }
    glBindTexture(GL_TEXTURE_2D,viewportTextureID);
    return true;
}
 
//-------------------------------------------------------------------------------