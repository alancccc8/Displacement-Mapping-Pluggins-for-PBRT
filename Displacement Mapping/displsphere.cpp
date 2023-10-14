/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.
*/
/*
    line 12 - 314 : Most codes are from pbrt-v3 source code.
    line 315 - 384 : Load OBJ data.
    line 385 - 428 : Load PGM data.
    line 429 - 512 : Updata vertices and normals, create new triangle meshes.
*/

// shapes/sphere.cpp*
#include "shapes/displSphere.h"
#include "sampling.h"
#include "paramset.h"
#include "efloat.h"
#include "stats.h"
#include "shapes/plymesh.h"
// load obj file
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <tuple>
// create new mesh
#include <iostream>
#include <cstdint>
#include "shapes/triangle.h"
#include "core/geometry.h"
#include "textures/constant.h"
#include "paramset.h"

namespace pbrt {

// Sphere Method Definitions
    Bounds3f DisplSphere::ObjectBound() const {
        return Bounds3f(Point3f(-radius, -radius, zMin),
                        Point3f(radius, radius, zMax));
    }

    bool DisplSphere::Intersect(const Ray &r, Float *tHit, SurfaceInteraction *isect,
                           bool testAlphaTexture) const {
        ProfilePhase p(Prof::ShapeIntersect);
        Float phi;
        Point3f pHit;
        // Transform _Ray_ to object space
        Vector3f oErr, dErr;
        Ray ray = (*WorldToObject)(r, &oErr, &dErr);

        // Compute quadratic sphere coefficients

        // Initialize _EFloat_ ray coordinate values
        EFloat ox(ray.o.x, oErr.x), oy(ray.o.y, oErr.y), oz(ray.o.z, oErr.z);
        EFloat dx(ray.d.x, dErr.x), dy(ray.d.y, dErr.y), dz(ray.d.z, dErr.z);
        EFloat a = dx * dx + dy * dy + dz * dz;
        EFloat b = 2 * (dx * ox + dy * oy + dz * oz);
        EFloat c = ox * ox + oy * oy + oz * oz - EFloat(radius) * EFloat(radius);

        // Solve quadratic equation for _t_ values
        EFloat t0, t1;
        if (!Quadratic(a, b, c, &t0, &t1)) return false;

        // Check quadric shape _t0_ and _t1_ for nearest intersection
        if (t0.UpperBound() > ray.tMax || t1.LowerBound() <= 0) return false;
        EFloat tShapeHit = t0;
        if (tShapeHit.LowerBound() <= 0) {
            tShapeHit = t1;
            if (tShapeHit.UpperBound() > ray.tMax) return false;
        }

        // Compute sphere hit position and $\phi$
        pHit = ray((Float)tShapeHit);

        // Refine sphere intersection point
        pHit *= radius / Distance(pHit, Point3f(0, 0, 0));
        if (pHit.x == 0 && pHit.y == 0) pHit.x = 1e-5f * radius;
        phi = std::atan2(pHit.y, pHit.x);
        if (phi < 0) phi += 2 * Pi;

        // Test sphere intersection against clipping parameters
        if ((zMin > -radius && pHit.z < zMin) || (zMax < radius && pHit.z > zMax) ||
            phi > phiMax) {
            if (tShapeHit == t1) return false;
            if (t1.UpperBound() > ray.tMax) return false;
            tShapeHit = t1;
            // Compute sphere hit position and $\phi$
            pHit = ray((Float)tShapeHit);

            // Refine sphere intersection point
            pHit *= radius / Distance(pHit, Point3f(0, 0, 0));
            if (pHit.x == 0 && pHit.y == 0) pHit.x = 1e-5f * radius;
            phi = std::atan2(pHit.y, pHit.x);
            if (phi < 0) phi += 2 * Pi;
            if ((zMin > -radius && pHit.z < zMin) ||
                (zMax < radius && pHit.z > zMax) || phi > phiMax)
                return false;
        }

        // Find parametric representation of sphere hit
        Float u = phi / phiMax;
        Float theta = std::acos(Clamp(pHit.z / radius, -1, 1));
        Float v = (theta - thetaMin) / (thetaMax - thetaMin);

        // Compute sphere $\dpdu$ and $\dpdv$
        Float zRadius = std::sqrt(pHit.x * pHit.x + pHit.y * pHit.y);
        Float invZRadius = 1 / zRadius;
        Float cosPhi = pHit.x * invZRadius;
        Float sinPhi = pHit.y * invZRadius;
        Vector3f dpdu(-phiMax * pHit.y, phiMax * pHit.x, 0);
        Vector3f dpdv =
            (thetaMax - thetaMin) *
            Vector3f(pHit.z * cosPhi, pHit.z * sinPhi, -radius * std::sin(theta));

        // Compute sphere $\dndu$ and $\dndv$
        Vector3f d2Pduu = -phiMax * phiMax * Vector3f(pHit.x, pHit.y, 0);
        Vector3f d2Pduv =
            (thetaMax - thetaMin) * pHit.z * phiMax * Vector3f(-sinPhi, cosPhi, 0.);
        Vector3f d2Pdvv = -(thetaMax - thetaMin) * (thetaMax - thetaMin) *
                          Vector3f(pHit.x, pHit.y, pHit.z);

        // Compute coefficients for fundamental forms
        Float E = Dot(dpdu, dpdu);
        Float F = Dot(dpdu, dpdv);
        Float G = Dot(dpdv, dpdv);
        Vector3f N = Normalize(Cross(dpdu, dpdv));
        Float e = Dot(N, d2Pduu);
        Float f = Dot(N, d2Pduv);
        Float g = Dot(N, d2Pdvv);

        // Compute $\dndu$ and $\dndv$ from fundamental form coefficients
        Float invEGF2 = 1 / (E * G - F * F);
        Normal3f dndu = Normal3f((f * F - e * G) * invEGF2 * dpdu +
                                 (e * F - f * E) * invEGF2 * dpdv);
        Normal3f dndv = Normal3f((g * F - f * G) * invEGF2 * dpdu +
                                 (f * F - g * E) * invEGF2 * dpdv);

        // Compute error bounds for sphere intersection
        Vector3f pError = gamma(5) * Abs((Vector3f)pHit);

        // Initialize _SurfaceInteraction_ from parametric information
        *isect = (*ObjectToWorld)(SurfaceInteraction(pHit, pError, Point2f(u, v),
                                                     -ray.d, dpdu, dpdv, dndu, dndv,
                                                     ray.time, this));

        // Update _tHit_ for quadric intersection
        *tHit = (Float)tShapeHit;
        return true;
    }

    bool DisplSphere::IntersectP(const Ray &r, bool testAlphaTexture) const {
        ProfilePhase p(Prof::ShapeIntersectP);
        Float phi;
        Point3f pHit;
        // Transform _Ray_ to object space
        Vector3f oErr, dErr;
        Ray ray = (*WorldToObject)(r, &oErr, &dErr);

        // Compute quadratic sphere coefficients

        // Initialize _EFloat_ ray coordinate values
        EFloat ox(ray.o.x, oErr.x), oy(ray.o.y, oErr.y), oz(ray.o.z, oErr.z);
        EFloat dx(ray.d.x, dErr.x), dy(ray.d.y, dErr.y), dz(ray.d.z, dErr.z);
        EFloat a = dx * dx + dy * dy + dz * dz;
        EFloat b = 2 * (dx * ox + dy * oy + dz * oz);
        EFloat c = ox * ox + oy * oy + oz * oz - EFloat(radius) * EFloat(radius);

        // Solve quadratic equation for _t_ values
        EFloat t0, t1;
        if (!Quadratic(a, b, c, &t0, &t1)) return false;

        // Check quadric shape _t0_ and _t1_ for nearest intersection
        if (t0.UpperBound() > ray.tMax || t1.LowerBound() <= 0) return false;
        EFloat tShapeHit = t0;
        if (tShapeHit.LowerBound() <= 0) {
            tShapeHit = t1;
            if (tShapeHit.UpperBound() > ray.tMax) return false;
        }

        // Compute sphere hit position and $\phi$
        pHit = ray((Float)tShapeHit);

        // Refine sphere intersection point
        pHit *= radius / Distance(pHit, Point3f(0, 0, 0));
        if (pHit.x == 0 && pHit.y == 0) pHit.x = 1e-5f * radius;
        phi = std::atan2(pHit.y, pHit.x);
        if (phi < 0) phi += 2 * Pi;

        // Test sphere intersection against clipping parameters
        if ((zMin > -radius && pHit.z < zMin) || (zMax < radius && pHit.z > zMax) ||
            phi > phiMax) {
            if (tShapeHit == t1) return false;
            if (t1.UpperBound() > ray.tMax) return false;
            tShapeHit = t1;
            // Compute sphere hit position and $\phi$
            pHit = ray((Float)tShapeHit);

            // Refine sphere intersection point
            pHit *= radius / Distance(pHit, Point3f(0, 0, 0));
            if (pHit.x == 0 && pHit.y == 0) pHit.x = 1e-5f * radius;
            phi = std::atan2(pHit.y, pHit.x);
            if (phi < 0) phi += 2 * Pi;
            if ((zMin > -radius && pHit.z < zMin) ||
                (zMax < radius && pHit.z > zMax) || phi > phiMax)
                return false;
        }
        return true;
    }

    Float DisplSphere::Area() const { return phiMax * radius * (zMax - zMin); }

    Interaction DisplSphere::Sample(const Point2f &u, Float *pdf) const {
        Point3f pObj = Point3f(0, 0, 0) + radius * UniformSampleSphere(u);
        Interaction it;
        it.n = Normalize((*ObjectToWorld)(Normal3f(pObj.x, pObj.y, pObj.z)));
        if (reverseOrientation) it.n *= -1;
        // Reproject _pObj_ to sphere surface and compute _pObjError_
        pObj *= radius / Distance(pObj, Point3f(0, 0, 0));
        Vector3f pObjError = gamma(5) * Abs((Vector3f)pObj);
        it.p = (*ObjectToWorld)(pObj, pObjError, &it.pError);
        *pdf = 1 / Area();
        return it;
    }

    Interaction DisplSphere::Sample(const Interaction &ref, const Point2f &u,
                               Float *pdf) const {
        Point3f pCenter = (*ObjectToWorld)(Point3f(0, 0, 0));

        // Sample uniformly on sphere if $\pt{}$ is inside it
        Point3f pOrigin =
            OffsetRayOrigin(ref.p, ref.pError, ref.n, pCenter - ref.p);
        if (DistanceSquared(pOrigin, pCenter) <= radius * radius) {
            Interaction intr = Sample(u, pdf);
            Vector3f wi = intr.p - ref.p;
            if (wi.LengthSquared() == 0)
                *pdf = 0;
            else {
                // Convert from area measure returned by Sample() call above to
                // solid angle measure.
                wi = Normalize(wi);
                *pdf *= DistanceSquared(ref.p, intr.p) / AbsDot(intr.n, -wi);
            }
            if (std::isinf(*pdf)) *pdf = 0.f;
            return intr;
        }

        // Sample sphere uniformly inside subtended cone

        // Compute coordinate system for sphere sampling
        Float dc = Distance(ref.p, pCenter);
        Float invDc = 1 / dc;
        Vector3f wc = (pCenter - ref.p) * invDc;
        Vector3f wcX, wcY;
        CoordinateSystem(wc, &wcX, &wcY);

        // Compute $\theta$ and $\phi$ values for sample in cone
        Float sinThetaMax = radius * invDc;
        Float sinThetaMax2 = sinThetaMax * sinThetaMax;
        Float invSinThetaMax = 1 / sinThetaMax;
        Float cosThetaMax = std::sqrt(std::max((Float)0.f, 1 - sinThetaMax2));

        Float cosTheta  = (cosThetaMax - 1) * u[0] + 1;
        Float sinTheta2 = 1 - cosTheta * cosTheta;

        if (sinThetaMax2 < 0.00068523f /* sin^2(1.5 deg) */) {
            /* Fall back to a Taylor series expansion for small angles, where
               the standard approach suffers from severe cancellation errors */
            sinTheta2 = sinThetaMax2 * u[0];
            cosTheta = std::sqrt(1 - sinTheta2);
        }

        // Compute angle $\alpha$ from center of sphere to sampled point on surface
        Float cosAlpha = sinTheta2 * invSinThetaMax +
            cosTheta * std::sqrt(std::max((Float)0.f, 1.f - sinTheta2 * invSinThetaMax * invSinThetaMax));
        Float sinAlpha = std::sqrt(std::max((Float)0.f, 1.f - cosAlpha*cosAlpha));
        Float phi = u[1] * 2 * Pi;

        // Compute surface normal and sampled point on sphere
        Vector3f nWorld =
            SphericalDirection(sinAlpha, cosAlpha, phi, -wcX, -wcY, -wc);
        Point3f pWorld = pCenter + radius * Point3f(nWorld.x, nWorld.y, nWorld.z);

        // Return _Interaction_ for sampled point on sphere
        Interaction it;
        it.p = pWorld;
        it.pError = gamma(5) * Abs((Vector3f)pWorld);
        it.n = Normal3f(nWorld);
        if (reverseOrientation) it.n *= -1;

        // Uniform cone PDF.
        *pdf = 1 / (2 * Pi * (1 - cosThetaMax));

        return it;
    }

    Float DisplSphere::Pdf(const Interaction &ref, const Vector3f &wi) const {
        Point3f pCenter = (*ObjectToWorld)(Point3f(0, 0, 0));
        // Return uniform PDF if point is inside sphere
        Point3f pOrigin =
            OffsetRayOrigin(ref.p, ref.pError, ref.n, pCenter - ref.p);
        if (DistanceSquared(pOrigin, pCenter) <= radius * radius)
            return Shape::Pdf(ref, wi);

        // Compute general sphere PDF
        Float sinThetaMax2 = radius * radius / DistanceSquared(ref.p, pCenter);
        Float cosThetaMax = std::sqrt(std::max((Float)0, 1 - sinThetaMax2));
        return UniformConePdf(cosThetaMax);
    }

    Float DisplSphere::SolidAngle(const Point3f &p, int nSamples) const {
        Point3f pCenter = (*ObjectToWorld)(Point3f(0, 0, 0));
        if (DistanceSquared(p, pCenter) <= radius * radius)
            return 4 * Pi;
        Float sinTheta2 = radius * radius / DistanceSquared(p, pCenter);
        Float cosTheta = std::sqrt(std::max((Float)0, 1 - sinTheta2));
        return (2 * Pi * (1 - cosTheta));
    }


    //define a function to load obj. file and acquire the vertex, face information

    /*
        reference: http://www.opengl-tutorial.org/beginners-tutorials/tutorial-7-model-loading/
        reference: https://blog.csdn.net/qq_43703185/article/details/119143768
        reference: https://blog.csdn.net/weixin_44176696/article/details/110884072
        reference: https://blog.csdn.net/Unclings/article/details/123708411
        reference: shapes\plymesh.cpp
    */
    void loadObjData(const std::string& filename,
                     std::vector<Point3f>& vertices,
                     std::vector<Normal3f>& normals,
                     std::vector<Point2f>& uvs,
                     std::vector<int>& vertexIndices,
                     std::vector<int>& uvIndices,
                     std::vector<int>& normalIndices
                     ){
        // Object: objFile, Usage: read data from obj file
        std::ifstream objFile(filename);
        // error handling: cannot open file
        if (!objFile.is_open()) {
            std::cout << "Failed to open file: " << filename << std::endl;
            return;
        }
        std::string line;
        // Read data line by line, store each line in a line variable until the end of the objFile.
        while (std::getline(objFile, line)){
            // Read a string from 'line' and store it in 'prefix' 
            std::istringstream lineStream(line);
            std::string prefix;
            lineStream >> prefix;
            // Determine information types by prefix
            if (prefix == "v"){
                // vertex XYZ Coord -> vertices
                Point3f vertex;
                // error handling: The data is not in the correct format or has been read to the end of the file
                /*if (!(lineStream >> vertex.x >> vertex.y >> vertex.z)) {
                    std::cerr << "Error reading vertex data at line " << lineNumber << ": " << line << std::endl;
                    break;
                }*/
                lineStream >> vertex.x >> vertex.y >> vertex.z;
                vertices.emplace_back(vertex);
            } else if(prefix == "vn"){
                // vertex normal -> normals
                Normal3f normal;
                lineStream >> normal.x >> normal.y >> normal.z;
                normals.emplace_back(normal);
            }  else if (prefix == "vt") {
                // vertex UV coord -> uvs
                Point2f uv;
                lineStream >> uv.x >> uv.y;
                uvs.push_back(uv);
            }  else if (prefix == "f") {
                // Each row has three sets of data, in oder of vetexIndex, uvIndex, normalIndex; format: 1/1/1 2/2/2 3/3/3
                int vertexIndex, uvIndex, normalIndex;
                char slash;
                for (int i = 0; i < 3; ++i) {
                    lineStream >> vertexIndex >> slash >> uvIndex >> slash >> normalIndex;
                    // store all indices; the index of the array starts from 0
                    vertexIndices.emplace_back(vertexIndex-1);
                    uvIndices.emplace_back(uvIndex-1);
                    normalIndices.emplace_back(normalIndex-1);
                }
            }
        }
        objFile.close();
    }

    //define a function to load pgm. file and acquire the gray scale 
    void loadPgmData(const std::string& filename, 
                     int& width, int& height, 
                     std::vector<int>& grayValues) {
        // Object: pgmFile, Usage: read gray scale data from pgm file
        std::ifstream pgmFile(filename);
        if (!pgmFile.is_open()){
            std::cout << "Failed to open file: " << filename << std::endl;
            return;
        }
        std::string line;
        // read the header on the first line
        pgmFile >> line;
        // error handling: incorrect pgm format
        if (line != "P2" && line != "P5") {
            std::cout << "Invalid PGM format." << std::endl;
            return;
        }
        // read the secon line: width and height
        pgmFile >> width >> height;
        // read the maxvalue
        int maxValue;
        pgmFile >> maxValue;
        // Get and ignore line breaks between lines
        //pgmFile.get();
        // Resize 'grayValues' vector to store all the pixel values.
        grayValues.resize(width * height);
        // Get grayValues according to the pgm file type
        if (line == "P2") {
            // ASCII PGM
            for (int i = 0; i < width * height; ++i) {
                int value;
                pgmFile >> value;
                grayValues[i] = static_cast<int>(value);
            }
        } else {
            // Binary PGM
            for (int i = 0; i < width * height; ++i) {  
                grayValues[i] = static_cast<int>(pgmFile.get());
            }
        }
        pgmFile.close();
    }


    // According to the ObJfile and PGMfile, create the new shape
    std::vector <std::shared_ptr<Shape>> CreateDisplSphereShape(const Transform *o2w,
                                             const Transform *w2o,
                                             bool reverseOrientation,
                                             const ParamSet &params) {
        Float radius = params.FindOneFloat("radius", 1.f);
        Float zmin = params.FindOneFloat("zmin", -radius);
        Float zmax = params.FindOneFloat("zmax", radius);
        Float phimax = params.FindOneFloat("phimax", 360.f);

        std::string defaultObjFilePath = "C://Users/15469/source/ball.obj";
        std::string defaultPgmFilePath = "C://Users/15469/source/texture.pgm";
        std::string objFile = params.FindOneFilename("objfile", defaultObjFilePath);
        std::string pgmFile = params.FindOneFilename("pgmfile", defaultPgmFilePath);

        // Customised displacement scaling factor
        Float displacementScale = params.FindOneFloat("displacementScale", 0.5f);
        // Define new objects for loadObjData()
        std::vector<Point3f> vertices;
        std::vector<Normal3f> normals;
        std::vector<Point2f> uvs;
        std::vector<int> vertexIndices;
        std::vector<int> uvIndices;
        std::vector<int> normalIndices;
        loadObjData(objFile, vertices, normals, uvs,
                    vertexIndices, uvIndices,normalIndices);
        // Define new objects for loadPgmData()
        int width;
        int height; 
        std::vector<int> grayValues;
        loadPgmData(pgmFile, width, height, grayValues);
        // Track vertices have been updated, avoid repeated updates of shared vertices
        std::vector<bool> processed(vertices.size(), false);
        // Iterate over all original triangle surface, get its three vertices, generate new vertices
        for (size_t i = 0; i < vertexIndices.size(); i+=3){
        // Vertex
            for (int j = 0; j < 3; ++j) {
                int vertexIndex = vertexIndices[i+j];; 
                int uvIndex = uvIndices[i+j];
                // Check if vertex has been processed
                if (!processed[vertexIndex]){
                    // If unprocessed, get the new vertex index
                    Point2f uv = uvs[uvIndex];
                    // Map uv coord to the pgm image to get grayValue
                    int x = static_cast<int>(uv.x * (width - 1));
                    int y = static_cast<int>(uv.y * (height - 1));
                    int grayValue = grayValues[y * width + x];
                    // Calculate the displacement value
                    Float displacement = (grayValue / 255.0f) * displacementScale;
                    // Displace vertex along the normal
                    int normalIndex = normalIndices[i+j];
                    Normal3f normal = normals[normalIndex];
                    Point3f normal3f(normal.x, normal.y, normal.z);
                    vertices[vertexIndex] += normal3f * displacement;
                    // Update the vertex processed status
                    processed[vertexIndex] = true;  
                }
            }
        }
        // Update the normals, proecess one triangle at a time
        for (size_t i = 0; i < vertexIndices.size(); i += 3) {
            // get the three vertices of the triangle
            Point3f v0 = vertices[vertexIndices[i]];
            Point3f v1 = vertices[vertexIndices[i + 1]];
            Point3f v2 = vertices[vertexIndices[i + 2]];
            // calculate the normal of the triangle
            Normal3f faceNormal = Normal3f(Normalize(Cross(v1 - v0, v2 - v0)));
            // update the normal of each vertex, determinded by adjacent Triangles
            for (int j = 0; j < 3; ++j) {
                int normalIndex = normalIndices[i + j];
                // Accumulating adjacent triangle normals
                normals[normalIndex] += faceNormal;
            }
        }
        // Normalize the new normals
        for (Normal3f &normal : normals) {
            normal = Normalize(normal);
        }
        // Create triangle mesh
        return CreateTriangleMesh(o2w, w2o, reverseOrientation,
                                vertexIndices.size() / 3, vertexIndices.data(),
                                vertices.size(), vertices.data(), 
                                nullptr, normals.data(), 
                                uvs.data(),nullptr, nullptr,
                                nullptr);    
    } 
} // namespace pbrt

