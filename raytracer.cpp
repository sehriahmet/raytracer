#include <iostream>
#include "parser.h"
#include "ppm.h"

#include <math.h>
#include <cmath>
#include <algorithm>

parser::Scene scene;

using vec3f = parser::Vec3f;
using vec3i = parser::Vec3i;
using vec4f = parser::Vec4f;
using PointLight = parser::PointLight;
using Material = parser::Material;
using Sphere = parser::Sphere;
using Face = parser::Face;
using Triangle = parser::Triangle;
using Mesh = parser::Mesh;

// Declare functions before usage
vec3f multiplicationScalarf(vec3f a, float s);
vec3f addVectorsf(vec3f a, vec3f b);
vec3f cross(vec3f a, vec3f b);
vec3f substractVectorsf(vec3f a, vec3f b);

vec3f e;
vec3f upVectoru, upVectorv, upVectorw;

typedef struct {
    vec3f origin, direction;
} ray;

vec3f computeColor(ray r, int depth);

// burada ilk camera elementini aldim ama birden fazl aoldugu case'i sonradan yapmak lazim
// Main function or initialization section

vec4f nearPlane; // left, right, bottom, top
vec3f gaze;
float dist;
int width;
int height;
int nx;
int ny;
vec3i backgroundColori;
vec3f backgroundColor;

void initializeCameraVectors() {
    // Here, we assume the scene has already been parsed
    const auto& camera = scene.cameras[0];

    // initialize background color also
    backgroundColori = scene.background_color;
    backgroundColor.x = float(backgroundColori.x);
    backgroundColor.y = float(backgroundColori.y);
    backgroundColor.z = float(backgroundColori.z);

    upVectorv = camera.up;
    upVectorw = multiplicationScalarf(camera.gaze, -1);
    upVectoru = cross(upVectorv, upVectorw);
    
    nearPlane = camera.near_plane;
    gaze = camera.gaze;
    dist = camera.near_distance;
    width = camera.image_width;
    height = camera.image_height;
    e = camera.position;

    nx = camera.image_width;
    ny = camera.image_height;

}

vec3f multiplicationScalarf(vec3f a, float s) {
    vec3f result;

    result.x = a.x * s;
    result.y = a.y * s;
    result.z = a.z * s;

    return result;
}

vec3f addVectorsf (vec3f a, vec3f b) {
    vec3f result;

    result.x = a.x + b.x;
    result.y = a.y + b.y;
    result.z = a.z + b.z;

    return result;
}

vec3f substractVectorsf (vec3f a, vec3f b) {
    vec3f result;
    
    result.x = a.x - b.x;
    result.y = a.y - b.y;
    result.z = a.z - b.z;

    return result;
}

vec3f cross(vec3f a, vec3f b) {
    vec3f result;

    result.x = a.y * b.z - a.z * b.y;
    result.y = a.z * b.x - a.x * b.z;
    result.z = a.x * b.y - a.y * b.x;

    return result;
}

float dot (vec3f a, vec3f b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

vec3f normalize (vec3f v) {
    vec3f result;
    
    float length;
    length = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);

    result.x = v.x/length;
    result.y = v.y/length;
    result.z = v.z/length;

    return result;
}

float length(vec3f v) {
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

vec3f reflect (vec3f light, vec3f normal) {
    return substractVectorsf(multiplicationScalarf(normal, 2 * dot(light, normal)), light);
}

ray generateRay (int i, int j) {
    ray result;

    float su, sv;

    su = (i+0.5) * (nearPlane.y-nearPlane.x)/width; // right-left
    sv = (j+0.5) * (nearPlane.w-nearPlane.z)/height; // top-bottom

    vec3f m, q, s;

    // m ve q generateRay disina alinabilir.
    m = addVectorsf(e, multiplicationScalarf(gaze, dist)); // middle point of the near plane 
    q = addVectorsf(m, addVectorsf(multiplicationScalarf(upVectoru, nearPlane.x), multiplicationScalarf(upVectorv, nearPlane.w))); // q = m + lu + tv // top-left corner
    s = addVectorsf(q, addVectorsf(multiplicationScalarf(upVectoru, su), multiplicationScalarf(upVectorv,-sv))); // pixel position

    result.origin = e;
    result.direction = substractVectorsf(s, e);

    return result;
}

float intersectSphere (ray r, Sphere s) {
    float A, B, C; // constants for quadratic equations

    float delta;

    vec3f c;

    // burada sphere vertex id'den cikarip bulmak gerekiyo
    c = scene.vertex_data[s.center_vertex_id-1];

    // printf("scene vertex datadan alinanlar: %d %f %f %f \n", s.center_vertex_id, c.x, c.y ,c.z);
    float t, t1, t2;

    C = (r.origin.x - c.x)*(r.origin.x- c.x) + (r.origin.y - c.y)*(r.origin.y- c.y) + (r.origin.z - c.z)*(r.origin.z- c.z) - s.radius*s.radius;
    B = 2*r.direction.x*(r.origin.x - c.x)+ 2*r.direction.y*(r.origin.y - c.y) + 2*r.direction.z*(r.origin.z - c.z);
    A = r.direction.x * r.direction.x + r.direction.y * r.direction.y + r.direction.z * r.direction.z;

    delta = B*B-4*A*C;

    if (delta<0) return -1;

    else if (delta == 0) {
        t = -B / (2*A);
    }

    else {
        delta = sqrt(delta);
        A = 2*A;
        t1 = (-B + delta) / A;
        t2 = (-B - delta) / A;

        if (t1<t2) t=t1; else t=t2;
    }
    // if ( t> 0) printf("intersectSphere: %f \n", t);
    return t;
}

// Function to check intersection of a ray with a triangle using barycentric coordinates
float intersectTriangle(ray r, Triangle intersectedTriangle) {
    // Retrieve triangle vertices from the scene's vertex data
    vec3f v0 = scene.vertex_data[intersectedTriangle.indices.v0_id - 1];
    vec3f v1 = scene.vertex_data[intersectedTriangle.indices.v1_id - 1];
    vec3f v2 = scene.vertex_data[intersectedTriangle.indices.v2_id - 1];

    // Calculate the edges of the triangle
    vec3f edge1 = substractVectorsf(v1, v0);
    vec3f edge2 = substractVectorsf(v2, v0);

    // Calculate the normal and determinant for the plane intersection
    vec3f h = cross(r.direction, edge2);
    float a = dot(edge1, h);

    const float EPSILON = 1e-5;
    if (fabs(a) < EPSILON) return -1; // Ray is parallel to the triangle

    // Calculate f, s, u, v
    float f = 1.0 / a;
    vec3f s = substractVectorsf(r.origin, v0);
    float u = f * dot(s, h);
    if (u < 0.0 || u > 1.0) return -1; // Intersection outside the triangle

    vec3f q = cross(s, edge1);
    float v = f * dot(r.direction, q);
    if (v < 0.0 || u + v > 1.0) return -1; // Intersection outside the triangle

    // Calculate t to determine the intersection point distance
    float t = f * dot(edge2, q);

    // Valid intersection if t is positive
    return (t > EPSILON) ? t : -1;
}

// mesh icin intersect fonksiyonu triangle'lari kullanarak bulunuyor
float intersectMesh(ray r, Mesh mesh) {
    float closestT = -1;
    for (int i = 0; i < mesh.faces.size(); ++i) {
        Face indices = mesh.faces[i];

        Triangle tri;
        tri.material_id = mesh.material_id;
        tri.indices = indices;

        float t = intersectTriangle(r, tri);
        if (t > 0 && (closestT < 0 || t < closestT)) {
            closestT = t;
        }
    }

    return closestT;
}


// TO-DO: in case of some bugs, this shadow function should be checked because it is setting 0 to color if in the shadow 
// shadow part
bool isInShadow(vec3f intersectionPoint, vec3f lightPosition) {
    vec3f shadowRay = substractVectorsf(lightPosition, intersectionPoint);
    float distanceToLight = length(shadowRay);
    shadowRay = normalize(shadowRay);

    ray shadowRayStruct;
    shadowRayStruct.origin = intersectionPoint;
    shadowRayStruct.direction = shadowRay;

    for (int i = 0; i < scene.spheres.size(); i++) {
        float t = intersectSphere(shadowRayStruct, scene.spheres[i]);
        if (t > 0 && t < distanceToLight) {
            return true; 
        }
    }

    for (int i = 0; i < scene.triangles.size(); i++) {
        float t = intersectTriangle(shadowRayStruct, scene.triangles[i]);
        if (t > 0 && t < distanceToLight) {
            return true; 
        }
    }

    for (int i = 0; i < scene.meshes.size(); i++) {
        float t = intersectMesh(shadowRayStruct, scene.meshes[i]);
        if (t > 0 && t < distanceToLight) {
            return true; 
        }
    }

    return false;
}

// Phong shading function
vec3f calculateColor(int materialId, vec3f intersectionPoint, vec3f normal, ray myRay, int depth) {
    vec3f color = {0, 0, 0};
    
    Material material = scene.materials[materialId-1];

    vec3f ambientLight = scene.ambient_light;
    // ambient component
    color.x += material.ambient.x * ambientLight.x / 255;
    color.y += material.ambient.y * ambientLight.y / 255;
    color.z += material.ambient.z * ambientLight.z / 255;

    // printf("ambientlight.x= %.2f\n", ambientLight.x);
    PointLight light;
    for (int i = 0; i < scene.point_lights.size(); i++){

        light = scene.point_lights[i];
        if (isInShadow(intersectionPoint, light.position)) {
            continue;
        }

        // Diffuse component
        vec3f L = normalize(substractVectorsf(light.position, intersectionPoint)); // Incoming light direction
        float distance = length(substractVectorsf(light.position, intersectionPoint)); // Distance to light
        float cosThetaPrime = fmax(0.0, dot(L, normal)); // Cosine of the angle between light direction and normal
        // float diffuse = (material.diffuse.x * cosThetaPrime * light.intensity.x) / (distance * distance) / 255;
        vec3f received_irradience;
        received_irradience.x = light.intensity.x / (distance * distance);
        received_irradience.y = light.intensity.y / (distance * distance);
        received_irradience.z = light.intensity.z / (distance * distance);

        color.x += material.diffuse.x * cosThetaPrime * received_irradience.x / 255;
        color.y += material.diffuse.y * cosThetaPrime * received_irradience.y / 255;
        color.z += material.diffuse.z * cosThetaPrime * received_irradience.z / 255;

        // Specular component
        vec3f w0 = normalize(substractVectorsf(scene.cameras[0].position, intersectionPoint));
        vec3f h = normalize(addVectorsf(L, w0));
        float cosAlphaPrime = fmax(0.0, dot(normal, h));

        vec3f reflectDir = substractVectorsf(multiplicationScalarf(normal, 2 * dot(L, normal)), L);

        float specular = powf(cosAlphaPrime, material.phong_exponent);
        color.x += material.specular.x * specular * received_irradience.x / 255;
        color.y += material.specular.y * specular * received_irradience.y / 255;
        color.z += material.specular.z * specular * received_irradience.z / 255;

    }

    // reflect ray baslangic konumu ilk nokta

    if (material.is_mirror) {
        vec3f reflectionDir = substractVectorsf(myRay.direction, multiplicationScalarf(normal, 2 * dot(myRay.direction, normal)));
        reflectionDir = normalize(reflectionDir);

        ray reflectionRay;
        reflectionRay.origin = addVectorsf(intersectionPoint, multiplicationScalarf(reflectionDir, scene.shadow_ray_epsilon));
        reflectionRay.direction = reflectionDir;

        
        vec3f reflectedColor = computeColor(reflectionRay, depth-1);

        if (reflectedColor.x <= 1 || reflectedColor.y <= 1 || reflectedColor.z <= 1){
            color.x += material.mirror.x * reflectedColor.x;
            color.y += material.mirror.y * reflectedColor.y;
            color.z += material.mirror.z * reflectedColor.z;
        }
    }

    // range is in 0, 1
    color.x = std::min(1.0f, std::max(0.0f, color.x));
    color.y = std::min(1.0f, std::max(0.0f, color.y));
    color.z = std::min(1.0f, std::max(0.0f, color.z));
    
    return color;
}

int minSphereI, minTriangleI, minMeshI;
vec3f computeColor (ray myRay, int depth) {

    if (depth<=0) return backgroundColor;

    int i;
    vec3f c;
    float minT_sphere = 90000; // some large number
    float minT_triangle = 90000;
    float minT_mesh = 90000;
    float t_sphere, t_triangle, t_mesh;
    vec3f L, N, P;

    c = backgroundColor;
    minSphereI = -1;
    minTriangleI = -1;
    minMeshI = -1;

    for (i=0; i<scene.spheres.size(); i++) {

        t_sphere = intersectSphere(myRay, scene.spheres[i]);

        if (t_sphere<minT_sphere && t_sphere>=0.001) {

            Sphere sphere = scene.spheres[i];

            minT_sphere = t_sphere;
            minSphereI = i;
        }
    }

    for (i = 0; i < scene.triangles.size(); i++) {
        t_triangle = intersectTriangle(myRay, scene.triangles[i]);

        if (t_triangle < minT_triangle && t_triangle >= 0.001) {
            minT_triangle = t_triangle;
            minTriangleI = i;
        }
    }

    for (i = 0; i < scene.meshes.size(); i++) {
        t_mesh = intersectMesh(myRay, scene.meshes[i]);

        if (t_mesh < minT_mesh && t_mesh >= 0.001) {
            minT_mesh = t_mesh;
            minMeshI = i;
        }
    }

    if (minSphereI != -1 && (minT_sphere < minT_triangle || minTriangleI == -1) && (minT_sphere < minT_mesh || minMeshI == -1)) {
        Sphere sphere = scene.spheres[minSphereI];
        vec3f sphereColor = {0, 0, 0};

        // printf("spherecolor.x = %f, .y = %f, .z = %f\n", sphereColor.x, sphereColor.y, sphereColor.z);

        // intersection point P
        vec3f P = addVectorsf(myRay.origin, multiplicationScalarf(myRay.direction, minT_sphere));

        // light direction vector L
        vec3f L = substractVectorsf(scene.point_lights[0].position, P);
        L = normalize(L);

        // normal vector N at the intersection point
        vec3f N = substractVectorsf(P, scene.vertex_data[sphere.center_vertex_id-1]);
        N = normalize(N);

        c = calculateColor(sphere.material_id, P, N, myRay, depth);
    }
    
    else if (minTriangleI != -1 && (minT_triangle < minT_mesh || minMeshI == -1)) {
        Triangle triangle = scene.triangles[minTriangleI];

        P = addVectorsf(myRay.origin, multiplicationScalarf(myRay.direction, minT_triangle));

        L = substractVectorsf(scene.point_lights[0].position, P);
        L = normalize(L);

        vec3f v0 = scene.vertex_data[triangle.indices.v0_id - 1];
        vec3f v1 = scene.vertex_data[triangle.indices.v1_id - 1];
        vec3f v2 = scene.vertex_data[triangle.indices.v2_id - 1];
        vec3f edge1 = substractVectorsf(v1, v0);
        vec3f edge2 = substractVectorsf(v2, v0);
        N = normalize(cross(edge1, edge2));

        c = calculateColor(triangle.material_id, P, N, myRay, depth);
    }

    else if (minMeshI != -1) {
        Mesh mesh = scene.meshes[minMeshI];

        P = addVectorsf(myRay.origin, multiplicationScalarf(myRay.direction, minT_mesh));

        L = substractVectorsf(scene.point_lights[0].position, P);
        L = normalize(L);

        Face closestFace = mesh.faces[minMeshI];
        vec3f v0 = scene.vertex_data[closestFace.v0_id - 1];
        vec3f v1 = scene.vertex_data[closestFace.v1_id - 1];
        vec3f v2 = scene.vertex_data[closestFace.v2_id - 1];
        vec3f edge1 = substractVectorsf(v1, v0);
        vec3f edge2 = substractVectorsf(v2, v0);
        N = normalize(cross(edge1, edge2));

        c = calculateColor(mesh.material_id, P, N, myRay, depth);
    }

    return c;
}


int main(int argc, char* argv[])
{
    scene.loadFromXml(argv[1]);
    // printf("sphere[0] material id: %d, material: %f\n", scene.spheres[0].material_id, scene.materials[scene.spheres[0].material_id-1].ambient.x);
    // printf("%f\n", scene.vertex_data[0].x);
    // printf("%zu \n", scene.spheres.size());
    initializeCameraVectors();

    unsigned char* image = new unsigned char [width * height * 3];
    
    int countPixel = 0; // for writing ppm in the image side

    // ray tracing loop
    for (int j = 0; j < height; ++j)
    {
        for (int i = 0; i < width; ++i)
        {
            // treat i as column value, j as row value 
            ray generatedRay = generateRay(i, j);
            
            vec3f pixel;
            pixel = addVectorsf(generatedRay.origin, generatedRay.direction);
            
            vec3f rayColor;
            rayColor = computeColor(generatedRay, scene.max_recursion_depth);
            // if (rayColor.x > 0.0) printf("raycolors: %f %f %f \n", rayColor.x, rayColor.y, rayColor.z);
            
            image[countPixel++] = (int)(rayColor.x*255+0.5);
            image[countPixel++] = (int)(rayColor.y*255+0.5);
            image[countPixel++] = (int)(rayColor.z*255+0.5);

            // printf("pixels: %.4f %.4f %.4f \n", pixel.x, pixel.y, pixel.z);

        }
    }

    write_ppm("test.ppm", image, width, height);

}
