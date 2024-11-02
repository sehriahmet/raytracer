#include <iostream>
#include "parser.h"
#include "ppm.h"

#include <math.h>
#include <cmath>
#include <algorithm>

#include <chrono> // for measuring time

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
using Camera = parser::Camera;

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
Camera camera;

void initializeCameraVectors(Camera camera) {
    // Here, we assume the scene has already been parsed
    // const auto& camera = scene.cameras[0];

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

    // TODO: epsilon degeri hata payini ayarliyo aslinda bazi taranamayan ucgenler bu yuzden gozukmuyor olabilir.
    const float EPSILON = 1e-7; // 1e-5 
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

// TO-DO: in case of some bugs, this shadow function should be checked because it is setting 0 to color if in the shadow 
// shadow part
bool isInShadow(vec3f intersectionPoint, vec3f lightPosition) {
    // Slightly offset the shadow ray to avoid self-intersection ("shadow acne")
    vec3f shadowRay = substractVectorsf(lightPosition, intersectionPoint);
    float distanceToLight = length(shadowRay);
    shadowRay = normalize(shadowRay);

    // Small epsilon to avoid self-intersection
    vec3f shadowRayOrigin = addVectorsf(intersectionPoint, multiplicationScalarf(shadowRay, scene.shadow_ray_epsilon));

    ray shadowRayStruct;
    shadowRayStruct.origin = shadowRayOrigin;
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
        Mesh mesh = scene.meshes[i];
        for (const auto& face : mesh.faces) {
            Triangle tri;
            tri.material_id = mesh.material_id;
            tri.indices = face;

            float t = intersectTriangle(shadowRayStruct, tri);
            if (t > 0 && t < distanceToLight) {
                return true;
            }
        }
    }

    return false;
}

// Phong shading function
vec3f calculateColor(int materialId, vec3f intersectionPoint, vec3f normal, ray myRay, int depth) {
    vec3f color = {0, 0, 0};
    Material material = scene.materials[materialId - 1];

    // Ambient component
    vec3f ambientLight = scene.ambient_light;
    color.x += material.ambient.x * ambientLight.x / 255;
    color.y += material.ambient.y * ambientLight.y / 255;
    color.z += material.ambient.z * ambientLight.z / 255;

    // Loop over all point lights
    // printf("start\n");
    for (const auto &light : scene.point_lights) {

        if (isInShadow(intersectionPoint, light.position)) {
            continue;
        }
        // printf("light_positionx: %f\n", light.position.x);
        // Diffuse component
        vec3f L = normalize(substractVectorsf(light.position, intersectionPoint)); // Light direction
        float distance = length(substractVectorsf(light.position, intersectionPoint));
        float cosTheta = fmax(0.0f, dot(L, normal));

        vec3f received_irradiance;
        received_irradiance.x = light.intensity.x / (distance * distance);
        received_irradiance.y = light.intensity.y / (distance * distance);
        received_irradiance.z = light.intensity.z / (distance * distance);

        color.x += material.diffuse.x * cosTheta * received_irradiance.x / 255;
        color.y += material.diffuse.y * cosTheta * received_irradiance.y / 255;
        color.z += material.diffuse.z * cosTheta * received_irradiance.z / 255;

        // Specular component
        vec3f viewDir = normalize(substractVectorsf(camera.position, intersectionPoint)); // scene.cameras[0].position
        vec3f halfDir = normalize(addVectorsf(L, viewDir));
        float cosAlpha = fmax(0.0f, dot(normal, halfDir));
        float specular = powf(cosAlpha, material.phong_exponent);

        color.x += material.specular.x * specular * received_irradiance.x / 255;
        color.y += material.specular.y * specular * received_irradiance.y / 255;
        color.z += material.specular.z * specular * received_irradiance.z / 255;
    }

    // Reflection component
    if (material.is_mirror && depth > 0) {
        vec3f reflectionDir = substractVectorsf(myRay.direction, multiplicationScalarf(normal, 2 * dot(myRay.direction, normal)));
        reflectionDir = normalize(reflectionDir);

        ray reflectionRay;
        reflectionRay.origin = addVectorsf(intersectionPoint, multiplicationScalarf(reflectionDir, scene.shadow_ray_epsilon));
        reflectionRay.direction = reflectionDir;

        vec3f reflectedColor = computeColor(reflectionRay, depth - 1);
        color.x += material.mirror.x * reflectedColor.x;
        color.y += material.mirror.y * reflectedColor.y;
        color.z += material.mirror.z * reflectedColor.z;
    }

    // Clamp color values to [0, 1]
    color.x = std::min(1.0f, std::max(0.0f, color.x));
    color.y = std::min(1.0f, std::max(0.0f, color.y));
    color.z = std::min(1.0f, std::max(0.0f, color.z));
    
    return color;
}

vec3f computeColor(ray myRay, int depth) {
    if (depth <= 0) {
        return {scene.background_color.x / 255.0f, scene.background_color.y / 255.0f, scene.background_color.z / 255.0f};
    }

    float minT = std::numeric_limits<float>::max(); // float minT = 9999999.9

    vec3f intersectionPoint, normal;
    int materialId = -1;
    bool hasIntersection = false;

    // sphere intersection
    for (int i = 0; i < scene.spheres.size(); i++) {
        float t = intersectSphere(myRay, scene.spheres[i]);
        if (t < minT && t >= 0.001) {
            minT = t;
            intersectionPoint = addVectorsf(myRay.origin, multiplicationScalarf(myRay.direction, t));
            normal = normalize(substractVectorsf(intersectionPoint, scene.vertex_data[scene.spheres[i].center_vertex_id - 1]));
            materialId = scene.spheres[i].material_id;
            hasIntersection = true;
        }
    }

    // triangle intersection
    for (int i = 0; i < scene.triangles.size(); i++) {
        float t = intersectTriangle(myRay, scene.triangles[i]);
        if (t < minT && t >= 0.001) {
            minT = t;
            intersectionPoint = addVectorsf(myRay.origin, multiplicationScalarf(myRay.direction, t));

            vec3f v0 = scene.vertex_data[scene.triangles[i].indices.v0_id - 1];
            vec3f v1 = scene.vertex_data[scene.triangles[i].indices.v1_id - 1];
            vec3f v2 = scene.vertex_data[scene.triangles[i].indices.v2_id - 1];
            normal = normalize(cross(substractVectorsf(v1, v0), substractVectorsf(v2, v0)));
            materialId = scene.triangles[i].material_id;
            hasIntersection = true;
        }
    }

    // mesh intersection
    for (int i = 0; i < scene.meshes.size(); i++) {
        Mesh mesh = scene.meshes[i];
        
        for (const auto& face : mesh.faces) {
            vec3f v0 = scene.vertex_data[face.v0_id - 1];
            vec3f v1 = scene.vertex_data[face.v1_id - 1];
            vec3f v2 = scene.vertex_data[face.v2_id - 1];
            Triangle tri;
            tri.material_id = mesh.material_id;
            tri.indices = face;
            float t = intersectTriangle(myRay, tri);

            if (t >= 0.001 && t < minT) {
                minT = t;
                intersectionPoint = addVectorsf(myRay.origin, multiplicationScalarf(myRay.direction, t));
                
                // Calculate the normal at the intersection point
                normal = normalize(cross(substractVectorsf(v1, v0), substractVectorsf(v2, v0)));
                materialId = mesh.material_id;
                hasIntersection = true;
            }
        }
    }

    // calculate the color
    if (hasIntersection) {
        return calculateColor(materialId, intersectionPoint, normal, myRay, depth);
    }

    return {scene.background_color.x / 255.0f, scene.background_color.y / 255.0f, scene.background_color.z / 255.0f};
    
}

int main(int argc, char* argv[])
{
    auto start = std::chrono::high_resolution_clock::now();

    scene.loadFromXml(argv[1]);

    // printf("%zu\n", scene.cameras.size());

    for (int k = 0; k <= scene.cameras.size()-1; k++) {

        // time measurement for the current camera
        auto camera_start = std::chrono::high_resolution_clock::now();

        // printf("%d\n", k);
        camera = scene.cameras[k];
        initializeCameraVectors(camera);

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
        
        // printf("%s\n", camera.image_name.c_str());
        write_ppm(camera.image_name.c_str(), image, width, height);
        // write_ppm("test.ppm", image, width, height);
        delete[] image;

        // End time measurement for the current camera
        auto camera_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> camera_elapsed_seconds = camera_end - camera_start;
        printf("Elapsed time for camera -> %s: %.4f s\n", camera.image_name.c_str(), camera_elapsed_seconds.count());

    }

    // end time measurement
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> total_elapsed_seconds = end - start;
    printf("Total elapsed time: %.4f s\n", total_elapsed_seconds.count());

    return 0;

}
