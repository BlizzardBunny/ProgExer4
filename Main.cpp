#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

using namespace glm;

struct Ray
{
    glm::vec3 origin; // Ray origin
    glm::vec3 direction; // Ray direction
};

struct Material
{
    glm::vec3 ambient; // Ambient
    glm::vec3 diffuse; // Diffuse
    glm::vec3 specular; // Specular
    float shininess; // Shininess
};

struct SceneObject
{
    Material material; // Material

    /**
     * Template function for calculating the intersection of this object with the provided ray.
     * @param[in]   incomingRay             Ray that will be checked for intersection with this object
     * @param[out]  outIntersectionPoint    Point of intersection (in case there is an intersection)
     * @param[out]  outIntersectionNormal   Normal vector at the point of intersection (in case there is an intersection)
     * @return If there is an intersection, returns the distance from the ray origin to the intersection point. Otherwise, returns a negative number.
     */
    virtual float Intersect(const Ray& incomingRay, glm::vec3& outIntersectionPoint, glm::vec3& outIntersectionNormal) = 0;
};

// Subclass of SceneObject representing a Sphere scene object
struct Sphere : public SceneObject
{
    glm::vec3 center; // center
    float radius; // radius

    /**
     * @brief Ray-sphere intersection
     * @param[in]   incomingRay             Ray that will be checked for intersection with this object
     * @param[out]  outIntersectionPoint    Point of intersection (in case there is an intersection)
     * @param[out]  outIntersectionNormal   Normal vector at the point of intersection (in case there is an intersection)
     * @return If there is an intersection, returns the distance from the ray origin to the intersection point. Otherwise, returns a negative number.
     */
    virtual float Intersect(const Ray& incomingRay, glm::vec3& outIntersectionPoint, glm::vec3& outIntersectionNormal)
    {
        //#5
        
        // In case there is an intersection, place the intersection point and intersection normal
        // that you calculated to the outIntersectionPoint and outIntersectionNormal variables.
        //
        // When you use this function from the outside, you can pass in the variables by reference.
        //
        // Example:
        // Ray ray = ...;
        // glm::vec3 point, normal;
        // float t = sphere->Intersect(ray, point, normal);
        //
        // (At this point, point and normal will now contain the intersection point and intersection normal)

        vec3 m = incomingRay.origin - center;
        float b = dot(m, incomingRay.direction);
        float c = dot(m, m) - (radius * radius);

        if ((b * b) - c < 0)
        {
            return -1.0f;
        }
        else if ((b * b) - c == 0)
        {
            float t = -b;
            outIntersectionPoint = incomingRay.origin + (t * incomingRay.direction);
            outIntersectionNormal = normalize(outIntersectionPoint - center);
            return t;
        }
        else
        {
            float t1 = -b + sqrt((b * b) - c);
            float t2 = -b - sqrt((b * b) - c);

            if (t1 < 0 && t2 < 0) //both negative
            {
                return -1.0f;
            }
            else
            {
                float t;
                if (t1 >= 0 && t2 < 0) //t1 is positive
                {
                    t = t1;
                }
                else if (t1 < 0 && t2 >= 0) //t2 is positive
                {
                    t = t2;
                }
                else //both are positive
                {
                    t = min(t1, t2);
                }

                outIntersectionPoint = incomingRay.origin + (t * incomingRay.direction);
                outIntersectionNormal = normalize(outIntersectionPoint - center);

                return t;
            }
        }
    }
};

// Subclass of SceneObject representing a Triangle scene object
struct Triangle : public SceneObject
{
    glm::vec3 A; // First point
    glm::vec3 B; // Second point
    glm::vec3 C; // Third point

    /**
     * @brief Ray-Triangle intersection
     * @param[in]   incomingRay             Ray that will be checked for intersection with this object
     * @param[out]  outIntersectionPoint    Point of intersection (in case there is an intersection)
     * @param[out]  outIntersectionNormal   Normal vector at the point of intersection (in case there is an intersection)
     * @return If there is an intersection, returns the distance from the ray origin to the intersection point. Otherwise, returns a negative number.
     */
    virtual float Intersect(const Ray& incomingRay, glm::vec3& outIntersectionPoint, glm::vec3& outIntersectionNormal)
    {
        //#9

        // The same idea for the outIntersectionPoint and outIntersectionNormal applies here
        vec3 n = cross(B - A, C - A);
        float f = dot(-incomingRay.direction, n);
        vec3 e = cross(-incomingRay.direction, incomingRay.origin - A);
        
        if (f <= 0)
        {
            return -1.0f;
        }
        
        float t = dot((incomingRay.origin - A), n) / f;
        float u = dot((C - A), e) / f;
        float v = dot(-(B - A), e) / f;

        if (t > 0)
        {
            if ((u + v <= 1) && (u >= 0) && (v >= 0))
            {
                outIntersectionPoint = incomingRay.origin + (t * incomingRay.direction);
                outIntersectionNormal = normalize(cross(B - A, C - A));
                return t;
            }
            else
            {
                return -1.0f;
            }
        }
        else
        {
            return -1.0f;
        }
        
    }
};

struct Camera
{
    glm::vec3 position; // Position
    glm::vec3 lookTarget; // Look target
    glm::vec3 globalUp; // Global up-vector
    float fovY; // Vertical field of view
    float focalLength; // Focal length

    int imageWidth; // image width
    int imageHeight; // image height
};

struct Light
{
    glm::vec4 position; // Light position (w = 1 if point light, w = 0 if directional light)

    glm::vec3 ambient; // Light's ambient intensity
    glm::vec3 diffuse; // Light's diffuse intensity
    glm::vec3 specular; // Light's specular intensity

    // --- Attenuation variables ---
    float constant; // Constant factor
    float linear; // Linear factor
    float quadratic; // Quadratic factor
};

struct IntersectionInfo
{
    Ray incomingRay; // Ray used to calculate the intersection
    float t; // Distance from the ray's origin to the point of intersection (if there was an intersection).
    SceneObject* obj; // Object that the ray intersected with. If this is equal to nullptr, then no intersection occured.
    glm::vec3 intersectionPoint; // Point where the intersection occured (if there was an intersection)
    glm::vec3 intersectionNormal; // Normal vector at the point of intersection (if there was an intersection)
};

struct Scene
{
    std::vector<SceneObject*> objects; // List of all objects in the scene
    std::vector<Light*> lights; // List of all lights in the scene
};

struct Image
{
    std::vector<unsigned char> data; // Image data
    int width; // Image width
    int height; // Image height

    /**
     * @brief Constructor
     * @param[in] w Width
     * @param[in] h Height
     */
    Image(const int& w, const int& h)
        : width(w)
        , height(h)
    {
        data.resize(w * h * 3, 0);
    }

    /**
     * @brief Converts the provided color value from [0, 1] to [0, 255]
     * @param[in] c Color value in [0, 1] range
     * @return Color value in [0, 255] range
     */
    unsigned char ToChar(float c)
    {
        c = glm::clamp(c, 0.0f, 1.0f);
        return static_cast<unsigned char>(c * 255);
    }

    /**
     * @brief Sets the color at the specified pixel location
     * @param[in] x     X-coordinate of the pixel
     * @param[in] y     Y-coordinate of the pixel
     * @param[in] color Pixel color
     */
    void SetColor(const int& x, const int& y, const glm::vec3& color)
    {
        int index = (y * width + x) * 3;
        data[index] = ToChar(color.r);
        data[index + 1] = ToChar(color.g);
        data[index + 2] = ToChar(color.b);
    }
};

/**
 * @brief Gets the ray that goes from the camera's position to the specified pixel at (x, y)
 * @param[in] camera Camera data
 * @param[in] x X-coordinate of the pixel (upper-left corner of the pixel)
 * @param[in] y Y-coordinate of the pixel (upper-left corner of the pixel)
 * @return Ray that passes through the pixel at (x, y)
 */
Ray GetRayThruPixel(const Camera &camera, const int& pixelX, const int& pixelY)
{
    //#4
    Ray ray;
    ray.origin = camera.position; 

    float aspect = (float) camera.imageWidth / camera.imageHeight;
    float hViewport = 2.0f * camera.focalLength * tan(radians(camera.fovY / 2));
    float wViewport = aspect * hViewport;

    vec3 lookDirection = normalize(camera.lookTarget - camera.position);

    vec3 u = cross(lookDirection, camera.globalUp) / length(cross(lookDirection, camera.globalUp));
    vec3 v = cross(u, lookDirection) / length(cross(u, lookDirection));

    vec3 L = camera.position + (lookDirection * camera.focalLength) - (u * (wViewport / 2)) - (v * (hViewport / 2));

    float s = ((pixelX + 0.5f) * wViewport) / camera.imageWidth;
    float t = ((pixelY + 0.5f) * hViewport) / camera.imageHeight;

    vec3 P = L + (u * s) + (v * t);

    ray.direction = (P - ray.origin)/length(P - ray.origin);

    return ray;
}

/**
 * @brief Cast a ray to the scene.
 * @param[in] ray   Ray to cast to the scene
 * @param[in] scene Scene object
 * @return Returns an IntersectionInfo object that will contain the results of the raycast
 */
IntersectionInfo Raycast(const Ray& ray, const Scene &scene)
{
    //#6
    IntersectionInfo ret;
    ret.incomingRay = ray;
    ret.obj = nullptr;

    //// Fields that need to be populated:
    //ret.intersectionPoint = glm::vec3(0.0f); // Intersection point
    //ret.intersectionNormal = glm::vec3(0.0f); // Intersection normal
    //ret.t = 0.0f; // Distance from ray origin to intersection point
    //ret.obj = nullptr; // First object hit by the ray. Set to nullptr if the ray does not hit anything
    float dist = -1.0f;
    for (int i = 0; i < scene.objects.size(); i++)
    {
        vec3 point, normal;
        float t = scene.objects[i]->Intersect(ray, point, normal);
        
        float tempDist = length(point - ray.origin);

        if ((tempDist < dist) || (dist < 0))
        {            
            if (t >= 0)
            {
                dist = tempDist;
                ret.intersectionPoint = point;
                ret.intersectionNormal = normal;
                ret.t = t;
                ret.obj = scene.objects[i];
            }
        }
    }

    return ret;
}

/**
 * @brief Perform a ray-trace to the scene
 * @param[in] ray       Ray to trace
 * @param[in] scene     Scene data
 * @param[in] camera    Camera data
 * @param[in] maxDepth  Maximum depth of the trace
 * @return Resulting color after the ray bounced around the scene
 */
glm::vec3 RayTrace(const Ray& ray, const Scene& scene, const Camera& camera, int maxDepth = 1)
{
    //#7
    glm::vec3 finalColor(0.0f);

    IntersectionInfo object = Raycast(ray, scene);

    
    if (object.obj != nullptr)
    {
        //#13 & #14
        glm::vec3 otherCompSum(0.0f); //sum of diffuse and specular components
        glm::vec3 ambientSum(0.0f); //sum of ambient components
        for (int i = 0; i < scene.lights.size(); i++)
        {
            ambientSum += scene.lights[i]->ambient;

            vec3 n = object.intersectionNormal;

            vec3 position = vec3(scene.lights[i]->position.x, scene.lights[i]->position.y, scene.lights[i]->position.z);
            vec3 l;

            if (scene.lights[i]->position.w > 0)
            {
                l = normalize(position - object.intersectionPoint);
            }
            else
            {
                l = scene.lights[i]->position;
            }
                        
            vec3 diffuseComponent = max(dot(n, l), 0.0f) * object.obj->material.diffuse * scene.lights[i]->diffuse;

            vec3 r = reflect(-l, n);
            vec3 v = camera.position - object.intersectionPoint;

            vec3 specularComponent = max(dot(r, v), 0.0f) * object.obj->material.specular * scene.lights[i]->specular;

            float d = length(position - object.intersectionPoint);
            float attenuation = (scene.lights[i]->quadratic * d * d) + (scene.lights[i]->linear * d) + scene.lights[i]->constant;

            otherCompSum += (diffuseComponent + specularComponent) * attenuation;

            //#15
            vec3 s = scene.lights[i]->position;
            float ObjectDistance = sqrt((pow(object.intersectionPoint.x - camera.position.x, 2)) + (pow(object.intersectionPoint.y - camera.position.y, 2)) + (pow(object.intersectionPoint.z - camera.position.z, 2)));
            float LightDistance = sqrt((pow(s.x - camera.position.x, 2)) + (pow(s.y - camera.position.y, 2)) + (pow(s.z - camera.position.z, 2)));

            if ((object.obj = nullptr) && (ObjectDistance < LightDistance))
            {
                glm::vec3 ambientSum(0.0f);
                glm::vec3 bias(0.001f); //value of bias used for shadow acne
                for (int i = 0; i < scene.lights.size(); i++)
                {
                    vec3 n = object.intersectionNormal;
                    object.intersectionPoint += n * bias;
                    ambientSum += scene.lights[i]->ambient;
                    finalColor = ambientSum;
                }
            }

            //#16
            if (maxDepth > 1)
            {
                while (maxDepth > 1)
                {
                    float k = material.shininess/128; //Reflectivity
                    Ray NewRay;
                    IntersectionInfo object = Raycast(NewRay, scene);
                    glm::vec3 NewColor = RayTrace(NewRay, scene, camera, maxDepth);
                    NewColor = NewColor * k;
                    finalColor += NewColor;
                    maxDepth--;
                }
            }
        }

        finalColor = (ambientSum / (float)(scene.lights.size())) + otherCompSum;
    }
    return finalColor;
}

//#8 run app

/**
 * Main function
 */
int main()
{
    Scene scene;
    Camera camera; 
    int maxDepth;

    /*
    //#2
    camera.imageWidth = 640;
    camera.imageHeight = 480;
    camera.position = vec3(0.0f, 0.0f, 3.0f);
    camera.lookTarget = vec3(0.0f, 0.0f, 0.0f);
    camera.globalUp = vec3(0.0f, 1.0f, 0.0f);
    camera.fovY = 45.0f;
    camera.focalLength = 1.0f;

    //#3
    Sphere* sphere = new Sphere();
    sphere->center = vec3(0.0f, 0.0f, 0.0f);
    sphere->radius = 1.0f;
    sphere->material.diffuse = vec3(1.0f, 0.0f, 0.0f);
    scene.objects.push_back(sphere);

    //#10
    Triangle* triangle = new Triangle();
    triangle->A = vec3(2.0f, -2.0f, -2.0f);
    triangle->B = vec3(0.0f, 2.0f, -2.0f);
    triangle->C = vec3(-2.0f, -2.0f, -2.0f);
    triangle->material.diffuse = vec3(0.0f, 1.0f, 0.0f);
    scene.objects.push_back(triangle);*/


    //#11
    std::cout << "Enter file name: " << std::endl;
    std::string filename;
    std::cin >> filename;
    std::ifstream sceneTest;
    sceneTest.open(filename.c_str()); 
    if (!sceneTest) { // file couldn't be opened
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }    
    
    std::vector<std::string> values;
    std::string temp;
    while (sceneTest >> temp)
    {
        values.push_back(temp);
    }
    
    try
    {
        int i = -1;
        camera.imageWidth = std::stoi(values[++i]);
        camera.imageHeight = std::stoi(values[++i]);

        camera.position.x = std::stof(values[++i]);
        camera.position.y = std::stof(values[++i]);
        camera.position.z = std::stof(values[++i]);

        camera.lookTarget.x = std::stof(values[++i]);
        camera.lookTarget.y = std::stof(values[++i]);
        camera.lookTarget.z = std::stof(values[++i]);

        camera.globalUp.x = std::stof(values[++i]);
        camera.globalUp.y = std::stof(values[++i]);
        camera.globalUp.z = std::stof(values[++i]);

        camera.fovY = std::stof(values[++i]);
        camera.focalLength = std::stof(values[++i]);

        maxDepth = std::stoi(values[++i]);
        int N = std::stoi(values[++i]);
        for (int j = 0; j < N; j++)
        {
            std::string objectName = values[++i];

            SceneObject* object = nullptr;
            if (objectName == "sphere")
            {
                Sphere* sphere = new Sphere();

                sphere->center.x = std::stof(values[++i]);
                sphere->center.y = std::stof(values[++i]);
                sphere->center.z = std::stof(values[++i]);

                sphere->radius = std::stof(values[++i]);

                object = sphere;
            }
            else if (objectName == "tri")
            {
                Triangle* triangle = new Triangle();

                triangle->A.x = std::stof(values[++i]);
                triangle->A.y = std::stof(values[++i]);
                triangle->A.z = std::stof(values[++i]);

                triangle->B.x = std::stof(values[++i]);
                triangle->B.y = std::stof(values[++i]);
                triangle->B.z = std::stof(values[++i]);

                triangle->C.x = std::stof(values[++i]);
                triangle->C.y = std::stof(values[++i]);
                triangle->C.z = std::stof(values[++i]);

                object = triangle;
            }

            if (object != nullptr)
            {
                object->material.ambient.r = std::stof(values[++i]);
                object->material.ambient.g = std::stof(values[++i]);
                object->material.ambient.b = std::stof(values[++i]);

                object->material.diffuse.r = std::stof(values[++i]);
                object->material.diffuse.g = std::stof(values[++i]);
                object->material.diffuse.b = std::stof(values[++i]);

                object->material.specular.r = std::stof(values[++i]);
                object->material.specular.g = std::stof(values[++i]);
                object->material.specular.b = std::stof(values[++i]);

                object->material.shininess = std::stof(values[++i]);
            }

            scene.objects.push_back(object);
        }

        int L = std::stoi(values[++i]);
        for (int j = 0; j < L; j++)
        {
            Light* light = new Light();

            light->position.x = std::stof(values[++i]);
            light->position.y = std::stof(values[++i]);
            light->position.z = std::stof(values[++i]);

            light->ambient.x = std::stof(values[++i]);
            light->ambient.y = std::stof(values[++i]);
            light->ambient.z = std::stof(values[++i]);

            light->diffuse.x = std::stof(values[++i]);
            light->diffuse.y = std::stof(values[++i]);
            light->diffuse.z = std::stof(values[++i]);

            light->specular.x = std::stof(values[++i]);
            light->specular.y = std::stof(values[++i]);
            light->specular.z = std::stof(values[++i]);

            light->constant = std::stof(values[++i]);
            light->linear = std::stof(values[++i]);
            light->quadratic = std::stof(values[++i]);

            scene.lights.push_back(light);
        }        

    }
    catch (...)
    {
        std::cerr << "File not formatted properly" << std::endl;
        exit(1);
    }

    Image image(camera.imageWidth, camera.imageHeight);
    for (int y = 0; y < image.height; ++y)
    {
        for (int x = 0; x < image.width; ++x)
        {
            Ray ray = GetRayThruPixel(camera, x, image.height - y - 1);

            glm::vec3 color = RayTrace(ray, scene, camera, maxDepth);
            image.SetColor(x, y, color);
        }

        std::cout << "Row: " << std::setfill(' ') << std::setw(4) << (y + 1) << " / " << std::setfill(' ') << std::setw(4) << image.height << "\r" << std::flush;
    }
    std::cout << std::endl;
    
    std::string imageFileName = "scene.png"; // You might need to make this a full path if you are on Mac
    stbi_write_png(imageFileName.c_str(), image.width, image.height, 3, image.data.data(), 0);
    
    for (size_t i = 0; i < scene.objects.size(); ++i)
    {
        delete scene.objects[i];
    }

    return 0;
}

