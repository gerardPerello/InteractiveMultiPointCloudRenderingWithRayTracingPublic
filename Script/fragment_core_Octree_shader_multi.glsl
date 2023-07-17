#version 440

#define RAND_MAX 0x7fff

const int MAX_STACK_SIZE = 64;
const int MAX_HITS = 16;

struct Material
{
	vec3 ambient;
	vec3 diffuse;
	vec3 specular;
	float shininess;
	float opacity;
  vec3 kt;
	sampler2D diffuseTex;
	sampler2D specularTex;
	
};

struct PointLight
{
	vec3 position;
	vec3 Ia;
	vec3 Id;
	vec3 Is;
	float intensity;
	float constant;
	float linear;
	float quadratic;
	
};

struct Ray
{
	vec3 origin;
	vec3 direction;	
};

struct Point 
{
  vec3 position;
  vec3 color;
  vec3 normal;
  vec3 label;
};

struct HitData
{
	int sphereID;
	float t;
	vec3 hitPosition;
	vec3 hitNormal;
	bool hitted;
    bool plane;
	vec4 hitColor;
    bool error;
    Point hitPoint;
};


struct FlattenedMultiNodeOctreeInfo {
    int pointIndex;
    int pointCount;
    int depends;
    int leaf;
};

struct FlattenedMultiNode {
    vec3 minBound;
    vec3 maxBound;
    ivec4 childrenInfo1;
    ivec4 childrenInfo2;
    FlattenedMultiNodeOctreeInfo octreeInfo[3];
};


//Se usa
in vec3 vs_position;

out vec4 fs_color;

//Uniforms
uniform Material material;
//uniform PointLight pointlight;
#define MAX_POINT_LIGHTS 3

uniform PointLight pointLights[MAX_POINT_LIGHTS];
uniform int nLights;
//UniformsUnics
uniform vec3 camPosition;
uniform mat4 ViewMatrix;
uniform mat4 ProjectionMatrix;
uniform float uEpsilon;
uniform float uMaxRayDistance;
uniform vec2 window_size;
uniform vec3 globalLight;
uniform vec2 randomSeed;
uniform vec3 planeOrigin;
uniform float planeSize;
uniform vec3 range;

//Uniforms Blending
uniform float radius;
uniform float radiusError;
uniform float radiusExtra;
uniform float alphaGround;
uniform float alphaError;
uniform float alphaExtra;

//UnfiromsMode
uniform int mode;
uniform int octreeError;
uniform int octreeExtra;
uniform int octreeGround;
uniform int shadowMode;
uniform int orderMode;
uniform int shadingMode;
uniform int planesMode;

//buffers
uniform samplerBuffer u_multiNodesTextureFloat;
uniform isamplerBuffer u_multiNodesTextureInt; 
uniform samplerBuffer u_multiPointsTexture1;
uniform samplerBuffer u_multiPointsTexture2;
uniform samplerBuffer u_multiPointsTexture3;






//Metodo que se encarga de mezclar dos colores.
/*
En esta versión, alfa1 tiene prioridad en el cálculo del nuevo alfa:

float newAlpha = min(1.0f, alfa1 + alfa2 * (1.0f - alfa1));

También modificamos la línea donde se calcula el nuevo color:

vec3 newColor = (color1 * alfa1 + color2 * alfa2 * (1.0f - alfa1)) / newAlpha;

Ahora, el peso de color1 es alfa1, y el peso de color2 es alfa2 * (1.0 - alfa1). Esto significa que si alfa1 es más grande, color1 tendrá más influencia en el resultado final.

Con estas modificaciones, la función ahora dará prioridad al color1 y al alfa1, como deseas.
*/
vec4 mixColor(vec3 color1, float alfa1, vec3 color2, float alfa2)
{
    float newAlpha = min(1.0f, alfa1 + alfa2 * (1.0f - alfa1));
    vec3 newColor = (color1 * alfa1 + color2 * alfa2 * (1.0f - alfa1)) / newAlpha;
    return vec4(newColor, newAlpha);
}


// Función drand48: devuelve un valor aleatorio entre -1 y 1 usando una técnica de generación de ruido
// Agrega un parámetro adicional 'vec2 seed' para mejorar la aleatoriedad
float drand48(vec2 co) {
  // Genera un valor aleatorio utilizando las coordenadas de entrada, el seed y constantes específicas
  return 2 * fract(sin(dot(co.xy + randomSeed, vec2(12.9898, 78.233))) * 43758.5453) - 1;
}

// Función random_in_unit_disk: genera un punto aleatorio dentro de un disco unitario
// Agrega un parámetro adicional 'vec2 seed' para mejorar la aleatoriedad
vec3 random_in_unit_disk(vec2 co) {
  vec3 p;
  int n = 0;
  do {
    // Genera un punto aleatorio utilizando la función drand48 y el seed
    p = vec3(drand48(co.xy), drand48(co.yx), 0);
    n++;
  // Verifica si el punto está dentro del disco unitario y si el número de intentos es menor a 3
  } while (dot(p, p) >= 1.0 && n < 3);
  return p;
}

// Función pointAtParameter: devuelve un punto en el rayo a una distancia t desde el origen
vec3 pointAtParameter(Ray ray, float t) {
  return ray.origin + t * ray.direction;
}

//TODO CAMBIAR AIXO
// Función get_ray: crea un rayo desde la posición de la cámara que pasa a través de un punto en la pantalla
Ray get_ray(vec2 screenPos) {
  // Convierte la posición de la pantalla a coordenadas de clip
  vec4 ray_clip = vec4(screenPos, -1.0, 1.0);
  
  // Convierte las coordenadas de clip a coordenadas de espacio ocular
  vec4 ray_eye = ProjectionMatrix * ray_clip;
  ray_eye = vec4(ray_eye.xy, -1.0, 0.0);
  
  // Convierte las coordenadas del espacio ocular a coordenadas del espacio mundial
  vec3 ray_wor = (ViewMatrix * ray_eye).xyz;
  
  // Normaliza la dirección del rayo y agrega un pequeño valor para evitar problemas numéricos
  ray_wor = normalize(ray_wor + vec3(0.0001, 0.0001, 0.0001));
  
  // Devuelve el rayo con la posición de la cámara como origen y la dirección calculada
  return Ray(camPosition, ray_wor);
}


// La función plane_hit verifica si un rayo interseca alguno de los tres planos
// definidos por su punto de intersección y tamaño.
bool plane_hit(Ray r, float t_min, float t_max, vec3 intersectionPoint, float planeSize, vec3 range, out HitData hit) {
  // Array para almacenar los datos de las intersecciones con los tres planos
  HitData hits[3];
  
  // Normales de los tres planos
  vec3 normals[3] = vec3[](vec3(1, 0, 0), vec3(0, 1, 0), vec3(0, 0, 1));
  
  // Variable para registrar si el rayo intersecta al menos un plano
  bool any_hit = false;
  float t;

  // Itera a través de los tres planos
  for (int i = 0; i < 3; i++) {
    // Calcula el valor t de la intersección con el plano actual
    t = (intersectionPoint[i] - r.origin[i]) / r.direction[i];

    // Verifica si t está dentro del rango válido
    if (t >= t_min && t <= t_max) {
      // Calcula la posición de intersección
      vec3 hitPosition = pointAtParameter(r, t);
      // Verifica si la posición de intersección está dentro de los límites del plano y dentro del rango especificado
      
    if (hitPosition[(i + 1) % 3] > intersectionPoint[(i + 1) % 3] &&
        hitPosition[(i + 2) % 3] > intersectionPoint[(i + 2) % 3] &&
        hitPosition[(i + 1) % 3] < intersectionPoint[(i + 1) % 3] + range[(i + 1) % 3] &&
        hitPosition[(i + 2) % 3] < intersectionPoint[(i + 2) % 3] + range[(i + 2) % 3] &&
        abs(hitPosition[i] - intersectionPoint[i]) <= planeSize) {
        
        // Almacena los datos de la intersección en el array
        hits[i].t = t;
        hits[i].hitPosition = hitPosition;
        hits[i].hitNormal = normals[i];
        hits[i].hitColor = vec4(1.0,1.0,1.0,1.0);//mixColor(hit.hitColor.xyz, hit.hitColor.w, vec3(1,1,1), 1.f);
        hits[i].hitted = true;
        hits[i].plane = true;
        
        // Indica que se encontró al menos una intersección
        any_hit = true;
      }
    }
  }

  // Si se encontró al menos una intersección
  if (any_hit) {
    float min_t = t_max;
    int min_index = -1;

    // Encuentra la intersección más cercana
    for (int i = 0; i < 3; i++) {
      if (hits[i].hitted && hits[i].t < min_t) {
        min_t = hits[i].t;
        min_index = i;
      }
    }

    // Si se encuentra una intersección más cercana, actualiza los datos de salida
    if (min_index >= 0) {
      hit = hits[min_index];
      return true;
    }
  }

  // Si no se encontró ninguna intersección, retorna false
  return false;
}

FlattenedMultiNode getMultiNode(int index) {
    vec3 minBound = texelFetch(u_multiNodesTextureFloat, index * 2).xyz;
    vec3 maxBound = texelFetch(u_multiNodesTextureFloat, index * 2 + 1).xyz;

    ivec4 pointInfo1 = texelFetch(u_multiNodesTextureInt, index * 5);
    ivec4 pointInfo2 = texelFetch(u_multiNodesTextureInt, index * 5 + 1);
    ivec4 pointInfo3 = texelFetch(u_multiNodesTextureInt, index * 5 + 2);

    ivec4 childrenInfo1 = texelFetch(u_multiNodesTextureInt, index * 5 + 3);
    ivec4 childrenInfo2 = texelFetch(u_multiNodesTextureInt, index * 5 + 4);

    FlattenedMultiNode multiNode;
    multiNode.minBound = minBound;
    multiNode.maxBound = maxBound;
    multiNode.childrenInfo1 = childrenInfo1;
    multiNode.childrenInfo2 = childrenInfo2;

    for(int i = 0; i < 3; ++i) {
        multiNode.octreeInfo[i].pointIndex = i == 0 ? pointInfo1.x : (i == 1 ? pointInfo2.x : pointInfo3.x);
        multiNode.octreeInfo[i].pointCount = i == 0 ? pointInfo1.y : (i == 1 ? pointInfo2.y : pointInfo3.y);
        multiNode.octreeInfo[i].depends = i == 0 ? pointInfo1.z : (i == 1 ? pointInfo2.z : pointInfo3.z);
        multiNode.octreeInfo[i].leaf = i == 0 ? pointInfo1.w : (i == 1 ? pointInfo2.w : pointInfo3.w);
    }

    return multiNode;
}

Point getMultiPoint1(int index){
  vec3 position = texelFetch(u_multiPointsTexture1, index * 4).xyz;
  vec3 color = texelFetch(u_multiPointsTexture1, index * 4 + 1).xyz;
  vec3 normal = texelFetch(u_multiPointsTexture1, index * 4 + 2).xyz;
  vec3 label = texelFetch(u_multiPointsTexture1, index * 4 + 3).xyz;
  Point point;
  point.position = position;
  point.color = color;
  point.normal = normal;
  point.label = label;
  return point;
}
Point getMultiPoint2(int index){
  vec3 position = texelFetch(u_multiPointsTexture2, index * 4).xyz;
  vec3 color = texelFetch(u_multiPointsTexture2, index * 4 + 1).xyz;
  vec3 normal = texelFetch(u_multiPointsTexture2, index * 4 + 2).xyz;
  vec3 label = texelFetch(u_multiPointsTexture2, index * 4 + 3).xyz;
  Point point;
  point.position = position;
  point.color = color;
  point.normal = normal;
  point.label = label;
  return point;
}
Point getMultiPoint3(int index){
  vec3 position = texelFetch(u_multiPointsTexture3, index * 4).xyz;
  vec3 color = texelFetch(u_multiPointsTexture3, index * 4 + 1).xyz;
  vec3 normal = texelFetch(u_multiPointsTexture3, index * 4 + 2).xyz;
  vec3 label = texelFetch(u_multiPointsTexture3, index * 4 + 3).xyz;
  Point point;
  point.position = position;
  point.color = color;
  point.normal = normal;
  point.label = label;
  return point;
}


bool rayIntersectsMultiNode(Ray ray, FlattenedMultiNode node, float t_min, float t_max, out float t_entry, out float t_exit, float epsilon) {
    vec3 expandedMinBound = node.minBound - vec3(epsilon, epsilon, epsilon);
    vec3 expandedMaxBound = node.maxBound + vec3(epsilon, epsilon, epsilon);

    vec3 rayEntry = (expandedMinBound - ray.origin) / ray.direction;
    vec3 rayExit = (expandedMaxBound - ray.origin) / ray.direction;

    vec3 t_min_v = min(rayEntry, rayExit);
    vec3 t_max_v = max(rayEntry, rayExit);

    t_entry = max(max(t_min_v.x, t_min_v.y), max(t_min_v.z, t_min));
    t_exit = min(min(t_max_v.x, t_max_v.y), min(t_max_v.z, t_max));

    return t_entry <= t_exit;
}

bool rayIntersectsMultiNode2(Ray ray, FlattenedMultiNode node, float t_min, float t_max, float epsilon) {
    vec3 expandedMinBound = node.minBound - vec3(epsilon, epsilon, epsilon);
    vec3 expandedMaxBound = node.maxBound + vec3(epsilon, epsilon, epsilon);

    vec3 rayEntry = (expandedMinBound - ray.origin) / ray.direction;
    vec3 rayExit = (expandedMaxBound - ray.origin) / ray.direction;

    vec3 t_min_v = min(rayEntry, rayExit);
    vec3 t_max_v = max(rayEntry, rayExit);

    float t_entry = max(max(t_min_v.x, t_min_v.y), max(t_min_v.z, t_min));
    float t_exit = min(min(t_max_v.x, t_max_v.y), min(t_max_v.z, t_max));

    return t_entry <= t_exit;
}


bool intersectSphere(Ray ray, vec3 center, float radius, inout float t) {
    vec3 oc = ray.origin - center;
    float a = dot(ray.direction, ray.direction);
    float b = 2.0 * dot(oc, ray.direction);
    float c = dot(oc, oc) - radius * radius;
    float discriminant = b * b - 4.0 * a * c;

    if (discriminant < 0.0) {
        return false;
    } else {
        float t0 = (-b - sqrt(discriminant)) / (2.0 * a);
        float t1 = (-b + sqrt(discriminant)) / (2.0 * a);

        if (t0 < t1 && t0 > 0.0) {
            t = t0;
        } else if (t1 > 0.0) {
            t = t1;
        } else {
            return false;
        }
        return true;
    }
}

bool intersectCube(Ray ray, vec3 center, float sideLength, inout float t) {
    float halfSide = sideLength * 0.5;
    vec3 minBound = center - vec3(halfSide);
    vec3 maxBound = center + vec3(halfSide);

    vec3 rayEntry = (minBound - ray.origin) / ray.direction;
    vec3 rayExit = (maxBound - ray.origin) / ray.direction;

    vec3 tmin = min(rayEntry, rayExit);
    vec3 tmax = max(rayEntry, rayExit);

    float t_entry = max(max(tmin.x, tmin.y), tmin.z);
    float t_exit = min(min(tmax.x, tmax.y), tmax.z);

    if (t_entry <= t_exit && t_exit >= 0.0) {
        //t = t_entry;
        /* vec3 epsilon = vec3(0.0001);
        vec3 intersection = ray.origin + ray.direction * t;
        
        if (abs(intersection.x - minBound.x) < epsilon.x) {
            normal = vec3(-1, 0, 0);
        } else if (abs(intersection.x - maxBound.x) < epsilon.x) {
            normal = vec3(1, 0, 0);
        } else if (abs(intersection.y - minBound.y) < epsilon.y) {
            normal = vec3(0, -1, 0);
        } else if (abs(intersection.y - maxBound.y) < epsilon.y) {
            normal = vec3(0, 1, 0);
        } else if (abs(intersection.z - minBound.z) < epsilon.z) {
            normal = vec3(0, 0, -1);
        } else if (abs(intersection.z - maxBound.z) < epsilon.z) {
            normal = vec3(0, 0, 1);
        } */

        return true;
    } else {
        return false;
    }
}

void selectColor(Point point,out vec4 currentHitColor)
{
    if(mode == 1) {
        currentHitColor = (point.label.z == 1.f) ? vec4(point.color, alphaGround) : vec4(1.f,0.f,0.f, alphaError);
    } else if(mode == 2) {
        currentHitColor = (point.label.z == 1.f) ? vec4(point.color, alphaGround) : vec4(point.normal, alphaError);
    } else if(mode == 3) {
        currentHitColor = (point.label.z == 1.f) ? vec4(point.color, alphaGround) : vec4(point.color, alphaError);
    }
}

void setColorHitNew(out HitData hit, Point point, Ray ray, float t_entry) {

    //Escojemos el color dependiendo del modo en el que estamos. 
    //De momento, en color se guardan el original
    //En normal se guarda el predicted.
    vec4 currentHitColor;
    selectColor(point, currentHitColor);
    

    //En caso de que aun no hayamos chocado con nada, guardamos los datos del primer hit para saber donde hemos chocado
    //Luego, seguimos avanzando.
    if(hit.hitted == false)
    {
        hit.t = t_entry;
        hit.hitPosition = pointAtParameter(ray, t_entry);
        hit.hitNormal = normalize(hit.hitPosition - point.position);
        hit.hitted = true;
        hit.hitColor = currentHitColor;
        hit.hitPoint = point;
        
    }
    else
    {
        hit.hitColor = mixColor(hit.hitColor.xyz, hit.hitColor.w, currentHitColor.xyz, currentHitColor.w);
    }
    //En caso de que ya hayamos chocado con algo, mezclamos los colores.
    
}


void getVisitOrder(vec3 rayDirOctant, out int visitOrder[8])
{
    //Los ordenamos al reves porque luego los añadimos al stack al reves.
    if (rayDirOctant.x >= 0) {
        if (rayDirOctant.y >= 0) {
            if (rayDirOctant.z >= 0) {
                int order[] = int[](7, 6, 5, 4, 3, 2, 1, 0);
                visitOrder = order;
            } else {
                int order[] = int[](3, 2, 1, 0, 7, 6, 5, 4);
                visitOrder = order;
            }
        } else {
            if (rayDirOctant.z >= 0) {
                int order[] = int[](5, 4, 7, 6, 1, 0, 3, 2);
                visitOrder = order;
            } else {
                int order[] = int[](1, 0, 3, 2, 5, 4, 7, 6);
                visitOrder = order;
            }
        }
    } 
    else {
        if (rayDirOctant.y >= 0) {
            if (rayDirOctant.z >= 0) {
                int order[] = int[](6, 7, 4, 5, 2, 3, 0, 1);
                visitOrder = order;
            } else {
                int order[] = int[](2, 3, 0, 1, 6, 7, 4, 5);
                visitOrder = order;
            }
        } else {
            if (rayDirOctant.z >= 0) {
                int order[] = int[](4, 5, 6, 7, 0, 1, 2, 3);
                visitOrder = order;
            } else {
                int order[] = int[](0, 1, 2, 3, 4, 5, 6, 7);
                visitOrder = order;
            }
        }
    }
}

struct SphereHitData
{
    HitData hitData;
    float t;
    Point hitPoint;
};
SphereHitData hits[MAX_HITS];

void insertSorted(out SphereHitData hits[MAX_HITS],out int numHits, SphereHitData newHit) {
    int i = numHits - 1;
    while (i >= 0 && hits[i].t > newHit.t-uEpsilon) {
        hits[i + 1] = hits[i];
        i = i - 1;
    }
    hits[i + 1] = newHit;
    numHits = numHits + 1;
}


void visitChildren(FlattenedMultiNode currentNode, out int stack[MAX_STACK_SIZE], out int stackIndex, Ray ray, float t_min, float t_max)
{
    // Visitar los hijos en el orden adecuado
    vec3 rayDirOctant = sign(ray.direction);
    int visitOrder[8];

    //Da el orden al reves porque realmente los estamos metiendo en un stack asi que el primero que metemos es el ultimo que visitamos.
    getVisitOrder(rayDirOctant, visitOrder);

    for (int i = 0; i < 8; ++i) {
        int childIndex = visitOrder[i];
        int childNodeIndex;
        if (childIndex < 4) {
            childNodeIndex = currentNode.childrenInfo1[childIndex];
        }
        else {
            childNodeIndex = currentNode.childrenInfo2[childIndex - 4];
        }
        //Guardamos solo los nodos con los que chocamos en el orden correspondiente.
        FlattenedMultiNode childNode = getMultiNode(childNodeIndex);
        if (rayIntersectsMultiNode2(ray, childNode, t_min, t_max, radius*3)) {
            if (childNodeIndex > 0 && stackIndex < MAX_STACK_SIZE) {
                stack[stackIndex++] = childNodeIndex;
            }
        }
    }
}


bool evaluateMultiNodeWithNoOrder(FlattenedMultiNode currentNode,Ray ray, float t_min, float t_max,float t_entry, float t_exit, out HitData hit, bool dependsGround, bool dependsError, bool dependsExtra, bool shadowing)
{
    
    // Comprobar si el rayo intersecta con alguno de los puntos
    int numHits = 0;
    if(dependsGround){
        for (int i = currentNode.octreeInfo[0].pointIndex; i < currentNode.octreeInfo[0].pointIndex + currentNode.octreeInfo[0].pointCount; ++i) {
            Point point = getMultiPoint1(i);
            if (intersectSphere(ray, point.position, abs(radius), t_entry)) {
                if(shadowing)
                {
                    hit.t = t_entry;
                    hit.hitPosition = pointAtParameter(ray, t_entry);
                    hit.hitNormal = normalize(hit.hitPosition - point.position);
                    hit.hitted = true;
                    selectColor(point, hit.hitColor);
                    return true;
                }
                else{
                    setColorHitNew(hit, point, ray, t_entry);
                    if(hit.hitColor.w == 1.0f){ 
                    return true;}
                }
            }
        }

        
    }
    if(dependsError){
        for (int i = currentNode.octreeInfo[1].pointIndex; i < currentNode.octreeInfo[1].pointIndex + currentNode.octreeInfo[1].pointCount; ++i) {
            Point point = getMultiPoint2(i);
            if (intersectSphere(ray, point.position, abs(radiusError), t_entry)) {
                if(shadowing)
                {
                    hit.t = t_entry;
                    hit.hitPosition = pointAtParameter(ray, t_entry);
                    hit.hitNormal = normalize(hit.hitPosition - point.position);
                    hit.hitted = true;
                    selectColor(point, hit.hitColor);
                    return true;
                }
                else{
                    setColorHitNew(hit, point, ray, t_entry);
                    if(hit.hitColor.w == 1.0f){ 
                    return true;}
                }
            }
        }

        
    }

    if(dependsExtra){
        for (int i = currentNode.octreeInfo[2].pointIndex; i < currentNode.octreeInfo[2].pointIndex + currentNode.octreeInfo[2].pointCount; ++i) {
            Point point = getMultiPoint3(i);
            if (intersectSphere(ray, point.position, abs(radiusExtra), t_entry)) {
                if(shadowing)
                {
                    hit.t = t_entry;
                    hit.hitPosition = pointAtParameter(ray, t_entry);
                    hit.hitNormal = normalize(hit.hitPosition - point.position);
                    hit.hitted = true;
                    selectColor(point, hit.hitColor);
                    return true;
                }
                else{
                    setColorHitNew(hit, point, ray, t_entry);
                    if(hit.hitColor.w == 1.0f){ 
                    return true;}
                }
            }
        }

        
    }
    
    return hit.hitted;
    
    
}

bool evaluateMultiNodeWithOrder(FlattenedMultiNode currentNode,Ray ray, float t_min, float t_max,float t_entry, float t_exit, out HitData hit, bool dependsGround, bool dependsError, bool dependsExtra, bool shadowing)
{
    
    // Comprobar si el rayo intersecta con alguno de los puntos
    int numHits = 0;
    if(dependsGround){
        for (int i = currentNode.octreeInfo[0].pointIndex; i < currentNode.octreeInfo[0].pointIndex + currentNode.octreeInfo[0].pointCount; ++i) {
            Point point = getMultiPoint1(i);
            if (intersectSphere(ray, point.position, abs(radius), t_entry)) {
                HitData newHit;
                SphereHitData newSphereHit = SphereHitData(newHit, t_entry, point);
                hits[numHits++] = newSphereHit;
                insertSorted(hits, numHits, newSphereHit);
            }
        }

        
    }
    if(dependsError){
        for (int i = currentNode.octreeInfo[1].pointIndex; i < currentNode.octreeInfo[1].pointIndex + currentNode.octreeInfo[1].pointCount; ++i) {
            Point point = getMultiPoint2(i);
            if (intersectSphere(ray, point.position, abs(radiusError), t_entry)) {
                HitData newHit;
                SphereHitData newSphereHit = SphereHitData(newHit, t_entry, point);
                hits[numHits++] = newSphereHit;
                insertSorted(hits, numHits, newSphereHit);
            }
        }

        
    }

    if(dependsExtra){
        for (int i = currentNode.octreeInfo[2].pointIndex; i < currentNode.octreeInfo[2].pointIndex + currentNode.octreeInfo[2].pointCount; ++i) {
            Point point = getMultiPoint3(i);
            if (intersectSphere(ray, point.position, abs(radiusExtra), t_entry)) {
                HitData newHit;
                SphereHitData newSphereHit = SphereHitData(newHit, t_entry, point);
                hits[numHits++] = newSphereHit;
                insertSorted(hits, numHits, newSphereHit);
            }
        }

        
    }
    if(shadowing && (numHits > 0)){

                hit.t = hits[0].t;
                hit.hitPosition = pointAtParameter(ray, hits[0].t);
                hit.hitNormal = normalize(hit.hitPosition - hits[0].hitPoint.position);
                hit.hitted = true;
                selectColor(hits[0].hitPoint, hit.hitColor);
                return true;
            }
    for (int i = 0; i < numHits; ++i) {
            
            setColorHitNew(hit, hits[i].hitPoint, ray, hits[i].t);
            if(hit.hitColor.w == 1.0f){ 
                    return true;} 
    } 
    return hit.hitted;
    
    
}

bool evaluateMultiNode(FlattenedMultiNode currentNode,Ray ray, float t_min, float t_max,float t_entry, float t_exit, out HitData hit, bool dependsGround, bool dependsError, bool dependsExtra, bool shadowing)
{
    
    if(orderMode == 1){
        return evaluateMultiNodeWithOrder(currentNode, ray, t_min, t_max, t_entry, t_exit, hit, dependsGround, dependsError, dependsExtra, shadowing);
    }
    else{
        return evaluateMultiNodeWithNoOrder(currentNode, ray, t_min, t_max, t_entry, t_exit, hit, dependsGround, dependsError, dependsExtra, shadowing);
    }
    
    
}
bool processMultiOctreeStack(out int stack[MAX_STACK_SIZE], out int stackIndex, Ray ray, float t_min, float t_max,float t_entry, float t_exit, out HitData hit, bool shadowing) {
    
    //Cojemos el indice del nodo que toca
    int currentNodeIndex = stack[--stackIndex];
    //Cojemos el multiNode
    FlattenedMultiNode currentNode = getMultiNode(currentNodeIndex);

    bool depGround = false;
    bool depError = false;
    bool leafGround = false;
    bool leafError = false;
    bool depExtra = false;
    bool leafExtra = false;

    
    //Miramos de que octrees depende
    if(octreeGround==1 && currentNode.octreeInfo[0].depends==1){
        depGround = true;
        if(currentNode.octreeInfo[0].leaf == 1){
            leafGround = true;
            
        }
        //return true;
    }
    

    if(octreeError==1 && currentNode.octreeInfo[1].depends==1){
        depError = true;
        if(currentNode.octreeInfo[1].leaf == 1){
            leafError = true;
        }
        
    }

    if(octreeExtra==1 && currentNode.octreeInfo[2].depends==1){
        depExtra = true;
        if(currentNode.octreeInfo[2].leaf == 1){
            leafExtra = true;
        }
        
    }
     
    if(depGround && !depError && !depExtra){
        if(leafGround){
            return evaluateMultiNode(currentNode, ray, t_min, t_max, t_entry, t_exit, hit, depGround, depError, depExtra, shadowing);
        }
        else
        {
            visitChildren(currentNode, stack, stackIndex, ray, t_min, t_max);
        }
    } 
     if(!depGround && depError && !depExtra){
        if(leafError){
           
            return evaluateMultiNode(currentNode, ray, t_min, t_max, t_entry, t_exit, hit, depGround, depError, depExtra, shadowing);        }
        else
        {
            visitChildren(currentNode, stack, stackIndex, ray, t_min, t_max);
        }
    }
    if(!depGround && !depError && depExtra){
        if(leafExtra){
           
            return evaluateMultiNode(currentNode, ray, t_min, t_max, t_entry, t_exit, hit, depGround, depError, depExtra, shadowing);        }
        else
        {
            visitChildren(currentNode, stack, stackIndex, ray, t_min, t_max);
        }
    }
    if(depGround && depError && !depExtra){
        if(leafGround && leafError){
            return evaluateMultiNode(currentNode, ray, t_min, t_max, t_entry, t_exit, hit, depGround, depError, depExtra, shadowing);        }
        else
        {
            visitChildren(currentNode, stack, stackIndex, ray, t_min, t_max);
        }
    }

    if(depGround && !depError && depExtra){
        if(leafGround && leafExtra){
            return evaluateMultiNode(currentNode, ray, t_min, t_max, t_entry, t_exit, hit, depGround, depError, depExtra, shadowing);        }
        else
        {
            visitChildren(currentNode, stack, stackIndex, ray, t_min, t_max);
        }
    } 
     if(!depGround && depError &&  depExtra){
        if(leafError && leafExtra){
           
            return evaluateMultiNode(currentNode, ray, t_min, t_max, t_entry, t_exit, hit, depGround, depError, depExtra, shadowing);        }
        else
        {
            visitChildren(currentNode, stack, stackIndex, ray, t_min, t_max);
        }
    }
    if(depGround && depError && depExtra){
        if(leafGround && leafError && leafExtra){
            return evaluateMultiNode(currentNode, ray, t_min, t_max, t_entry, t_exit, hit, depGround, depError, depExtra, shadowing);        }
        else
        {
            visitChildren(currentNode, stack, stackIndex, ray, t_min, t_max);
        }
    }

    return hit.hitted;
    


} 


bool searchMultiOctree(Ray ray, float t_min, float t_max, out HitData hit, bool shadowing) {
    int stack[MAX_STACK_SIZE];

    // Agregar el nodo raíz de ambos octrees al inicio de sus respectivas pilas
    int stackIndex = 0;
    stack[stackIndex++] = 0;

    float t_entry, t_exit;
    bool inBox = rayIntersectsMultiNode(ray, getMultiNode(0), t_min, t_max, t_entry, t_exit, max(radius, max(radiusError, radiusExtra)));
    t_max = t_exit;
    bool foundHit = false;
    
    if(!inBox)
    {
        return false;
    }
   while (stackIndex > 0 && hit.hitColor.w < 1.f) { 
        processMultiOctreeStack(stack, stackIndex, ray, t_min, t_max, t_entry, t_exit, hit, shadowing);
    }
    
    return hit.hitted;
}

bool world_hit(Ray r, float t_min, float t_max, out HitData hit) {
  HitData temp_hit;
  bool hit_anything = false;
  float closest_so_far = t_max;

  if(searchMultiOctree(r, t_min, closest_so_far, temp_hit, false))
  {
    hit_anything = true;
    hit = temp_hit;
    closest_so_far = temp_hit.t;
  }
  
  if(planesMode == 1){
    if (plane_hit(r, t_min, closest_so_far,planeOrigin,planeSize,range, temp_hit)) {
        hit_anything = true;
        hit = temp_hit;
    }
  }

  

  return hit_anything;
}

bool world_has_hit(Ray r, float t_min, float t_max) {
  HitData temp_hit;
  bool hit_anything = false;
  float closest_so_far = t_max;
  if(searchMultiOctree(r, t_min, closest_so_far, temp_hit,true))
  {
	  return true;
  }
  if(planesMode == 1){
      if (plane_hit(r, t_min, closest_so_far,planeOrigin,planeSize,range, temp_hit)) {
	return true;
    }
  }


  return false;
}

vec3 computeShadow(PointLight light, vec3 p, int maxT, bool soft){ 
    vec3 shadow = vec3(0.f,0.f,0.f);
	float eps = 0.05;
	
	vec3 kt;
	HitData hitData;
	float u;
	float v;
	vec3 newLightPosition;
	vec3 dir;
	float dmax;
	vec3 lightWidth = vec3(1.f,0.f,0.f);
	vec3 lighHeight = vec3(0.f,1.f,0.f);
	Ray r;

	if (!soft)
	{
		 dir = normalize(light.position - p);
		 dmax = distance(p, light.position);
		 Ray r = Ray(p, dir);

		if(world_has_hit(r, eps, dmax))
		{
			return vec3(0.f,0.f,0.f);
		}
		return vec3(1.f,1.f,1.f);

	}
	
	for (int i = 0; i< maxT; i++)
	{
		u = drand48(light.position.xy + i) / 2;
      	v = drand48(light.position.xz + i) / 2;
		newLightPosition = light.position + u * lightWidth + v * lighHeight;

		dir = normalize(newLightPosition - p);
		
		dmax =  distance(p, newLightPosition);
		r = Ray(p, dir);
		
		if(world_has_hit(r, eps, dmax) == false)
		{   
            //return vec3(0.f,0.f,0.f);
			shadow += vec3(1.f,1.f,1.f);
		}
		
	}

	

	return shadow / maxT;
    
} 

vec4 phong(Ray r, HitData hitData, bool shadowActive, bool localOcclusionActive, bool globalOcclusionActive, vec3 kd, vec3 ks, vec3 ka)
{
	
	//Inicialitzem les llums
	vec3 ambientFinal = globalLight;// * globalOccRatio;
	vec3 diffuseFinal = vec3(0.f);
	vec3 specularFinal = vec3(0.f);

	
	//Vector Hit to LlumPosition.
    for(int i=0; i<nLights; i++)
    {
        PointLight pointlight = pointLights[i];
        vec3 L = normalize(pointlight.position - hitData.hitPosition);
        //Invers de la direccio del raig.
        vec3 V = -r.direction;


        //Attenuation
        float distlightHit = length(pointlight.position - hitData.hitPosition);
        //Constant linear quadratic
        float attenuation = 1.f / (pointlight.constant + pointlight.linear * distlightHit + pointlight.quadratic * (distlightHit * distlightHit));

        //SHADOW
        vec3 shadow;
        if(shadowActive && shadowMode==1)
        {
            //return vec4(1.f,1.f,1.f,1.f);
            shadow = computeShadow(pointlight,hitData.hitPosition,128,false);
        }
        else{
            shadow = vec3(1.f,1.f,1.f);
        }


        //DIFUSA
        float lightFactor = dot(L, hitData.hitNormal);
        if(lightFactor < 0) {lightFactor = 0;}
        vec3 Id = pointlight.Id;
        vec3 Kd = kd;
        diffuseFinal += Id*Kd*lightFactor*attenuation*shadow;

        //SPECULAR
        vec3 Is = pointlight.Is;
        vec3 Ks = ks;
        if ((abs(Kd.x-Kd.y) < 0.2) && (Kd.z > Kd.x) && (Kd.z > Kd.y)){
                Ks += vec3(0.1, 0.1, 0.1);
        }
        vec3 H = normalize(L+V);
        lightFactor = dot(hitData.hitNormal,H);
        if(lightFactor<0){lightFactor=0;}
        lightFactor = pow(lightFactor,0.1f);
        specularFinal += Ks*Is*lightFactor*attenuation*shadow;

        //Ambient
        vec3 Ia = pointlight.Ia;
        vec3 Ka = ka;
        ambientFinal += Ia*Ka*attenuation;
    }
        
	

	vec3 col  = diffuseFinal + specularFinal+ambientFinal;

	return vec4(col,hitData.hitColor.w);
}

vec4 color(Ray r, int depth)
{
     
	//Visible color

	 HitData hitData;
    vec4 color = vec4(0.f,0.f,0.f,0.f);
	if (world_hit(r, 0.0001f, uMaxRayDistance, hitData))
	{
		
		
        if(hitData.hitted == true)
        {   
            vec3 kd = material.diffuse;
            vec3 ks = material.specular;
            vec3 ka = material.ambient;
            if(hitData.plane){}
            else
            {
                kd = hitData.hitColor.xyz;
                ks = hitData.hitColor.xyz;
                ka = hitData.hitColor.xyz;
            }
            if(shadingMode == 1){
                return phong(r, hitData, hitData.plane, false, false,kd,ks,ka);
            }
            else{
                return vec4(hitData.hitColor.xyz,1.f);
            }
            
        }
        
        else
        {
            return vec4(material.ambient,0.1f);
        }
    }
}

Point pointClick(Ray r, int depth)
{
     
	//Visible color

	 HitData hitData;
    vec4 color = vec4(0.f,0.f,0.f,0.f);
	if (world_hit(r, 0.0001f, uMaxRayDistance, hitData))
	{
		return hitData.hitPoint;
    }
}

void main()
{
  
    
	vec4 col = vec4(0.f,0.f,0.f,0.f);

	float u, v;
    Ray r;

	const int nsamples = 1;
    for (int s = 0; s < nsamples; s++) {
      u = ((gl_FragCoord.x + drand48(col.xy + s)) / window_size.x) * 2.0 - 1.0;
      v = ((gl_FragCoord.y + drand48(col.xz + s)) / window_size.y) * 2.0 - 1.0;
      r = get_ray(vec2(u, v));
      col += color(r,3);
    }
    col /= nsamples;
    col = vec4(sqrt(col.x),sqrt(col.y),sqrt(col.z),sqrt(col.w));
    
    fs_color = col;
}
