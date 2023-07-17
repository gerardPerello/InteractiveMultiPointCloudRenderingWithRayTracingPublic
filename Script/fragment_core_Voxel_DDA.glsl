#version 440

#define RAND_MAX 0x7fff

struct Ray
{
	vec3 origin;
	vec3 direction;	
	//vec4 colorIntensity;
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
};

struct FlattenedNode {
  vec3 minBound;
  vec3 maxBound;
  ivec4 childrenInfo1; 
  ivec4 childrenInfo2; 
  int pointIndex;
  int pointCount;
};

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


//Se usa
in vec3 vs_position;

out vec4 fs_color;
//Uniforms
uniform Material material;
//uniform PointLight pointlight;
#define MAX_POINT_LIGHTS 3

uniform PointLight pointLights[MAX_POINT_LIGHTS];
uniform int nLights;
//Uniforms
uniform vec3 camPosition;
uniform mat4 ViewMatrix;
uniform mat4 ProjectionMatrix;
uniform sampler3D texture3DD;
uniform vec3 uVoxelDimensions;
uniform vec3 uMinVertex;
uniform vec3 uMaxVertex;
uniform float uEpsilon;
uniform float uMaxRayDistance;
uniform vec2 window_size;
uniform vec3 globalLight;
uniform vec2 randomSeed;
uniform vec3 planeOrigin;
uniform float planeSize;
uniform vec3 range;
uniform ivec3 grid_dim;

//UnfiromsMode
uniform int mode;
uniform int shadowMode;
uniform int orderMode;
uniform int shadingMode;
uniform int planesMode;


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

vec4 mixColor(vec3 color1, float alfa1, vec3 color2, float alfa2)
{
    
    float newAlpha = min(1.0f, alfa1 + alfa2 * (1.0f - alfa1));
    vec3 newColor = (color1 * alfa1 + color2 * alfa2 * (1.0f - alfa1)) / newAlpha;
    return vec4(newColor, newAlpha);
}

// Función pointAtParameter: devuelve un punto en el rayo a una distancia t desde el origen
vec3 pointAtParameter(Ray ray, float t) {
  return ray.origin + t * ray.direction;
}

// Función get_ray: crea un rayo desde la posición de la cámara que pasa a través de un punto en la pantalla
Ray get_ray(vec2 screenPos) {
  // Convierte la posición de la pantalla a coordenadas de clip
  vec4 ray_clip = vec4(screenPos, -1.0, 1.0);
  
  // Convierte las coordenadas de clip a coordenadas de espacio ocular
  vec4 ray_eye = (ProjectionMatrix) * ray_clip;
  ray_eye = vec4(ray_eye.xy, -1.0, 0.0);
  
  // Convierte las coordenadas del espacio ocular a coordenadas del espacio mundial
  vec3 ray_wor = ((ViewMatrix) * ray_eye).xyz;
  
  // Normaliza la dirección del rayo y agrega un pequeño valor para evitar problemas numéricos
  ray_wor = normalize(ray_wor + vec3(0.0001, 0.0001, 0));
  
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

vec3 calculateNormal(vec3 stepDir, vec3 t_next) {
    vec3 normal = vec3(0.0, 0.0, 0.0);

    if (t_next.x < t_next.y) {
        if (t_next.x < t_next.z) {
            normal.x = -stepDir.x;
        } else {
            normal.z = -stepDir.z;
        }
    } else {
        if (t_next.y < t_next.z) {
            normal.y = -stepDir.y;
        } else {
            normal.z = -stepDir.z;
        }
    }

    return normal;
}

bool slabs(vec3 boxMin, vec3 boxMax, Ray ray, float t_min, float t_max, out float t_entry, out float t_exit, float epsilon) {
    vec3 expandedMinBound = boxMin - vec3(epsilon, epsilon, epsilon);
    vec3 expandedMaxBound = boxMax + vec3(epsilon, epsilon, epsilon);

    vec3 rayEntry = (expandedMinBound - ray.origin) / ray.direction;
    vec3 rayExit = (expandedMaxBound - ray.origin) / ray.direction;

    vec3 t_min_v = min(rayEntry, rayExit);
    vec3 t_max_v = max(rayEntry, rayExit);

    t_entry = max(max(t_min_v.x, t_min_v.y), max(t_min_v.z, t_min));
    t_exit = min(min(t_max_v.x, t_max_v.y), min(t_max_v.z, t_max));

    return t_entry <= t_exit;
}

bool hit_texture(Ray ray, float t_min, float t_max, out HitData hit)
{
	
  vec3 rayDirection = normalize(ray.direction);

	// Cálculo de delta para avanzar en el rayo
	vec3 delta = uVoxelDimensions * (1.0 / abs(rayDirection));

  float t_entry;
  float t_exit;
  slabs(uMinVertex, uMaxVertex, ray, t_min, t_max, t_entry, t_exit, uEpsilon);

	// Verifica si el rayo no intersecta el volumen, la bounding box inicial
  if (t_entry > t_exit) {
      return false;
  }
	
	// Asegura que el tiempo de entrada esté dentro del rango permitido
    float t_current = t_entry;

	 // Encuentra el voxel actual y el paso en cada eje
    vec3 currentVoxel = floor((ray.origin + t_current * rayDirection - uMinVertex) / uVoxelDimensions);
    vec3 stepDir = sign(rayDirection);
    
    // Encuentra el siguiente límite del voxel y el tiempo de llegada al mismo
    vec3 nextBoundary = (currentVoxel + stepDir * 0.5) * uVoxelDimensions + uMinVertex;
    vec3 t_next = (nextBoundary - ray.origin) / rayDirection;

	vec3 intersectPoint;
	vec3 textureCoords;
	vec4 textureColor;

	// Itera sobre los voxels atravesados por el rayo
	while (t_exit >= t_current && t_current <= t_max && hit.hitColor.w < 1.f)
	{
		// Calcula el punto de intersección y las coordenadas de textura
		intersectPoint = pointAtParameter(ray, t_current);
		textureCoords =(currentVoxel + 0.5) / grid_dim; //((intersectPoint - uMinVertex) / (uMaxVertex - uMinVertex));
		textureColor = texture3D(texture3DD, textureCoords);
    
    if (textureColor.w > 0)
		{
        hit.t = t_current;
        if(hit.hitted == false)
        {
          hit.hitPosition = intersectPoint;
          hit.hitNormal = calculateNormal(stepDir, t_next);
          hit.hitted = true;
          hit.hitColor = vec4(textureColor);//vec4(maskValue.r,0.f,0.f,1.f);
          
        }
        mixColor(hit.hitColor.xyz, hit.hitColor.w, textureColor.xyz, textureColor.w);
    
		}	
		//Cual de los siguientes limites se determinara primero? Vamos a ese.
    if (t_next.x < t_next.y) {
        if (t_next.x < t_next.z) {
            t_current = t_next.x;
            t_next.x += delta.x;
            currentVoxel.x += stepDir.x;
        } else {
            t_current = t_next.z;
            t_next.z += delta.z;
            currentVoxel.z += stepDir.z;
        }
    } else {
        if (t_next.y < t_next.z) {
            t_current = t_next.y;
            t_next.y += delta.y;
            currentVoxel.y += stepDir.y;
        } else {
            t_current = t_next.z;
            t_next.z += delta.z;
            currentVoxel.z += stepDir.z;
        }
    }
	}
	return hit.hitted;

}


/* Check all objects in world for hit with ray */
bool world_hit(Ray r, float t_min, float t_max, out HitData hit) {
  HitData temp_hit;
  bool hit_anything = false;
  float closest_so_far = t_max;
  if(hit_texture(r, t_min, closest_so_far, temp_hit))
  {
    hit_anything = true;
    hit = temp_hit;
    closest_so_far = temp_hit.t;
  }
  if (planesMode == 1){
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
  if(hit_texture(r, t_min, closest_so_far, temp_hit))
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

vec4 color(Ray r)
{
     
	  HitData hitData;
    vec4 color = vec4(0.f,0.f,0.f,0.f);
	if (world_hit(r, 0.0001f, uMaxRayDistance, hitData))
	{
		
		
        if(hitData.hitted == true)
        {   
            vec3 kd = material.diffuse;
            vec3 ks = material.specular;
            vec3 ka = material.ambient;
            if(hitData.plane){
                
            }
                
            
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
      col += color(r);
    }
    col /= nsamples;
    col = vec4(sqrt(col.x),sqrt(col.y),sqrt(col.z),sqrt(col.w));

    fs_color = col;
    //fs_color = vec4(point.x,point.y, point.z, 1.f);
}
