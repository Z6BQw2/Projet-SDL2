#include <iostream>
#include <vector>
#include <cuda_runtime.h>
#include "sdl_class.h"

// Config
#define W 800
#define H 600
#define INF 1e20f

// Structures compatibles CPU/GPU pour simplifier le transfert
struct Vec3 {
    float x, y, z;
    
    __host__ __device__ Vec3 operator+(const Vec3& v) const { return {x+v.x, y+v.y, z+v.z}; }
    __host__ __device__ Vec3 operator-(const Vec3& v) const { return {x-v.x, y-v.y, z-v.z}; }
    __host__ __device__ Vec3 operator*(float f) const { return {x*f, y*f, z*f}; }
    
    __host__ __device__ float dot(const Vec3& v) const { return x*v.x + y*v.y + z*v.z; }
    
    __host__ __device__ Vec3 normalize() const {
        float l = sqrtf(x*x + y*y + z*z);
        return {x/l, y/l, z/l};
    }
};

struct Ray {
    Vec3 orig;
    Vec3 dir;
};

struct CudaSphere {
    Vec3 center;
    float radius;
    Vec3 color;
    bool is_light; 
};

// Intersection simple Rayon-Sphere
__device__ float hit_sphere(const CudaSphere& s, const Ray& r) {
    Vec3 oc = r.orig - s.center;
    float a = r.dir.dot(r.dir);
    float b = 2.0f * oc.dot(r.dir);
    float c = oc.dot(oc) - s.radius * s.radius;
    float delta = b*b - 4*a*c;
    
    if (delta < 0) return -1.0f;
    return (-b - sqrtf(delta)) / (2.0f*a);
}

__global__ void render_kernel(uchar4* buffer, CudaSphere* spheres, int nb_spheres, float time) {
    int x = threadIdx.x + blockIdx.x * blockDim.x;
    int y = threadIdx.y + blockIdx.y * blockDim.y;
    
    if (x >= W || y >= H) return;

    // Setup rayon caméra (fixe)
    float u = ((float)x / W * 2.0f - 1.0f) * ((float)W / H);
    float v = ((float)y / H * 2.0f - 1.0f);

    Ray r = { {0.0f, 0.0f, -20.0f}, Vec3{u, v, 1.0f}.normalize() };

    // Hack: Copie locale pour animation directe dans le kernel
    // Évite de renvoyer les données depuis le CPU à chaque frame
    CudaSphere local_s[3]; 
    for(int i=0; i<3; i++) local_s[i] = spheres[i];

    // Animation orbites
    local_s[1].center.x = cos(time) * 8.0f;
    local_s[1].center.z = sin(time) * 8.0f;

    local_s[2].center.x = cos(time * 1.5f + 1.0f) * 12.0f;
    local_s[2].center.z = sin(time * 1.5f + 1.0f) * 12.0f;

    // Recherche objet le plus proche
    float t_min = INF;
    int id = -1;

    for (int i = 0; i < nb_spheres; ++i) {
        float t = hit_sphere(local_s[i], r);
        if (t > 0.0f && t < t_min) {
            t_min = t;
            id = i;
        }
    }

    Vec3 col = {0.0f, 0.0f, 0.0f};

    if (id != -1) {
        if (local_s[id].is_light) {
            col = local_s[id].color;
        } else {
            // Calcul éclairage (Phong simple)
            Vec3 p = r.orig + r.dir * t_min;
            Vec3 n = (p - local_s[id].center).normalize();
            
            // Light source = obj 0
            Vec3 light_pos = local_s[0].center;
            Vec3 l_dir = (light_pos - p).normalize();

            // Shadow ray
            bool shadow = false;
            Ray s_ray = { p + n * 0.01f, l_dir };
            
            // Check collision vers la lumière
            float dist_light = (light_pos - p).dot(l_dir);
            for(int j=0; j<nb_spheres; j++) {
                if (j == id || local_s[j].is_light) continue;
                
                float t = hit_sphere(local_s[j], s_ray);
                if (t > 0.0f && t < dist_light) {
                    shadow = true;
                    break;
                }
            }

            if (!shadow) {
                float diff = max(0.0f, n.dot(l_dir));
                col = local_s[id].color * diff;
            }
        }
    }

    // Écriture texture
    int idx = y * W + x;
    buffer[idx] = make_uchar4(
        (unsigned char)min(255.0f, col.x * 255.0f),
        (unsigned char)min(255.0f, col.y * 255.0f),
        (unsigned char)min(255.0f, col.z * 255.0f),
        255
    );
}

int main() {
    Sdl engine("Bonus CUDA Raytracing", W, H);

    // Setup scène
    std::vector<CudaSphere> h_spheres;
    
    // 0: Soleil
    h_spheres.push_back({{0.0f, 0.0f, 0.0f}, 4.0f, {1.0f, 0.6f, 0.1f}, true});
    // 1: Planete 1
    h_spheres.push_back({{0.0f, 0.0f, 0.0f}, 1.5f, {1.0f, 0.1f, 0.1f}, false});
    // 2: Planete 2
    h_spheres.push_back({{0.0f, 0.0f, 0.0f}, 2.0f, {0.1f, 0.8f, 0.8f}, false});

    // Alloc buffers GPU
    uchar4* d_buffer;
    cudaMalloc(&d_buffer, W * H * sizeof(uchar4));

    CudaSphere* d_spheres;
    cudaMalloc(&d_spheres, h_spheres.size() * sizeof(CudaSphere));
    
    // Transfert initial (Animation gérée côté GPU)
    cudaMemcpy(d_spheres, h_spheres.data(), h_spheres.size() * sizeof(CudaSphere), cudaMemcpyHostToDevice);

    uchar4* h_buffer = (uchar4*)malloc(W * H * sizeof(uchar4));

    // Config kernel
    dim3 block(16, 16);
    dim3 grid((W + 15) / 16, (H + 15) / 16);

    float time = 0.0f;
    bool running = true;

    while (running) {
        if (engine.poll_quit()) running = false;

        render_kernel<<<grid, block>>>(d_buffer, d_spheres, h_spheres.size(), time);
        
        // Synchro avant copie
        cudaDeviceSynchronize();
        cudaMemcpy(h_buffer, d_buffer, W * H * sizeof(uchar4), cudaMemcpyDeviceToHost);

        engine.clear();
        engine.draw_buffer(h_buffer, W * sizeof(uchar4));
        engine.present();

        time += 0.01f;
        SDL_Delay(16);
    }

    // Cleanup
    cudaFree(d_buffer);
    cudaFree(d_spheres);
    free(h_buffer);

    return 0;
}