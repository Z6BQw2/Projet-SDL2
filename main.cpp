#include <iostream>
#include <vector>
#include <cmath>
#include "shapes.h"     
#include "sdl_class.h"

// Taille fenetre
#define W 800
#define H 600

// Projection 3D vers 2D
Point2d project(Point3D<float> p, Point3D<float> cam) {
    float x = p.x - cam.x;
    float y = p.y - cam.y;
    float z = p.z - cam.z;

    // Focale arbitraire
    float f = 600.0f;

    // Formule de projection
    int sx = (int)((x * f) / z) + (W / 2);
    int sy = (int)((y * f) / z) + (H / 2);

    return Point2d(sx, sy);
}

// Fonction pour dessiner une face
void draw_face(Sdl& sdl, Quad3D<float> q, Point3D<float> cam, int r, int g, int b) {
    // On recupere les 2 triangles
    Triangle3D<float> tris[2];
    tris[0] = q.get_T1();
    tris[1] = q.get_T2();

    for (int i = 0; i < 2; i++) {
        Triangle3D<float> t = tris[i];

        // Calcul normale
        Point3D<float> v1 = t.p2 - t.p1;
        Point3D<float> v2 = t.p3 - t.p1;
        Point3D<float> norm = cross_product(v1, v2);

        // Backface culling
        Point3D<float> cam_vec = t.p1 - cam;
        if (scal(norm, cam_vec) <= 0) continue; 

        // Lumiere (flat shading)
        norm.norm();
        Point3D<float> light(0.2f, 0.5f, -1.0f); // direction lumiere
        light.norm();

        float inten = scal(norm, light);
        if (inten < 0.2f) inten = 0.2f;
        if (inten > 1.0f) inten = 1.0f;

        // Couleur finale
        sdl.set_color(r * inten, g * inten, b * inten);

        // Projection des points
        Point2d p1 = project(t.p1, cam);
        Point2d p2 = project(t.p2, cam);
        Point2d p3 = project(t.p3, cam);

        // Remplissage
        sdl.draw_filled_triangle(p1, p2, p3);
    }
}

int main(int argc, char* argv[]) {
    // Init SDL
    Sdl sdl("Moteur 3D", W, H);
    
    // Scene
    Scene3D<float> scene;

    // Ajout objets
    Pave3D<float> cube(Point3D<float>(-1.0f, -1.0f, 0.0f), 2.0f, 2.0f, 2.0f);
    scene.add_pave(cube);

    Sphere3D<float> sphere(Point3D<float>(3.0f, 0.0f, 0.0f), 1.5f, 15, 15);
    scene.add_sphere(sphere);

    // Camera
    scene.get_camera() = Point3D<float>(0, 0, -5);

    bool running = true;
    
    while (running) {
        // Events
        if (sdl.poll_quit()) running = false;

        // Animation
        std::vector<Pave3D<float>>& paves = scene.get_paves_mut();
        for(int i = 0; i < paves.size(); i++) {
            paves[i].rotate(Point3D<float>(0,0,1), 0.01f, 0.02f, 0.0f);
        }

        std::vector<Sphere3D<float>>& spheres = scene.get_spheres_mut();
        for(int i = 0; i < spheres.size(); i++) {
            spheres[i].rotate(Point3D<float>(3,0,0), 0.00f, 0.03f, 0.0f);
        }

        // Rendu
        sdl.clear();
        Point3D<float> cam = scene.get_camera();

        // Dessin des paves
        for (int i = 0; i < paves.size(); i++) {
            std::vector<Quad3D<float>> faces = paves[i].get_faces();
            for(int j = 0; j < faces.size(); j++) {
                draw_face(sdl, faces[j], cam, 0, 255, 0); // Vert
            }
        }

        // Dessin des spheres
        for (int i = 0; i < spheres.size(); i++) {
            std::vector<Quad3D<float>> faces = spheres[i].get_faces();
            for(int j = 0; j < faces.size(); j++) {
                draw_face(sdl, faces[j], cam, 0, 100, 255); // Bleu
            }
        }

        sdl.present();
        SDL_Delay(16);
    }

    return 0;
}