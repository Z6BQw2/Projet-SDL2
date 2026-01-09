#ifndef SHAPE_H
#define SHAPE_H

#include <cmath>
#include <iostream>
#include <vector>
#include <stdexcept>

/// @brief Classe Template pour un point 3D.
template <typename T>
class Point3D {
public: 
    T x, y, z;

    /// @brief Constructeurs.
    Point3D() : x(0), y(0), z(0) {}
    Point3D(T val) : x(val), y(val), z(val) {}
    Point3D(T x, T y) : x(x), y(y), z(0) {}
    Point3D(T f, T s, T t) : x(f), y(s), z(t) {}

    /// @brief Opérateur +=.
    Point3D<T>& operator+=(const Point3D& p){
        x += p.x; y += p.y; z += p.z;
        return *this;
    }

    /// @brief Opérateur -=.
    Point3D<T>& operator-=(const Point3D& p){
        x -= p.x; y -= p.y; z -= p.z;
        return *this; 
    }

    /// @brief Produit scalaire.
    T scal(const Point3D& p) const {
        return x*p.x + y*p.y + z*p.z;
    }

    /// @brief Normalisation du vecteur.
    Point3D<T>& norm(){
        T val = std::sqrt(scal(*this));
        if (val == 0) throw std::invalid_argument("Norme nulle");
        x /= val; y /= val; z /= val;
        return *this;
    }

    /// @brief Rotation et Translation combinées (Code origine).
    template <typename T2>
    Point3D<T>& translate_rotate(const Point3D<T>& center, const T2& angle_x, const T2& angle_y, const T2& angle_z){
        T tx = x - center.x;
        T ty = y - center.y;
        T tz = z - center.z;

        // Rotation X
        T temp_y = ty * std::cos(angle_x) - tz * std::sin(angle_x);
        T temp_z = ty * std::sin(angle_x) + tz * std::cos(angle_x);
        ty = temp_y; tz = temp_z;

        // Rotation Y
        T temp_x = tx * std::cos(angle_y) + tz * std::sin(angle_y);
        temp_z = -tx * std::sin(angle_y) + tz * std::cos(angle_y);
        tx = temp_x; tz = temp_z;

        // Rotation Z
        temp_x = tx * std::cos(angle_z) - ty * std::sin(angle_z);
        temp_y = tx * std::sin(angle_z) + ty * std::cos(angle_z);
        tx = temp_x; ty = temp_y;
        
        x = tx + center.x;
        y = ty + center.y;
        z = tz + center.z;
        return *this;
    }

    /// @brief Affichage console.
    void print_point(){
        std::cout << "x: " << x << " y: " << y << " z: " << z << std::endl;
    }
}; 

/// @brief Addition de points.
template <typename T>
Point3D<T> operator+(const Point3D<T>&p, const Point3D<T>&p2){
    Point3D<T> res = p; res += p2; return res;
}

/// @brief Multiplication scalaire.
template <typename T>
Point3D<T> operator*(const Point3D<T>&p, T val){
    Point3D<T> res = p;
    res.x *= val; res.y *= val; res.z *= val;
    return res;
}

/// @brief Soustraction.
template <typename T>
Point3D<T> operator-(const Point3D<T>&p, const Point3D<T>&p2){
    Point3D<T> res = p; res -= p2; return res;
}

/// @brief Wrapper produit scalaire.
template <typename T>
T scal(const Point3D<T>&p, const Point3D<T>&p2){ return p.scal(p2); }

/// @brief Égalité avec epsilon.
template <typename T>
bool operator==(const Point3D<T>&p, const Point3D<T>&p2){
    return (std::abs(p.x - p2.x) < T(1e-5) && std::abs(p.y - p2.y) < T(1e-5) && std::abs(p.z - p2.z) < T(1e-5));
}

/// @brief Inégalité.
template <typename T>
bool operator!=(const Point3D<T>&p, const Point3D<T>&p2){ return !(p == p2); }

/// @brief Produit vectoriel (Cross Product).
template <typename T>
Point3D<T> cross_product(const Point3D<T>& a, const Point3D<T>& b) {
    return Point3D<T>(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}

/// @brief Classe Triangle (3 points).
template <typename T>
class Triangle3D{
public:
    Point3D<T> p1, p2, p3;

    /// @brief Constructeurs.
    Triangle3D() : p1(), p2(Point3D<T>{T(0), T(2), T(0)}), p3(Point3D<T>{T(0), T(0), T(2)}) {}
    Triangle3D(Point3D<T> val) : p1(val), p2(val), p3(val) {}
    Triangle3D(Point3D<T> x, Point3D<T> y) : p1(x), p2(y), p3() {}
    Triangle3D(Point3D<T> f, Point3D<T> s, Point3D<T> t) : p1(f), p2(s), p3(t) {}

    /// @brief Centre de gravité.
    Point3D<T> center(){
        Point3D<T> G = p1 + p2 + p3;
        G.x /= 3; G.y /= 3; G.z /= 3;
        return G;
    }

    /// @brief Translation.
    Triangle3D<T>& translation(const Point3D<T>& p){
        p1 += p; p2 += p; p3 += p;
        return *this;
    }

    /// @brief Rotation.
    template <typename T2>
    Triangle3D<T>& rotate(const T2& angle_x, const T2& angle_y, const T2& angle_z){
        Point3D<T> center = this->center();
        p1 = p1.translate_rotate(center, angle_x, angle_y, angle_z);
        p2 = p2.translate_rotate(center, angle_x, angle_y, angle_z);
        p3 = p3.translate_rotate(center, angle_x, angle_y, angle_z);
        return *this;
    }

    /// @brief Redimensionnement.
    Triangle3D<T>& enlarge_reduce(T diff){
        Point3D<T> center = this->center();
        p1 += center + (p1 - center) * diff;
        p2 += center + (p2 - center) * diff;
        p3 += center + (p3 - center) * diff;
        return *this;
    }
};

/// @brief Classe Quad (2 triangles).
template <typename T>
class Quad3D{
private:
    Triangle3D<T> T1, T2;

public:
    /// @brief Constructeurs.
    Quad3D() : T1(), T2(){}
    Quad3D(const Triangle3D<T>& Tr) : T1(Tr), T2(Triangle3D<T>{Tr.p1, Tr.p2, Tr.p1 + Tr.p2 - Tr.p3}){}
    Quad3D(const Triangle3D<T>& Tr1, const Triangle3D<T>& Tr2) : T1(Tr1), T2(Tr2){}

    /// @brief Getters.
    const Triangle3D<T>& get_T1() const { return T1; }
    Triangle3D<T>& get_T1() { return T1; }
    const Triangle3D<T>& get_T2() const { return T2; }
    Triangle3D<T>& get_T2() { return T2; }

    /// @brief Déplacement d'un point spécifique.
    Quad3D<T>& move_point(T val, bool p1 = true, bool p2 = false, bool p3 = false){
        if(p1) { T1.p1 += val; T2.p1 += val; } 
        else if(p2) { T1.p2 += val; T2.p2 += val; } 
        else if(p3) { T1.p3 += val; } 
        else { T2.p3 += val; }
        return *this;
    }

    /// @brief Extension.
    Quad3D<T>& extend(T diff){
        T1 = T1.enlarge_reduce(diff);
        T2 = T2.enlarge_reduce(diff);
        return *this;
    }

    /// @brief Rotation.
    template <typename T_2>
    Quad3D<T>& rotate(Point3D<T> p, const T_2& angle_x = 0, const T_2& angle_y = 0, const T_2& angle_z = 0){
        T1.p1.translate_rotate(p, angle_x, angle_y, angle_z);
        T1.p2.translate_rotate(p, angle_x, angle_y, angle_z);
        T1.p3.translate_rotate(p, angle_x, angle_y, angle_z);
        T2.p1.translate_rotate(p, angle_x, angle_y, angle_z);
        T2.p2.translate_rotate(p, angle_x, angle_y, angle_z);
        T2.p3.translate_rotate(p, angle_x, angle_y, angle_z);
        return *this;
    }
};

/// @brief Classe Pavé (Cube).
template <typename T>
class Pave3D{
private:
    std::vector<Quad3D<T>> faces;

public:
    /// @brief Constructeur par défaut.
    Pave3D() : Pave3D(Point3D<T>(0, 0, 0), 1, 1, 1) {}

    /// @brief Constructeur manuel (Code origine).
    Pave3D(Point3D<T> origine, T largeur, T hauteur, T profondeur) {
        Point3D<T> p000 = origine;
        Point3D<T> p100 = origine + Point3D<T>(largeur, 0, 0);
        Point3D<T> p010 = origine + Point3D<T>(0, hauteur, 0);
        Point3D<T> p001 = origine + Point3D<T>(0, 0, profondeur);
        Point3D<T> p110 = origine + Point3D<T>(largeur, hauteur, 0);
        Point3D<T> p101 = origine + Point3D<T>(largeur, 0, profondeur);
        Point3D<T> p011 = origine + Point3D<T>(0, hauteur, profondeur);
        Point3D<T> p111 = origine + Point3D<T>(largeur, hauteur, profondeur);

        // Faces construites manuellement
        faces.push_back(Quad3D<T>(Triangle3D<T>(p000, p100, p110), Triangle3D<T>(p000, p110, p010)));
        faces.push_back(Quad3D<T>(Triangle3D<T>(p101, p001, p011), Triangle3D<T>(p101, p011, p111)));
        faces.push_back(Quad3D<T>(Triangle3D<T>(p001, p000, p010), Triangle3D<T>(p001, p010, p011)));
        faces.push_back(Quad3D<T>(Triangle3D<T>(p100, p101, p111), Triangle3D<T>(p100, p111, p110)));
        faces.push_back(Quad3D<T>(Triangle3D<T>(p010, p110, p111), Triangle3D<T>(p010, p111, p011)));
        faces.push_back(Quad3D<T>(Triangle3D<T>(p001, p101, p100), Triangle3D<T>(p001, p100, p000)));
    }

    /// @brief Accesseurs faces.
    const std::vector<Quad3D<T>>& get_faces() const { return faces; }
    std::vector<Quad3D<T>> get_faces_copy() const { return faces; }

    /// @brief Translation.
    void translate(const Point3D<T>& vec) {
        for (auto& face : faces) {
            face.get_T1().translation(vec);
            face.get_T2().translation(vec);
        }
    }

    /// @brief Rotation.
    template <typename T2>
    void rotate(Point3D<T> center, const T2& ax, const T2& ay, const T2& az) {
        for (auto& face : faces) {
            face.rotate(center, ax, ay, az);
        }
    }
    
    /// @brief Centre moyen.
    Point3D<T> center() {
        Point3D<T> sum(0, 0, 0);
        for(auto& face : faces) {
            sum += face.get_T1().center(); sum += face.get_T2().center();
        }
        return sum * (T(1) / T(12)); 
    }
};

/// @brief Classe Sphère.
template <typename T>
class Sphere3D {
private:
    std::vector<Quad3D<T>> faces; 
    Point3D<T> center_pos; 
    T radius; 

    /// @brief Helper coordonnées sphériques.
    Point3D<T> compute_point(int i, int j, int rings, int sectors) const {
        double r_step = (double)i / rings;
        double s_step = (double)j / sectors;
        double theta = r_step * 3.1415926535;
        double phi = s_step * 2 * 3.1415926535;

        T x = radius * std::sin(theta) * std::cos(phi);
        T y = radius * std::cos(theta); 
        T z = radius * std::sin(theta) * std::sin(phi);

        return center_pos + Point3D<T>(x, y, z);
    }

public:
    /// @brief Constructeurs.
    Sphere3D() : Sphere3D(Point3D<T>(0, 0, 0), 1.0) {}
    Sphere3D(Point3D<T> c, T r, int rings = 12, int sectors = 12) : center_pos(c), radius(r) {
        for (int i = 0; i < rings; ++i) {
            for (int j = 0; j < sectors; ++j) {
                Point3D<T> p1 = compute_point(i, j, rings, sectors);
                Point3D<T> p2 = compute_point(i + 1, j, rings, sectors);
                Point3D<T> p3 = compute_point(i + 1, j + 1, rings, sectors);
                Point3D<T> p4 = compute_point(i, j + 1, rings, sectors);
                faces.push_back(Quad3D<T>(Triangle3D<T>(p1, p2, p3), Triangle3D<T>(p1, p3, p4)));
            }
        }
    }

    /// @brief Getters.
    const std::vector<Quad3D<T>>& get_faces() const { return faces; }
    std::vector<Quad3D<T>>& get_faces_mut() { return faces; }

    /// @brief Translation.
    void translate(const Point3D<T>& v) {
        center_pos += v;
        for (auto& face : faces) {
            face.get_T1().translation(v);
            face.get_T2().translation(v);
        }
    }

    /// @brief Rotation.
    template <typename T2>
    void rotate(Point3D<T> p, const T2& ax, const T2& ay, const T2& az) {
        center_pos.translate_rotate(p, ax, ay, az);
        for (auto& face : faces) {
            face.rotate(p, ax, ay, az);
        }
    }
};

/// @brief Classe Scène (Conteneur).
template <typename T>
class Scene3D {
private:
    std::vector<Pave3D<T>> paves;
    std::vector<Sphere3D<T>> spheres;
    Point3D<T> camera_pos;
    Point3D<T> look_target;

public:
    Scene3D() : camera_pos(0, 0, -50), look_target(0, 0, 0) {}

    void add_pave(const Pave3D<T>& p) { paves.push_back(p); }
    void add_sphere(const Sphere3D<T>& s) { spheres.push_back(s); }

    const std::vector<Pave3D<T>>& get_paves() const { return paves; }
    std::vector<Pave3D<T>>& get_paves_mut() { return paves; }
    const std::vector<Sphere3D<T>>& get_spheres() const { return spheres; }
    std::vector<Sphere3D<T>>& get_spheres_mut() { return spheres; }
    Point3D<T>& get_camera() { return camera_pos; }

    void clear() { paves.clear(); spheres.clear(); }
};

/// @brief Point 2D écran.
class Point2d {
public:
    int x, y;
    float depth;
    Point2d() : x(0), y(0), depth(0.0f) {}
    Point2d(int _x, int _y) : x(_x), y(_y), depth(0.0f) {}
    Point2d(int _x, int _y, float _z) : x(_x), y(_y), depth(_z) {}
};

/// @brief Triangle 2D écran.
class Triangle2d {
public:
    Point2d p1, p2, p3;
    Triangle2d() {}
    Triangle2d(Point2d a, Point2d b, Point2d c) : p1(a), p2(b), p3(c) {}
};

#endif