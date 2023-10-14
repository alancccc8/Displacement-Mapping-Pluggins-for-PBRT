# Displacement-Mapping-Pluggins-for-PBRT
PBRT is a rendering engine widely used in the computer graphics area. This plugin adds a new 'displacement map' function based on PBRT-v3. Users can use the command console to easily add a displacement mapping to the Sphere shape, which realizes a more realistic effect.   
The following is an example using PBRT to render.
![image](https://github.com/alancccc8/Displacement-Mapping-Pluggins-for-PBRT/assets/93360333/39640991-0e65-499c-8e9e-e67b52cd3d44)
The displacemnet mapping not only changes the appearance of the model, but also changes the physic model, bring a more realistic effect. Shown as below.
![image](https://github.com/alancccc8/Displacement-Mapping-Pluggins-for-PBRT/assets/93360333/86316013-3030-4d78-876a-b481e3ce2cf6)
Intersection Test:
```jsx
Bounds3f ObjectBound() const;
    bool Intersect(const Ray &ray, Float *tHit, SurfaceInteraction *isect) const;
    bool IntersectP(const Ray &ray) const;
```
