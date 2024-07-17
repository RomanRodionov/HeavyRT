# cbvh_raytrace

This is the core repository for cooperation RT project. Originally commits are based on several open source repositories. They are:

* https://github.com/msu-graphics-group/LiteMath
* https://github.com/msu-graphics-group/scenes
* https://github.com/Ray-Tracing-Systems/HydraCore

## Repository content

The repository consists of:
* 'core' directory -- implements Ray Tracing on CPU
  * 'core/builders' -- CPU and GPU based builders
    * 'core/builders/cbvh.h' -- API for builders
    * 'core/builders/lbvh.h' -- API for lbvh builder
    * 'core/builders_tests' -- some internals tets for CPU builders
* 'examples' which uses RT implemantation for testing and other purposes;
  * 'examples/01_eye_rays' -- implementation of ray tracing for primary rays
  * 'examples/02_rtao'     -- implementation of ray tracing for Ray Traced Ambient Occlusion

## Add your own RT algorithm

1. Add your own implementation of 'ISceneObject' interface ('CrossRT.h') in some file (lets call it "MyRT.cpp");
2. Create factory function 'ISceneObject* CreateMyRT(const char* a_implName)' that should return a pointer to your implementation;
3. Add declaration of your function 'CreateMyRT' inside FactoryRT.cpp;
4. Add check for you class name inside 'CreateSceneRT' function in FactoryRT.cpp;
5. You are done!

## Add your own GI/RT sample

1. Add your own implementation of 'IRenderer' interface ('IRenderer.h') in some file (lets call it "MyRenderer.cpp");
2. Create factory function 'IRenderer* CreateMyRenderer(const char* a_implName)' that should return a pointer to your implementation;
3. Add declaration of your function 'CreateMyRenderer' inside RenderFactory.cpp;
4. Add check for you class name inside 'CreateRender' function in RenderFactory.cpp;
5. You are done!

### Scenes and shaders
* TBD

# Acknowlegments 

We sincerely thank Huawei for contribution to this Open Source project. Guys, you are the best!