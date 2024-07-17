#include "CrossRT.h"
#include <string>
#include <iostream>

ISceneObject* CreateEmbreeRT(const char* a_implName);
ISceneObject* MakeBruteForceRT(const char* a_implName);

ISceneObject* MakeBVH2CommonRT(const char* a_implName, const char* a_buildName);
ISceneObject *MakeBVH2CommonRTStacklessLBVH(const char *a_implName, const char* a_buildName);
ISceneObject *MakeBVH2CommonLoftRT(const char *a_implName, const char* a_buildName);
ISceneObject *MakeBVH2CommonLoft16RT(const char *a_implName, const char* a_buildName); 
ISceneObject *MakeBVH2_LRFT(const char *a_implName, const char* a_buildName);

ISceneObject* MakeBVH2FatRT(const char* a_implName, const char* a_buildName, const char* a_layoutName);
ISceneObject* MakeBVH2FatRTCompact(const char* a_implName, const char* a_buildName, const char* a_layoutName);
ISceneObject* MakeBVH2Fat16RT(const char* a_implName, const char* a_buildName, const char* a_layoutName);
ISceneObject* MakeBVH4CommonRT(const char* a_implName, const char* a_buildName);
//ISceneObject* MakeBVH4CommonLoftRT(const char* a_implName, const char* a_buildName);
ISceneObject* MakeBVH4HalfRT(const char* a_implName, const char* a_buildName);

ISceneObject *MakeBVH2Stackless(const char *a_implName, const char* a_buildName);
ISceneObject* MakeNanoRT(const char* a_implName);
ISceneObject* MakeNanoRTExt(const char* a_implName);
ISceneObject* MakeShtRT(const char* a_implName, const char* a_buildName);
ISceneObject* MakeShtRT64(const char* a_implName, const char* a_buildName);
ISceneObject* MakeSht4NodesRT(const char* a_implName, const char* a_buildName);
ISceneObject* MakeSht4NodesHalfRT(const char* a_implName, const char* a_buildName);

ISceneObject* CreateSceneRT(const char* a_implName, const char* a_buildName, const char* a_layoutName)
{
  const std::string className(a_implName); //
  if (className.find("NanoRTExt") != std::string::npos)
    return MakeNanoRTExt(a_implName);
  else if (className.find("NanoRT") != std::string::npos)
    return MakeNanoRT(a_implName);
  else if (className.find("BruteForce") != std::string::npos)
    return MakeBruteForceRT(a_implName);
  else if (className.find("BVH2CommonRTStacklessLBVH") != std::string::npos || className.find("BVH2CommonSTAS") != std::string::npos || className.find("BVH2_STAS") != std::string::npos)
    return MakeBVH2CommonRTStacklessLBVH(a_implName, a_buildName);
  else if (className.find("BVH2CommonLoft16") != std::string::npos || className.find("BVH2CommonLOFT16") != std::string::npos || className.find("BVH2_LOFT_16") != std::string::npos)
    return MakeBVH2CommonLoft16RT(a_implName, a_buildName);
  else if (className.find("BVH2CommonLoft") != std::string::npos || className.find("BVH2CommonLOFT") != std::string::npos || className.find("BVH2_LOFT") != std::string::npos)
    return MakeBVH2CommonLoftRT(a_implName, a_buildName);
   else if (className.find("BVH2_LRFT") != std::string::npos)
    return MakeBVH2_LRFT(a_implName, a_buildName);
  else if (className.find("BVH2Common") != std::string::npos || className.find("BVH2Default") != std::string::npos)
    return MakeBVH2CommonRT(a_implName, a_buildName);
  else if (className.find("BVH2FatCompact") != std::string::npos)
    return MakeBVH2FatRTCompact(a_implName, a_buildName, a_layoutName);
  else if (className.find("BVH2FatHalf") != std::string::npos || className.find("BVH2Fat16") != std::string::npos)
    return MakeBVH2Fat16RT(a_implName, a_buildName, a_layoutName);
  else if (className.find("BVH2Fat") != std::string::npos)
    return MakeBVH2FatRT(a_implName, a_buildName, a_layoutName);
  //else if(className.find("BVH4CommonLoft") != std::string::npos || className.find("BVH4Loft") != std::string::npos || className.find("BVH4_LOFT") != std::string::npos)
  //  return MakeBVH4CommonLoftRT(a_implName, a_buildName);
  else if(className.find("BVH4Half") != std::string::npos || className.find("BVH416") != std::string::npos || className.find("BVH4_16") != std::string::npos)
    return MakeBVH4HalfRT(a_implName, a_buildName);
  else if(className.find("BVH4Common") != std::string::npos || className.find("BVH4Default") != std::string::npos || className == "BVH4" || className == "BVH4_32")
    return MakeBVH4CommonRT(a_implName, a_buildName);
   else if(className.find("BVH2Stackless") != std::string::npos) 
    return MakeBVH2Stackless(a_implName, a_buildName);
  else if(className.find("Embree") != std::string::npos || className.find("embree") != std::string::npos)
    return CreateEmbreeRT(a_implName);
  else if (className.find("ShtRT64") != std::string::npos || className == "BVH2_SHRT64")
    return MakeShtRT64(a_implName,a_buildName);
  else if (className.find("ShtRT") != std::string::npos || className == "BVH2_SHRT32")
    return MakeShtRT(a_implName,a_buildName);
  else if (className.find("Sht4NodesHalfRT") != std::string::npos || className == "BVH4_SHRT_HALF")
    return MakeSht4NodesHalfRT(a_implName, a_buildName);
  else if (className.find("Sht4NodesRT") != std::string::npos || className == "BVH4_SHRT")
    return MakeSht4NodesRT(a_implName, a_buildName);
  else
  {
    std::cout << "[ERROR!] CreateSceneRT: unknown implementation!" << className.c_str() << std::endl;
    std::cout << "Use Embree" << std::endl;
    return CreateEmbreeRT(a_implName);
  }
}

void DeleteSceneRT(ISceneObject* a_impl) { delete a_impl; }

MetricStats ISceneObject::GetStats() 
{
  MetricStats stats = {};
  const double rays = std::max(double(m_stats.raysNumber), 1.0);
  stats.avgNC  = float(double(m_stats.NC)/rays);
  stats.avgLC  = float(double(m_stats.LC)/rays);
  stats.avgTC  = float(double(m_stats.TC)/rays);
  for (int i = 0; i < TREELET_ARR_SIZE; i++) {
    stats.avgLJC[i] = float(double(m_stats.LJC[i])/rays);
    stats.avgCMC[i] = float(double(m_stats.CMC[i])/rays);
    stats.avgWSS[i] = uint32_t(m_stats.WSS[i].size());
  }
  stats.avgBLB = float(double(m_stats.BLB)/rays);
  stats.avgSOC = float(double(m_stats.SOC)/rays);
  stats.avgSBL = float(double(m_stats.SBL)/rays);
  stats.bvhTotalSize  = m_stats.bvhTotalSize;
  stats.geomTotalSize = m_stats.geomTotalSize;
  return stats;
}

void ISceneObject::ResetStats()
{
  m_stats.clear();
}