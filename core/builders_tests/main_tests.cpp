#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdint>
#include <cassert>

#include "../builders/cbvh_core.h"
#include "../builders/cbvh_internal.h"

#include "bvhree_tests.h"
#include "cbvh_metrics.h"

#include "cmesh.h"


bool RunAllTestsForVSGFMesh(std::string bunny, cbvh::BVHPresets presets = cbvh::BVHPresets(), std::ofstream *pFout = nullptr)
{
  std::cout << "loading ... '" << bunny.c_str()  << "'" << std::endl;

  auto mesh = cmesh::LoadMeshFromVSGF(bunny.c_str());
  std::cout << "mesh " << bunny.c_str() << ":  verts/tris  = " << mesh.VerticesNum() << " / " << mesh.indices.size()/3 << std::endl;
  
  if(mesh.VerticesNum() == 0)
    return false;
  
  std::cout << "building ... " << std::endl;
  
  cbvh::BVHTree treeData;
  bool testsForBvh2 = true;
  if(presets.childrenNum == 2) 
  {
    cbvh::BVHTree notConverted, notConverted2;
    presets.desiredFormat = cbvh::FMT_BVH2Node32_Interval32_Dynamic;
    notConverted = cbvh::BuildBVH(mesh.vPos4f.data(), mesh.VerticesNum(), mesh.indices.data(), mesh.indices.size(), presets);
    
    presets.desiredFormat = cbvh::FMT_BVH2Node32_Interval32_Static;
    notConverted2 = cbvh::BuildBVH(mesh.vPos4f.data(), mesh.VerticesNum(), mesh.indices.data(), mesh.indices.size(), presets);
    
    std::cout << "testing (2-leaves-direct) ... " << std::endl;

    treeData     = cbvh::ConvertBVH2DynamicToBVH4Flat(notConverted);
    testsForBvh2 = cbvh_test::TestTrisIndideBoxes(notConverted, mesh.vPos4f.data());
    if(testsForBvh2)
      std::cout << "[TestBVH2TrisIndicesBoxes]:\t\tPASSED!" << std::endl; 
    else
      std::cout << "[TestBVH2TrisIndicesBoxes]:\t\tFAILED!" << std::endl;
    
    bool test1 = cbvh_test::TestLeafTrisSet(notConverted, mesh.indices.data(), mesh.indices.size(), mesh.vPos4f.data());
    if(test1)
      std::cout << "[TestBVH2LeafTrisSet_S]:\tPASSED!" << std::endl;
    else
      std::cout << "[TestBVH2LeafTrisSet_S]:\tFAILED!" << std::endl;
    
    bool test2 = cbvh_test::TestLeafTrisSet(notConverted2, mesh.indices.data(), mesh.indices.size(), mesh.vPos4f.data());
    if(test2)
      std::cout << "[TestBVH2LeafTrisSet_D]:\tPASSED!" << std::endl;
    else
      std::cout << "[TestBVH2LeafTrisSet_D]:\tFAILED!" << std::endl;
    
    testsForBvh2 = testsForBvh2 && test1 && test2;
  }
  else
    treeData = cbvh::BuildBVH(mesh.vPos4f.data(), mesh.VerticesNum(), mesh.indices.data(), mesh.indices.size(), presets);
  
  std::cout << "testing ... " << std::endl;
  bool allPassed = cbvh_test::RunAllTests(treeData, mesh.vPos4f.data(), mesh.VerticesNum(), mesh.indices.data(), mesh.indices.size()) && testsForBvh2;

  
  std::cout << "metrics ..." << std::endl;
  const float SAH = cbvh_metrics::MetricSAH   (treeData, mesh.vPos4f.data(), mesh.VerticesNum());
  std::cout << "SAH:\t" << SAH << std::endl;
  const float EPO = cbvh_metrics::MetricEPOTri(treeData, mesh.vPos4f.data(), mesh.VerticesNum());
  std::cout << "EPO:\t" << EPO << std::endl;

  auto rayMetrics = cbvh_test::ComputeRayTraversalMetrics(treeData, mesh.vPos4f.data(), mesh.vPos4f.size(), mesh.indices.data(), mesh.indices.size());
    
  std::cout << "LCV:\t"    << rayMetrics.varLeafesCount << std::endl;
  std::cout << "NC_avg:\t" << rayMetrics.avgNodesCount << std::endl;
  std::cout << "TC_avg:\t" << rayMetrics.avgPrimsCount << std::endl;
  std::cout << std::endl;

  if(pFout != nullptr)
  {
    auto& fout = *pFout;
    const size_t last_slash_idx = bunny.find_last_of("\\/");
    const size_t dot_idx        = bunny.find_last_of(".");

    if(last_slash_idx == std::string::npos || dot_idx == std::string::npos)
      return allPassed;

    const std::string sceneName = bunny.substr(last_slash_idx+1, dot_idx-last_slash_idx-1);
    fout << sceneName.c_str() << ";";
    if(presets.desiredFormat == cbvh::FMT_BVH2Node32_Interval32_Dynamic)
      fout << "lbvh;";
    else if(presets.btype == cbvh::BVH_CONSTRUCT_QUALITY)
      fout << "esc2;";
    else
      fout << "sweep;";

    fout << SAH*10.0f << ";";
    fout << EPO << ";";
    fout << rayMetrics.avgNodesCount << ";";
    fout << rayMetrics.avgPrimsCount << ";";
    fout << rayMetrics.varLeafesCount << ";";
    fout << std::endl;
  }

  return allPassed;
}

void RunAllTestsForBoxesFromMesh(std::string bunny, cbvh::BVHPresets presets = cbvh::BVHPresets())
{
  std::cout << "loading ... " << std::endl;

  auto mesh = cmesh::LoadMeshFromVSGF(bunny.c_str());
  std::cout << "mesh " << bunny.c_str() << ": inds.size()  = " << mesh.VerticesNum() << " / " << mesh.indices.size()/3 << std::endl;
  
  if(mesh.VerticesNum() == 0)
    return;
  
  auto boxes = cbvh_test::TriangleBoxesFromMesh(mesh.vPos4f.data(), mesh.VerticesNum(), mesh.indices.data(), mesh.indices.size());

  std::cout << "building ... " << std::endl;
  
  cbvh::BVHTree treeData;
  if(presets.childrenNum == 2) 
  {
    presets.desiredFormat = cbvh::FMT_BVH2Node32_Interval32_Dynamic;
    treeData = cbvh::ConvertBVH2DynamicToBVH4Flat(cbvh::BuildBVH((const cbvh::BVHNode*)boxes.data(), boxes.size(), presets));
  }
  else
    treeData = cbvh::BuildBVH((const cbvh::BVHNode*)boxes.data(), boxes.size(), presets);

  // now just get vertices indices from boxes
  //
  {
    std::vector<uint32_t> triIndices(mesh.indices.size());
    assert(boxes.size() == treeData.indicesReordered.size());

    for(size_t i=0;i<boxes.size();i++)
    {
      const uint32_t boxIndex = treeData.indicesReordered[i];
      const auto& box         = boxes[boxIndex];
      assert(box.getCount() == 1);
      uint32_t triId = box.getStart();
      triIndices[3*i+0] = mesh.indices[triId*3+0];
      triIndices[3*i+1] = mesh.indices[triId*3+1];
      triIndices[3*i+2] = mesh.indices[triId*3+2];
    }

    treeData.indicesReordered = triIndices;
  }

  std::cout << "testing ... " << std::endl;
  cbvh_test::RunAllTests(treeData, mesh.vPos4f.data(), mesh.VerticesNum(), mesh.indices.data(), mesh.indices.size());

  std::cout << "metrics ..." << std::endl;
  std::cout << "SAH: " << cbvh_metrics::MetricSAH(treeData, mesh.vPos4f.data(), mesh.VerticesNum()) << std::endl;
  std::cout << "EPO: " << cbvh_metrics::MetricEPOTri(treeData, mesh.vPos4f.data(), mesh.VerticesNum()) << std::endl << std::endl;
}

int main(int argc, const char** argv)
{
  assert(sizeof(cbvh::BVHNode)  == sizeof(LiteMath::float4)*2);
  assert(sizeof(cbvh::Interval) == sizeof(LiteMath::uint2));

  #ifdef WIN32
  std::string meshFolder = "../../../../cbvh_stf/resources/meshes/";
  #else
  std::string meshFolder = "../../../scenes/meshes/";
  #endif
  
  //cbvh::BVHPresets presets;
  //presets.childrenNum   = 2;
  //presets.btype         = cbvh::BVH_CONSTRUCT_FAST;
  //presets.desiredFormat = cbvh::FMT_BVH2Node32_Interval32_Dynamic;
  //RunAllTestsForVSGFMesh(meshFolder + "cone_low.vsgf", presets);

  //cbvh::BVHPresets presets;
  //presets.childrenNum   = 2;
  //presets.btype         = cbvh::BVH_CONSTRUCT_MEDIUM;
  //presets.desiredFormat = cbvh::FMT_BVH2Node32_Interval32_Static;
  //RunAllTestsForVSGFMesh(meshFolder + "sponza.vsgf", presets);

  //cbvh::BVHPresets presets;
  //presets.childrenNum   = 4;
  //presets.btype         = cbvh::BVH_CONSTRUCT_MEDIUM;
  //presets.desiredFormat = cbvh::FMT_BVH4Node32_Interval32_Static;
  //RunAllTestsForVSGFMesh(meshFolder + "sponza.vsgf", presets);

  
  cbvh::BVHPresets allTypesOfPresets[3] = {};
  {
    allTypesOfPresets[0].childrenNum = 2;
    allTypesOfPresets[0].btype       = cbvh::BVH_CONSTRUCT_FAST;
  
    allTypesOfPresets[1].childrenNum = 4;
    allTypesOfPresets[1].btype       = cbvh::BVH_CONSTRUCT_MEDIUM;
  
    allTypesOfPresets[2].childrenNum = 4;
    allTypesOfPresets[2].btype       = cbvh::BVH_CONSTRUCT_QUALITY;
  }
  
  std::ofstream fout("z_metrics.csv", std::ios::app);
  bool passed = true;
  for(int i=0;i<3;i++)
  {
    fout << "preset_id = " << i << ";" << std::endl;
    fout << "scene; builder; SAH; EPO; NC; TC; LCV;" << std::endl;

    std::cout << "begin tests for preset #(" << i << ")" << std::endl;
    passed = passed && RunAllTestsForVSGFMesh(meshFolder + "cone_low.vsgf",     allTypesOfPresets[i], &fout);
    passed = passed && RunAllTestsForVSGFMesh(meshFolder + "cylinder_low.vsgf", allTypesOfPresets[i], &fout);
    passed = passed && RunAllTestsForVSGFMesh(meshFolder + "helix_mid.vsgf",    allTypesOfPresets[i], &fout);
    passed = passed && RunAllTestsForVSGFMesh(meshFolder + "teapot.vsgf",       allTypesOfPresets[i], &fout);
    passed = passed && RunAllTestsForVSGFMesh(meshFolder + "bunny0.vsgf",       allTypesOfPresets[i], &fout);
    passed = passed && RunAllTestsForVSGFMesh(meshFolder + "sponza.vsgf",       allTypesOfPresets[i], &fout);
    passed = passed && RunAllTestsForVSGFMesh(meshFolder + "lucy.vsgf");
    std::cout << "end  tests for preset #(" << i << ")" << std::endl << std::endl;
    fout << std::endl << std::endl;
  }
  if(passed)
    std::cout << "ALL TESTS --> PASSED!" << std::endl;
  else
    std::cout << "SOME TEST --> FAILED!" << std::endl;

  return 0;
}

