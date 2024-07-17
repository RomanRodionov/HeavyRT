#include <fstream>
#include <iostream>
#include "gltf_loader.h"
#include "stb_image.h"


namespace gltf_loader
{
  bool ends_with (const std::string &fullString, const std::string &ending)
  {
    if (fullString.length() >= ending.length())
    {
      return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
    }
    else
    {
      return false;
    }
  }

  std::pair<LiteMath::float4x4, LiteMath::float4x4> makeDefaultCamera(float aspect)
  {
    auto worldView = lookAt(float3(0, 1, 1), float3(0, 0, -1), float3(0, 1, 0));

    auto proj = perspectiveMatrix(45, aspect, 0.01f, 1000.0f);

    return {inverse4x4(worldView), inverse4x4(proj)};
  }

  float vfov_from_hfov(float hfov, uint32_t w, uint32_t h)
  {
    return (2.0f * atanf((0.5f * float(h)) /
                         (0.5f * float(w) / tanf(DEG_TO_RAD * hfov / 2.0f)))) / DEG_TO_RAD;
  }

  float hfov_from_vfov(float hfov, uint32_t w, uint32_t h)
  {
    return (2.0f * atanf((0.5f * float(w)) /
                         (0.5f * float(h) / tanf(DEG_TO_RAD * hfov / 2.0f)))) / DEG_TO_RAD;
  }

  std::pair<LiteMath::float4x4, LiteMath::float4x4> makeCameraFromSceneBBox(uint32_t render_width, uint32_t render_height, Box4f box)
  {
    float aspect = float(render_width) / float(render_height);
    float3 sceneCenter = 0.5f * (to_float3(box.boxMin) + to_float3(box.boxMax));
    float fovY = 45;
    float fovX = hfov_from_vfov(fovY, render_width, render_height);

    float dist1 = std::fabs(box.boxMax.y - box.boxMin.y) * 0.5f / std::tan(fovY * 0.5f * DEG_TO_RAD);
    float dist2 = std::fabs(box.boxMax.x - box.boxMin.x) * 0.5f / std::tan(fovX * 0.5f * DEG_TO_RAD);
    float dist = std::max(dist1, dist2) * 1.25f;

    float3 eye = {sceneCenter.x, sceneCenter.y, std::max(box.boxMax.z, box.boxMin.z) + dist };

    auto worldView = lookAt(eye, sceneCenter, float3(0, 1, 0));
    auto proj = perspectiveMatrix(fovY, aspect, 0.01f, 1000.0f);
    
    std::cout << "[GLTF]: make camera: " << std::endl;
    std::cout << "[GLTF]: camPos    = (" << eye.x << "," << eye.y << "," << eye.z << ")" << std::endl;
    std::cout << "[GLTF]: camLookAT = (" << sceneCenter.x << "," << sceneCenter.y << "," << sceneCenter.z << ")" << std::endl;

    return {inverse4x4(worldView), inverse4x4(proj)};
  }

  LiteMath::float4x4 rotMatrixFromQuaternion(const LiteMath::float4 &q)
  {
    LiteMath::float4x4 rot{};

    rot.m_col[0][0] = 1 - 2 * (q[1] * q[1] + q[2] * q[2]);
    rot.m_col[0][1] = 2 * (q[0] * q[1] + q[2] * q[3]);
    rot.m_col[0][2] = 2 * (q[0] * q[2] - q[1] * q[3]);
    rot.m_col[0][3] = 0.0f;

    rot.m_col[1][0] = 2 * (q[0] * q[1] - q[2] * q[3]);
    rot.m_col[1][1] = 1 - 2 * (q[0] * q[0] + q[2] * q[2]);
    rot.m_col[1][2] = 2 * (q[1] * q[2] + q[0] * q[3]);
    rot.m_col[1][3] = 0.0f;

    rot.m_col[2][0] = 2 * (q[0] * q[2] + q[1] * q[3]);
    rot.m_col[2][1] = 2 * (q[1] * q[2] - q[0] * q[3]);
    rot.m_col[2][2] = 1 - 2 * (q[0] * q[0] + q[1] * q[1]);
    rot.m_col[2][3] = 0.0f;

    rot.m_col[3][0] = 0.0f;
    rot.m_col[3][1] = 0.0f;
    rot.m_col[3][2] = 0.0f;
    rot.m_col[3][3] = 1.0f;

    return rot;
  }

  LiteMath::float4x4 transformMatrixFromGLTFNode(const tinygltf::Node &node)
  {
    LiteMath::float4x4 nodeMatrix;

    if(node.matrix.size() == 16)
    {
      nodeMatrix.set_col(0, float4(node.matrix[0], node.matrix[1], node.matrix[2], node.matrix[3]));
      nodeMatrix.set_col(1, float4(node.matrix[4], node.matrix[5], node.matrix[6], node.matrix[7]));
      nodeMatrix.set_col(2, float4(node.matrix[8], node.matrix[9], node.matrix[10], node.matrix[11]));
      nodeMatrix.set_col(3, float4(node.matrix[12], node.matrix[13], node.matrix[14], node.matrix[15]));
    }
    else
    {
      if(node.scale.size() == 3)
      {
        LiteMath::float3 s = LiteMath::float3(node.scale[0], node.scale[1], node.scale[2]);
        nodeMatrix         = LiteMath::scale4x4(s) * nodeMatrix;
      }
      if(node.rotation.size() == 4)
      {
        LiteMath::float4 rot = LiteMath::float4(node.rotation[0], node.rotation[1], node.rotation[2], node.rotation[3]);
        nodeMatrix           = rotMatrixFromQuaternion(rot) * nodeMatrix;
      }
      if(node.translation.size() == 3)
      {
        LiteMath::float3 t = LiteMath::float3(node.translation[0], node.translation[1], node.translation[2]);
        nodeMatrix         = LiteMath::translate4x4(t) * nodeMatrix;
      }

    }

    return nodeMatrix;
  }


  MaterialDataGLTF materialDataFromGLTF(const tinygltf::Material &gltfMat, const std::vector<tinygltf::Texture> &gltfTextures)
  {
    MaterialDataGLTF mat = {};
    auto& baseColor = gltfMat.pbrMetallicRoughness.baseColorFactor;
    mat.baseColor   = LiteMath::float4(baseColor[0], baseColor[1], baseColor[2], baseColor[3]);
    mat.alphaCutoff = gltfMat.alphaCutoff;

    if(gltfMat.alphaMode == "OPAQUE")
      mat.alphaMode = 0;
    else if(gltfMat.alphaMode == "MASK")
      mat.alphaMode = 1;
    else if(gltfMat.alphaMode == "BLEND")
      mat.alphaMode = 2;

    mat.metallic    = gltfMat.pbrMetallicRoughness.metallicFactor;
    mat.roughness   = gltfMat.pbrMetallicRoughness.roughnessFactor;
    auto& emission = gltfMat.emissiveFactor;
    mat.emissionColor  = LiteMath::float3(emission[0], emission[1], emission[2]);

    mat.baseColorTexId = gltfMat.pbrMetallicRoughness.baseColorTexture.index;
    mat.emissionTexId  = gltfMat.emissiveTexture.index;
    mat.normalTexId    = gltfMat.normalTexture.index;
    mat.occlusionTexId = gltfMat.occlusionTexture.index;
    mat.metallicRoughnessTexId = gltfMat.pbrMetallicRoughness.metallicRoughnessTexture.index;

    return mat;
  }

  void getNumVerticesAndIndicesFromGLTFMesh(const tinygltf::Model &a_model, const tinygltf::Mesh &a_mesh,
                                            uint32_t& numVertices, uint32_t& numIndices)
  {
    auto numPrimitives   = a_mesh.primitives.size();
    for(size_t j = 0; j < numPrimitives; ++j)
    {
      const tinygltf::Primitive &glTFPrimitive = a_mesh.primitives[j];
      if(glTFPrimitive.attributes.find("POSITION") != glTFPrimitive.attributes.end())
      {
        numVertices += a_model.accessors[glTFPrimitive.attributes.find("POSITION")->second].count;
      }

      numIndices += static_cast<uint32_t>(a_model.accessors[glTFPrimitive.indices].count);
    }
  }

  cmesh::SimpleMesh2 simpleMeshFromGLTFMesh(const tinygltf::Model &a_model, const tinygltf::Mesh &a_mesh)
  {
    uint32_t numVertices = 0;
    uint32_t numIndices  = 0;
    auto numPrimitives   = a_mesh.primitives.size();

    getNumVerticesAndIndicesFromGLTFMesh(a_model, a_mesh, numVertices, numIndices);

    auto simpleMesh = cmesh::SimpleMesh2(numVertices, numIndices);

    uint32_t firstIndex  = 0;
    uint32_t vertexStart = 0;
    for(size_t j = 0; j < numPrimitives; ++j)
    {
      const tinygltf::Primitive &glTFPrimitive = a_mesh.primitives[j];

      // Vertices
      size_t vertexCount = 0;
      {
        const float *positionBuffer  = nullptr;
        const float *normalsBuffer   = nullptr;
        const float *texCoordsBuffer = nullptr;
        const float *tangentsBuffer  = nullptr;

        if(glTFPrimitive.attributes.find("POSITION") != glTFPrimitive.attributes.end())
        {
          const tinygltf::Accessor &accessor = a_model.accessors[glTFPrimitive.attributes.find("POSITION")->second];
          const tinygltf::BufferView &view   = a_model.bufferViews[accessor.bufferView];
          positionBuffer = reinterpret_cast<const float *>(&(a_model.buffers[view.buffer].data[accessor.byteOffset + view.byteOffset]));
          vertexCount                        = accessor.count;

          if(accessor.minValues.size() == 3 && accessor.maxValues.size() == 3)
          {
            simpleMesh.bbox.boxMin = float3(accessor.minValues[0], accessor.minValues[1], accessor.minValues[2]);
            simpleMesh.bbox.boxMax = float3(accessor.maxValues[0], accessor.maxValues[1], accessor.maxValues[2]);
          }
        }

        if(glTFPrimitive.attributes.find("NORMAL") != glTFPrimitive.attributes.end())
        {
          const tinygltf::Accessor &accessor = a_model.accessors[glTFPrimitive.attributes.find("NORMAL")->second];
          const tinygltf::BufferView &view   = a_model.bufferViews[accessor.bufferView];
          normalsBuffer = reinterpret_cast<const float *>(&(a_model.buffers[view.buffer].data[accessor.byteOffset + view.byteOffset]));
        }

        if(glTFPrimitive.attributes.find("TEXCOORD_0") != glTFPrimitive.attributes.end())
        {
          const tinygltf::Accessor &accessor = a_model.accessors[glTFPrimitive.attributes.find("TEXCOORD_0")->second];
          const tinygltf::BufferView &view   = a_model.bufferViews[accessor.bufferView];
          texCoordsBuffer = reinterpret_cast<const float *>(&(a_model.buffers[view.buffer].data[accessor.byteOffset + view.byteOffset]));
        }

        if(glTFPrimitive.attributes.find("TANGENT") != glTFPrimitive.attributes.end())
        {
          const tinygltf::Accessor &accessor = a_model.accessors[glTFPrimitive.attributes.find("TANGENT")->second];
          const tinygltf::BufferView &view   = a_model.bufferViews[accessor.bufferView];
          tangentsBuffer = reinterpret_cast<const float *>(&(a_model.buffers[view.buffer].data[accessor.byteOffset + view.byteOffset]));
        }

        for(size_t v = 0; v < vertexCount; v++)
        {
          simpleMesh.vPos4f[vertexStart + v][0] = positionBuffer[v * 3 + 0];
          simpleMesh.vPos4f[vertexStart + v][1] = positionBuffer[v * 3 + 1];
          simpleMesh.vPos4f[vertexStart + v][2] = positionBuffer[v * 3 + 2];
          simpleMesh.vPos4f[vertexStart + v][3] = 1.0f;

          simpleMesh.vNorm4f[vertexStart + v][0] = normalsBuffer ? normalsBuffer[v * 3 + 0] : 0.0f;
          simpleMesh.vNorm4f[vertexStart + v][1] = normalsBuffer ? normalsBuffer[v * 3 + 1] : 0.0f;
          simpleMesh.vNorm4f[vertexStart + v][2] = normalsBuffer ? normalsBuffer[v * 3 + 2] : 0.0f;
          simpleMesh.vNorm4f[vertexStart + v][3] = normalsBuffer ? 1.0f : 0.0f;

          simpleMesh.vTexCoord2f[vertexStart + v][0] = texCoordsBuffer ? texCoordsBuffer[v * 2 + 0] : 0.0f;
          simpleMesh.vTexCoord2f[vertexStart + v][1] = texCoordsBuffer ? texCoordsBuffer[v * 2 + 1] : 0.0f;

          simpleMesh.vTang4f[vertexStart + v][0] = tangentsBuffer ? tangentsBuffer[v * 4 + 0] : 0.0f;
          simpleMesh.vTang4f[vertexStart + v][1] = tangentsBuffer ? tangentsBuffer[v * 4 + 1] : 0.0f;
          simpleMesh.vTang4f[vertexStart + v][2] = tangentsBuffer ? tangentsBuffer[v * 4 + 2] : 0.0f;
          simpleMesh.vTang4f[vertexStart + v][3] = tangentsBuffer ? tangentsBuffer[v * 4 + 3] : 0.0f;
        }
      }

      // Indices
      {
        const tinygltf::Accessor &accessor     = a_model.accessors[glTFPrimitive.indices];
        const tinygltf::BufferView &bufferView = a_model.bufferViews[accessor.bufferView];
        const tinygltf::Buffer &buffer         = a_model.buffers[bufferView.buffer];

        auto indexCount = static_cast<uint32_t>(accessor.count);

        std::fill(simpleMesh.matIndices.begin() + firstIndex / 3,
                  simpleMesh.matIndices.begin() + (firstIndex + indexCount) / 3, glTFPrimitive.material);

        switch(accessor.componentType)
        {
          case TINYGLTF_PARAMETER_TYPE_UNSIGNED_INT:
          {
            auto *buf = new uint32_t[accessor.count];
            memcpy(buf, &buffer.data[accessor.byteOffset + bufferView.byteOffset], accessor.count * sizeof(uint32_t));
            for(size_t index = 0; index < accessor.count; index++)
            {
              simpleMesh.indices[firstIndex + index] = buf[index] + vertexStart;
            }
            break;
          }
          case TINYGLTF_PARAMETER_TYPE_UNSIGNED_SHORT:
          {
            auto *buf = new uint16_t[accessor.count];
            memcpy(buf, &buffer.data[accessor.byteOffset + bufferView.byteOffset], accessor.count * sizeof(uint16_t));
            for(size_t index = 0; index < accessor.count; index++)
            {
              simpleMesh.indices[firstIndex + index] = buf[index] + vertexStart;
            }
            break;
          }
          case TINYGLTF_PARAMETER_TYPE_UNSIGNED_BYTE:
          {
            auto *buf = new uint8_t[accessor.count];
            memcpy(buf, &buffer.data[accessor.byteOffset + bufferView.byteOffset], accessor.count * sizeof(uint8_t));
            for(size_t index = 0; index < accessor.count; index++)
            {
              simpleMesh.indices[firstIndex + index] = buf[index] + vertexStart;
            }
            break;
          }
          default:
//             "[LoadSceneGLTF]: Unsupported index component type");
            return { };
        }

        firstIndex  += indexCount;
        vertexStart += vertexCount;
      }
    }

    return simpleMesh;
  }

  LiteImage::Sampler samplerFromGLTF(const tinygltf::Sampler &gltfSampler)
  {
    LiteImage::Sampler sampler {};
    if (gltfSampler.minFilter == TINYGLTF_TEXTURE_FILTER_NEAREST)
      sampler.filter = LiteImage::Sampler::Filter::NEAREST;
    else if (gltfSampler.minFilter == TINYGLTF_TEXTURE_FILTER_LINEAR)
      sampler.filter = LiteImage::Sampler::Filter::LINEAR;

    if(gltfSampler.wrapS == TINYGLTF_TEXTURE_WRAP_REPEAT)
      sampler.addressU = LiteImage::Sampler::AddressMode::WRAP;
    else if(gltfSampler.wrapS == TINYGLTF_TEXTURE_WRAP_CLAMP_TO_EDGE)
      sampler.addressU = LiteImage::Sampler::AddressMode::CLAMP;
    else if(gltfSampler.wrapS == TINYGLTF_TEXTURE_WRAP_MIRRORED_REPEAT)
      sampler.addressU = LiteImage::Sampler::AddressMode::MIRROR;

    if(gltfSampler.wrapT == TINYGLTF_TEXTURE_WRAP_REPEAT)
      sampler.addressV = LiteImage::Sampler::AddressMode::WRAP;
    else if(gltfSampler.wrapT == TINYGLTF_TEXTURE_WRAP_CLAMP_TO_EDGE)
      sampler.addressV = LiteImage::Sampler::AddressMode::CLAMP;
    else if(gltfSampler.wrapT == TINYGLTF_TEXTURE_WRAP_MIRRORED_REPEAT)
      sampler.addressV = LiteImage::Sampler::AddressMode::MIRROR;

    return sampler;
  }


  #if defined(__ANDROID__)
  std::shared_ptr<LiteImage::ICombinedImageSampler> LoadTextureAndMakeCombined(AAssetManager* assetManager, const tinygltf::Image& a_image, const LiteImage::Sampler& a_sampler, const std::string& a_sceneDir)
  #else
  std::shared_ptr<LiteImage::ICombinedImageSampler> LoadTextureAndMakeCombined(const tinygltf::Image& a_image, const LiteImage::Sampler& a_sampler, const std::string& a_sceneDir)
  #endif
  {
    std::shared_ptr<LiteImage::ICombinedImageSampler> pResult = nullptr;

    auto texturePath = a_sceneDir + a_image.uri;
    int bpp = a_image.component * (a_image.bits / 8);

    unsigned char* loadFromMem = nullptr;
    size_t memSize = 0;

    #if defined(__ANDROID__)
    AAsset* asset = AAssetManager_open(assetManager, texturePath.c_str(), AASSET_MODE_STREAMING);
    if(asset == nullptr)
      return nullptr;
    size_t size = AAsset_getLength(asset);
    assert(size > 0);
    std::vector<unsigned char> dataVector(size);
    AAsset_read(asset, dataVector.data(), size);
    loadFromMem = dataVector.data();
    memSize     = size;
    #endif

    if(bpp == 16)
    {
      int w, h, channels, req_channels;
      req_channels = 4;
      float* pixels = nullptr;

      if(loadFromMem)
        pixels = (float*)stbi_load_from_memory(loadFromMem, memSize, &w, &h, &channels, req_channels);
      else
        pixels = stbi_loadf(texturePath.c_str(), &w, &h, &channels, req_channels);

      std::vector<float> data(w * h * req_channels);
      memcpy(data.data(), pixels, data.size());

      stbi_image_free(pixels);

      auto pTexture = std::make_shared< LiteImage::Image2D<float4> >(w, h, (const float4*)data.data());
      pResult = MakeCombinedTexture2D(pTexture, a_sampler);
    }
    else
    {
      int w, h, channels, req_channels;
      req_channels = 4;
      unsigned char *pixels = nullptr;

      if(loadFromMem)
        pixels = stbi_load_from_memory(loadFromMem, memSize, &w, &h, &channels, req_channels);
      else
        pixels = stbi_load(texturePath.c_str(), &w, &h, &channels, req_channels);

      std::vector<uint32_t> data;
      if(pixels != nullptr)
      {
        data.resize(w * h);
        #pragma omp parallel for shared(h, w, data, pixels) default(none)
        for(int y = 0; y < h; ++y)
        {
          for(int x = 0; x < w; ++x)
          {

            data[x + y * w] =  pixels[(x + y * w) * 4 + 0] | (pixels[(x + y * w) * 4 + 1] << 8) |
                               (pixels[(x + y * w) * 4 + 2] << 16) | (pixels[(x + y * w) * 4 + 3] << 24);
          }
        }
        //memcpy(data.data(), pixels, data.size());
        stbi_image_free(pixels);
      }
      else
      {
        data.resize(1);
        data[0] = 0xFFFFFFFF;
        w = 1;
        h = 1;
      }
      auto pTexture = std::make_shared< LiteImage::Image2D<uint32_t> >(w, h, data.data());
      pTexture->setSRGB(true);
      pResult = MakeCombinedTexture2D(pTexture, a_sampler);
    }

    return pResult;
  }
}