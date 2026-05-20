// DISCRETE EXTERIOR CALCULUS

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"


#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"


#ifdef _WIN32
#include <Windows.h>
#include <commdlg.h>
#endif


#include "colormap.h"
#include "solvers.h"
#include <algorithm>
#include <map>
#include "spectral-conformal-parameterization.h"
// ** taken from Polyscope gl_engine.h -- not sure how to avoid using OpenGL for texture loading
#include "stb_image.h"
#ifdef __APPLE__
#define GLFW_INCLUDE_GLCOREARB
#include "GLFW/glfw3.h"
#else
#include "glad/glad.h"
// glad must come first
#include "GLFW/glfw3.h"
#endif

using namespace geometrycentral;
using namespace geometrycentral::surface;

std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;

polyscope::SurfaceMesh* psMesh;
std::string currentFilePath = "/Users/liguoxiong/Desktop/GeometryLab/Discrete-Geometry-and-Topology/input/costa.obj";
bool showFileDialog = false;
char filePathBuffer[256];

// 展平算法相关变量
enum class FlattenMethod { Tutte, LSCM, SCP,ARAP};
FlattenMethod currentFlattenMethod = FlattenMethod::Tutte;
polyscope::SurfaceParameterizationQuantity* checkerboard;
Vector3 CoM;                        // original center of mass, for re-centering purposes
VertexData<Vector3> ORIGINAL;       // original vertex positions
VertexData<Vector3> SCP_MESH;       // vertex positions of SCP
VertexData<Vector2> SCP_FLATTENING; // SCP solution
bool DISPLAY_FLAT = false;

// Windows文件选择对话框函数
#ifdef _WIN32
std::string openFileBrowser() {
    OPENFILENAME ofn;
    char fileName[256] = "";
    char initialDir[MAX_PATH] = "../../../input/"; // 设置默认目录
    ZeroMemory(&ofn, sizeof(ofn));
    ofn.lStructSize = sizeof(ofn);
    ofn.lpstrFilter = "OBJ Files (*.obj)\0*.obj\0All Files (*.*)\0*.*\0";
    ofn.lpstrFile = fileName;
    ofn.nMaxFile = sizeof(fileName);
    ofn.lpstrTitle = "Select Model File";
    ofn.lpstrInitialDir = initialDir; // 设置初始目录
    ofn.Flags = OFN_DONTADDTORECENT | OFN_FILEMUSTEXIST;
    
    if (GetOpenFileName(&ofn)) {
        return std::string(fileName);
    }
    
    return "";
}
#endif
#ifdef __APPLE__
#include <cstdio>
#include <iostream>
#include <string>

std::string openFileBrowser() {
    // 使用AppleScript调用macOS文件选择对话框
    std::string command = "osascript -e 'tell application \"Finder\" to set selectedFile to choose file with prompt \"选择OBJ模型文件\" of type {\"obj\"}' -e 'get POSIX path of selectedFile' 2>/dev/null";
    FILE* pipe = popen(command.c_str(), "r");
    if (!pipe) return "";

    char buffer[1024];
    std::string result = "";
    while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
        result += buffer;
    }
    pclose(pipe);

    // 移除换行符
    if (!result.empty() && result.back() == '\n') {
        result.pop_back();
    }
    return result;
}
#endif

void flipZ() {
    // Rotate mesh 180 deg about up-axis on startup
    glm::mat4x4 rot = glm::rotate(glm::mat4x4(1.0f), static_cast<float>(PI), glm::vec3(0, 1, 0));
    for (Vertex v : mesh->vertices()) {
        Vector3 vec = geometry->inputVertexPositions[v];
        glm::vec4 rvec = {vec[0], vec[1], vec[2], 1.0};
        rvec = rot * rvec;
        geometry->inputVertexPositions[v] = {rvec[0], rvec[1], rvec[2]};
    }
    psMesh->updateVertexPositions(geometry->inputVertexPositions);
}

void showSelected() {
    // pass
}


void redraw() {
    polyscope::requestRedraw();
}

//// 主曲率计算函数（在顶点处）
//std::tuple<VertexData<double>, VertexData<double>, VertexData<Vector2>, VertexData<Vector2>>
//computePrincipalCurvatures(VertexPositionGeometry& geometry) {
//    SurfaceMesh& mesh = *geometry.mesh;
//
//    VertexData<double> k1(mesh), k2(mesh);
//    VertexData<Vector2> dir1(mesh), dir2(mesh);
//
//    for (Vertex v : mesh.vertices()) {
//        // 获取顶点位置
//        Vector3 p0 = geometry.inputVertexPositions[v];
//
//        // 获取切空间基（用于投影到局部坐标系）
//        auto basis = geometry.vertexTangentBasis(v);
//        Vector3 e1 = basis.first;
//        Vector3 e2 = basis.second;
//
//        // 收集 1-ring 邻域点（在切平面投影）
//        std::vector<Vector2> neighbors;
//
//        for (Vertex nv : v.adjacentVertices()) {
//            Vector3 p = geometry.inputVertexPositions[nv];
//            Vector3 diff = p - p0;
//            // 投影到切平面
//            double x = dot(diff, e1);
//            double y = dot(diff, e2);
//            neighbors.push_back(Vector2{x, y});
//        }
//
//        if (neighbors.size() < 3) {
//            k1[v] = k2[v] = 0.0;
//            dir1[v] = Vector2{1.0, 0.0};
//            dir2[v] = Vector2{0.0, 1.0};
//            continue;
//        }
//
//        // 构建协方差矩阵（拟合局部高度场 z = f(x,y) 的 Hessian）
//        // 简化方法：用 PCA 拟合二次曲面
//        // 更准确方法：最小二乘拟合 ax² + bxy + cy²
//
//        // 这里使用一个简化的“邻域协方差法”
//        Eigen::Matrix2d cov = Eigen::Matrix2d::Zero();
//
//        for (const Vector2& pt : neighbors) {
//            cov(0,0) += pt.x * pt.x;
//            cov(0,1) += pt.x * pt.y;
//            cov(1,0) += pt.x * pt.y;
//            cov(1,1) += pt.y * pt.y;
//        }
//        cov /= neighbors.size();
//
//        // 特征值分解
//        Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> solver(cov);
//        Eigen::Vector2d eigenvalues = solver.eigenvalues();
//        Eigen::Matrix2d eigenvectors = solver.eigenvectors();
//
//        // 曲率近似为特征值的某种缩放（实际应拟合高度场二阶导数）
//        // 这里仅为可视化，用相对大小即可
//        double lambda1 = eigenvalues(0);
//        double lambda2 = eigenvalues(1);
//
//        // 简单近似：假设曲率正比于特征值倒数或直接使用
//        // 更准确做法：拟合局部二次曲面 z = ax² + bxy + cy²，然后 Hessian = [2a, b; b, 2c]
//        // 我们这里用一个启发式缩放：
//        k1[v] = lambda1 > lambda2 ? 1.0 / (lambda1 + 1e-8) : 1.0 / (lambda2 + 1e-8);
//        k2[v] = lambda1 > lambda2 ? 1.0 / (lambda2 + 1e-8) : 1.0 / (lambda1 + 1e-8);
//
//        // 主方向
//        dir1[v] = Vector2{eigenvectors(0,0), eigenvectors(1,0)};
//        dir2[v] = Vector2{eigenvectors(0,1), eigenvectors(1,1)};
//
//        // 可选：根据曲率符号调整（这里省略，仅用于可视化相对大小）
//    }
//
//    return std::make_tuple(k1, k2, dir1, dir2);
//}


void loadMesh(const std::string& filepath) {
    // Unregister the previous mesh if it exists
    if (psMesh) {
        polyscope::removeStructure(psMesh);
    }
    
    std::tie(mesh, geometry) = readManifoldSurfaceMesh(filepath);
    psMesh = polyscope::registerSurfaceMesh("Primal mesh", geometry->inputVertexPositions, mesh->getFaceVertexList(),
                                            polyscopePermutations(*mesh));

 /*   auto [k1, k2, dir1, dir2] = computePrincipalCurvatures(*geometry);
    psMesh->addVertexScalarQuantity("k1 (approx)", k1)->setEnabled(true);
    psMesh->addVertexScalarQuantity("k2 (approx)", k2);
    VertexData<Vector3> dir1_world(*mesh);
    for (Vertex v : mesh->vertices()) {
        auto basis = geometry->vertexTangentBasis(v);
        dir1_world[v] = dir1[v].x * basis.first + dir1[v].y * basis.second;
    }
    psMesh->addVertexVectorQuantity("Principal Dir 1", dir1_world)
          ->setEnabled(true)
          ->setVectorRadius(0.001)
          ->setVectorLengthScale(0.1);*/


    psMesh->setEnabled(true);
    
    currentFilePath = filepath;
    strcpy(filePathBuffer, filepath.c_str());
    flipZ();
}

void addCheckerboard(VertexData<Vector2>& flattening) {

    std::vector<std::array<double, 2>> V_uv(mesh->nVertices());
    for (Vertex v : mesh->vertices()) {
        V_uv[v.getIndex()] = {flattening[v][0], flattening[v][1]};
    }
    checkerboard = psMesh->addVertexParameterizationQuantity("checkerboard", V_uv);
    checkerboard->setEnabled(true);
    checkerboard->setStyle(polyscope::ParamVizStyle::CHECKER);
    checkerboard->setCheckerColors(std::make_pair(glm::vec3{1.0, 0.45, 0.0}, glm::vec3{0.55, 0.27, 0.07}));
    checkerboard->setCheckerSize(0.01);
}
VertexData<Vector3> mapToMeshData(VertexData<Vector2>& flattening) {

    VertexData<Vector3> flat = geometry->inputVertexPositions;
    for (Vertex v : mesh->vertices()) {
        flat[v] = Vector3{flattening[v][0], flattening[v][1], 0.0};
    }
    return flat;
}
// 展平算法函数实现（占位符）
void flattenMesh(FlattenMethod method) {
    SpectralConformalParameterization SCP = SpectralConformalParameterization(mesh.get(), geometry.get());

    switch (method) {
        case FlattenMethod::Tutte:
            SCP_FLATTENING = SCP.tutte_flatten();
            break;
        case FlattenMethod::SCP:
            SCP_FLATTENING = SCP.arap_flatten();            
            break;
        case FlattenMethod::LSCM:
            SCP_FLATTENING = SCP.lscm_flatten();
            break;
    }
    
    
    SCP_MESH = mapToMeshData(SCP_FLATTENING);
    addCheckerboard(SCP_FLATTENING);

    if (DISPLAY_FLAT) 
    {
        // Display SCP parameterization
        geometry->inputVertexPositions = SCP_MESH;
        geometry->normalize(CoM, true);
        polyscope::view::flyToHomeView();
        polyscope::view::style = polyscope::view::NavigateStyle::Planar;
    } else {
        // Display original
        geometry->inputVertexPositions = ORIGINAL;
        polyscope::view::style = polyscope::view::NavigateStyle::Turntable;
    }
    psMesh->updateVertexPositions(geometry->inputVertexPositions);
    redraw();
}

void functionCallback() {
    if (ImGui::Button("Select Model File")) {
#ifdef _WIN32
        std::string selectedFile = openFileBrowser();
        if (!selectedFile.empty()) {
            try {
                ::loadMesh(selectedFile);
                currentFilePath = selectedFile;
                strcpy(filePathBuffer, currentFilePath.c_str());
            } catch (const std::exception& e) {
                // Handle error
            }
        }
#elif __APPLE__
        std::string selectedFile = openFileBrowser();
        if (!selectedFile.empty()) {
            try {
                ::loadMesh(selectedFile);
                currentFilePath = selectedFile;
                strcpy(filePathBuffer, currentFilePath.c_str());
            } catch (const std::exception& e) {
                // Handle error
            }
        }
#else
        showFileDialog = true;
#endif
    }
        // 展平算法下拉列表
    const char* flattenMethodNames[] = { "Tutte", "LSCM","SCP" };
    
    int currentMethodIndex = static_cast<int>(currentFlattenMethod);
    if (ImGui::Combo("Flatten Method", &currentMethodIndex, flattenMethodNames, IM_ARRAYSIZE(flattenMethodNames))) {
        currentFlattenMethod = static_cast<FlattenMethod>(currentMethodIndex);
    }
    
    // 展平按钮
    if (ImGui::Checkbox("Flat", &DISPLAY_FLAT)) {
        flattenMesh(currentFlattenMethod);
    }
    if (showFileDialog) {
        ImGui::OpenPopup("Select Model File");
        showFileDialog = false;
    }
    
    if (ImGui::BeginPopupModal("Select Model File", NULL, ImGuiWindowFlags_AlwaysAutoResize)) {
        ImGui::Text("Current file: %s", currentFilePath.c_str());
        ImGui::InputText("File Path", filePathBuffer, IM_ARRAYSIZE(filePathBuffer));
        
        if (ImGui::Button("Load")) {
            try {
                ::loadMesh(std::string(filePathBuffer));
                ImGui::CloseCurrentPopup();
            } catch (const std::exception& e) {
                // Handle error - in a real application you might want to show this to the user
            }
        }
        
        ImGui::SameLine();
        if (ImGui::Button("Cancel")) {
            ImGui::CloseCurrentPopup();
        }
        
        ImGui::EndPopup();
    }
    

}




int main(int argc, char** argv) {


  // Configure the argument parser
    args::ArgumentParser parser("15-458 HW6");
    args::Positional<std::string> inputFilename(parser, "mesh", "A mesh file.");

    // Parse args
    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }

    // If a mesh name was not given, use default mesh.
    std::string filepath = "../../../input/Nefertiti_face.obj";
    if (inputFilename) {
        filepath = args::get(inputFilename);
    }



    // strcpy(filePathBuffer, currentFilePath.c_str());
    ::loadMesh(filepath);//currentFilePath
    ORIGINAL = geometry->inputVertexPositions;

    // // 计算主曲率
    //auto [k1, k2, dir1, dir2] = ::computePrincipalCurvatures(*geometry);

    polyscope::init();
    polyscope::state::userCallback = functionCallback;


    polyscope::show();

    return EXIT_SUCCESS;
}