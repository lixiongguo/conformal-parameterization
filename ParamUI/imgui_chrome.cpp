#include "imgui_chrome.h"

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include <GLFW/glfw3.h>

#include <filesystem>

void imguiInit(GLFWwindow* window) {
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGui::StyleColorsLight();
    ImGuiIO& io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_NoMouseCursorChange;
    ImGui_ImplGlfw_InitForOpenGL(window, false);
    ImGui_ImplOpenGL3_Init("#version 330 core");

    float xscale = 1.f;
    float yscale = 1.f;
    glfwGetWindowContentScale(window, &xscale, &yscale);
    const float size = 16.f * 0.5f * (xscale + yscale);
    const char* paths[] = {
        "/System/Library/Fonts/Hiragino Sans GB.ttc",
        "/System/Library/Fonts/STHeiti Medium.ttc",
        "/Library/Fonts/Arial Unicode.ttf",
        "/System/Library/Fonts/Supplemental/Arial Unicode.ttf",
    };
    ImFontConfig config;
    config.PixelSnapH = true;
    static const ImWchar ranges[] = {
        0x0020, 0x00FF, 0x0370, 0x03FF, 0x2000, 0x206F, 0x2200, 0x22FF, 0x3000, 0x30FF, 0x4E00, 0x9FFF, 0xFF00, 0xFFEF, 0,
    };
    for (const char* path : paths) {
        std::error_code error;
        if (!std::filesystem::exists(path, error)) continue;
        if (io.Fonts->AddFontFromFileTTF(path, size, &config, ranges)) return;
    }
    io.Fonts->AddFontDefault();
}

void imguiShutdown() {
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
}

void imguiBegin() {
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
}

void imguiEnd() {
    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}

bool imguiCapturesMouse() {
    return ImGui::GetCurrentContext() != nullptr && ImGui::GetIO().WantCaptureMouse;
}

bool imguiCapturesKeyboard() {
    return ImGui::GetCurrentContext() != nullptr && ImGui::GetIO().WantCaptureKeyboard;
}

void imguiMouseButton(GLFWwindow* window, int button, int action, int mods) {
    ImGui_ImplGlfw_MouseButtonCallback(window, button, action, mods);
}

void imguiScroll(GLFWwindow* window, double xoffset, double yoffset) {
    ImGui_ImplGlfw_ScrollCallback(window, xoffset, yoffset);
}

void imguiKey(GLFWwindow* window, int key, int scancode, int action, int mods) {
    ImGui_ImplGlfw_KeyCallback(window, key, scancode, action, mods);
}

void imguiChar(GLFWwindow* window, unsigned int codepoint) {
    ImGui_ImplGlfw_CharCallback(window, codepoint);
}
