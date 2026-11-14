#pragma once

struct GLFWwindow;

void imguiInit(GLFWwindow* window);
void imguiShutdown();
void imguiBegin();
void imguiEnd();
bool imguiCapturesMouse();
bool imguiCapturesKeyboard();
void imguiMouseButton(GLFWwindow* window, int button, int action, int mods);
void imguiScroll(GLFWwindow* window, double xoffset, double yoffset);
void imguiKey(GLFWwindow* window, int key, int scancode, int action, int mods);
void imguiChar(GLFWwindow* window, unsigned int codepoint);
