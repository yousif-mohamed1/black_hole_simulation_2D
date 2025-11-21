#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <iostream>

int main()
{
    // --- Initialize GLFW ---
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW!\n";
        return -1;
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // MacOS requirement
#endif

    // --- Create a window ---
    GLFWwindow* window = glfwCreateWindow(800, 600, "OpenGL Test Window", NULL, NULL);
    if (!window) {
        std::cerr << "Failed to create GLFW window!\n";
        glfwTerminate();
        return -1;
    }

    glfwMakeContextCurrent(window);

    // Enable vsync (optional)
    glfwSwapInterval(1);

    // --- Initialize GLEW ---
    glewExperimental = GL_TRUE; // Needed for core profile to expose modern functions
    if (glewInit() != GLEW_OK) {
        std::cerr << "Failed to initialize GLEW!\n";
        glfwDestroyWindow(window);
        glfwTerminate();
        return -1;
    }

    // Query framebuffer size and set initial viewport
    int fbWidth = 0, fbHeight = 0;
    glfwGetFramebufferSize(window, &fbWidth, &fbHeight);
    glViewport(0, 0, fbWidth, fbHeight);

    std::cout << "GLFW + GLEW + GLM working successfully!\n";

    // Use GLM to test math
    glm::vec3 testVec(1.0f, 2.0f, 3.0f);
    std::cout << "GLM vector: "
        << testVec.x << ", "
        << testVec.y << ", "
        << testVec.z << "\n";

    // Render loop
    while (!glfwWindowShouldClose(window)) {
        // Handle dynamic resize
        glfwGetFramebufferSize(window, &fbWidth, &fbHeight);
        glViewport(0, 0, fbWidth, fbHeight);

        glClearColor(0.1f, 0.1f, 0.15f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
