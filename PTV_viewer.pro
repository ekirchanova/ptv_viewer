TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXXFLAGS += -O3 -ffast-math

SOURCES += \
    GUI/gui_control.cpp \
    GUI/gui_display.cpp \
    GUI/gui_fields_common.cpp \
    GUI/gui_fields_iso.cpp \
    GUI/gui_geom.cpp \
    GUI/gui_init.cpp \
    GUI/gui_text.cpp \
    GUI/gui_textures.cpp \
    GUI/gui_cameras.cpp \
    IO/IO_bmp.cpp \
    IO/IO_ppm.cpp \
    IO/IO_tecplot.cpp \
    IO/IO_tga.cpp \
    IO/IO_vtk.cpp \
    image/image_opencv.cpp \
    Main.cpp



unix{
LIBS+=  -lGL -lGLU -lglut -lm
}
win32{
LIBS += -lopenGL32 -lGLU32 -lm
LIBS += -L$$PWD/my_lib -lglut32
INCLUDEPATH += "$$PWD/my_include/cv"
LIBS +=  "$$PWD/my_lib/cv/libopencv_*.dll"
}
DISTFILES += \
    my_lib/glut32.lib \
    my_lib/libpng.lib \
    my_lib/glut32.dll \
    my_lib/libpng12.dll \
    my_lib/zlib1.dll \
    test_in/ImperxCamera0_18_10_123/image0.bmp \
    test_in/ImperxCamera1_18_10_123/image0.bmp \
    test_in/ImperxCamera2_18_10_123/image0.bmp \
    test_in/ImperxCamera3_18_10_123/image0.bmp \
    test_in/ImperxCamera0_18_10_123/test.txt

HEADERS += \
    GUI/GUI.h \
    my_include/gl.h \
    my_include/glext.h \
    my_include/glu.h \
    my_include/glut.h \
    my_include/png.h \
    my_include/pngconf.h \
    IO/IO.h \
    data/basic.h \
    image/image.h
