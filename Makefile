CC=g++
CPPFLAGS=-std=c++11 -fmax-errors=3 # -lasan -Wno-narrowing # -Wpedantic
LGRAPH="-lgraph -Wl,--rpath -Wl,/usr/local/lib"
TARGET=./chips

DIR="/media/$$USER/Elements/Anton/CODING/chips"
# DIR="."
X11MOUSE="-lX11"
IMAGE=image

GLUT=-lglut -lGL
OPEN=xdg-open
DEPENDS=*.hpp *.cpp
OBJECTS= \
	Chip.o Edge.o Node.o Face.o \
	canvas.include.o \
	include/canvas/Vector.o \
	include/canvas/color.o \
	include/canvas/kernel.o \
	include/canvas/BasicCanvas.o \
	include/canvas/LayeredCanvas.o \
	include/spline/Point.o \
	include/spline/Spline.o \
	include/spline/SplineConstruct.o \

DETATCH=gnome-terminal --working-directory=$(DIR) --

chips: $(DEPENDS) main.o $(OBJECTS)
	$(CC) $(CPPFLAGS) -o chips main.o $(OBJECTS)

run: chips
	$(TARGET) $(IMAGE) --chip --color
	$(OPEN) $(IMAGE).ppm

test: test.o $(OBJECTS)
	$(CC) $(CPPFLAGS) -o test test.o $(OBJECTS)

run_test: test
	./test | grep -v "request for out-of-frame pixel"

clear:
	rm -f *.ppm *.basic_canvas *.basic_canvas_dump

clean:
	rm -f *.output *.o *.gch

redo: clean
	cd include/canvas/; make clean
	cd include/spline/; make clean
	make chips

all_bg: chips
	./chips C_black --chip --black --color && convert C_black.ppm C_black.jpg
	./chips C_green --chip --green --color && convert C_green.ppm C_green.jpg
	./chips C_rainbow --chip --rainbow --color && convert C_rainbow.ppm C_rainbow.jpg
