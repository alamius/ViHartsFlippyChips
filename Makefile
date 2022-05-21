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
	main.o Chip.o \
	Edge.o Node.o Face.o \
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

chips: $(DEPENDS) $(OBJECTS)
	$(CC) -o chips $(OBJECTS)

test: chips
	$(TARGET) --chip --color $(IMAGE)
	convert $(IMAGE).ppm $(IMAGE).jpg
	rm ./$(IMAGE).ppm
	$(OPEN) $(IMAGE).jpg

clear:
	rm -f *.basic_canvas

clean:
	rm -f *.output *.o *.gch
	cd include/canvas/; make clean
	cd include/spline/; make clean
