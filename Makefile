CC=g++
CXXARGS=-std=c++11 -fmax-errors=3 -lasan -Wno-narrowing # -Wpedantic
LGRAPH="-lgraph -Wl,--rpath -Wl,/usr/local/lib"
TARGET=./chips

DIR="/media/$$USER/Elements/Anton/CODING/chips"
# DIR="."
X11MOUSE="-lX11"
IMAGE=image

INC=-I/media/$$USER/Elements/programs/include
GLUT=-lglut -lGL
OPEN=xdg-open
DEPENDS=*.hpp *.cpp
DETATCH=gnome-terminal --working-directory=$(DIR) --

chips: $(DEPENDS)
	$(CC) $(INC) main.cpp -o chips $(CXXARGS) # 2> gpp.output

test: chips
	$(TARGET) $(IMAGE) #&> $(TARGET).output
	convert $(IMAGE).ppm $(IMAGE).jpg
	rm ./$(IMAGE).ppm
	$(OPEN) $(IMAGE).jpg

clear:
	rm -f *.basic_canvas

clean:
	rm -f *.output
