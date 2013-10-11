_img.so: img_wrap.o
	gcc -dynamiclib -lpython *.o -o _img.so
img_wrap.o: img_wrap.c
	gcc -fPIC img_wrap.c  -o img_wrap.o -I/usr/include/python2.7/ -L/usr/lib/python2.7 -c
img_wrap.c: img.c img.i
	swig -python -o img_wrap.c img.i
clean:
	rm _img.so img_wrap.* img.py*
