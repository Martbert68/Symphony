CC=gcc
render: render.c  martin.o martin.h
	$(CC) -Ofast -ggdb render.c martin.o -o render -lX11 -lm -ljpeg 

tree: tree.c  martin.o martin.h
	$(CC) -Ofast -ggdb tree.c martin.o -o tree -lX11 -lm -ljpeg 

martin.o: martin.c martin.h
	gcc -Ofast martin.c -c -o martin.o 

out.mp4: jpegs/kmage0001.jpg
	#ffmpeg -y -r 30 -t 240 -i jpegs/jmage%04d.jpg -i Beginnings.flac out.mp4
	#ffmpeg -y -r 10 -t 136 -i jpegs/jmage%04d.jpg -vf "transpose=1" out.mp4
	#ffmpeg -y -r 60 -t 60 -i jpegs/kmage%04d.jpg out.mp4
	ffmpeg -r 60 -i jpegs/kmage%04d.jpg -vf scale=-2:720 -c:v libx264 -profile:v main -level:v 3.0 -x264-params scenecut=0:open_gop=0:min-keyint=72:keyint=72 -c:a aac -preset slow -crf 23 -r 30 -sn -f mp4 out.mp4
all: out.mp4

