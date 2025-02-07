# compiler and compilation flags
CC=gcc
CFLAGS = -D_GNU_SOURCE -Wall -Winline -std=c99 -g -O2
EXTRA = -DBUILD_DIVSUFSORT64 -DHAVE_CONFIG_H

OBJS = sssort.o trsort.o utils.o divsufsort.o

# divsufsort 3.0.1
divsufsort64.a: ${OBJS}
	ar rcs $@ ${OBJS}

# pattern rule for all objects files
%.o: %.c ${HEADERS} 
	$(CC) -c $(CFLAGS) -o $@ $(EXTRA) $<


clean: 
	rm -f *.o *.a
            





