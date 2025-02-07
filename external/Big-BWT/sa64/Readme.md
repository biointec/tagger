
## A libdivsufsort64 clone 

This is a stand alone version of libdivsufsort64 the 64-bit version of the [libdivsufsort](https://github.com/y-256/libdivsufsort) suffix sorting library by Yuta Mori. By default the standard library available on github builds a (shared) 32-bit version that works for inputs of size at most 2^31 -2; to avoid changing the cmake configuration file this version builds directly a static 64-bit version. 

It has been obtained from Yuta Mori's libdivsufsort library version 2.0.2-1 extracting the files

	config.h
	divsufsort.h	
	divsufsort64.h	
	divsufsort_private.h
	divsufsort.c 
	sssort.c
	trsort.c
	utils.c

and adding a makefile that builds the static libdivsufsort64.a library


