ld2hdf5
=======

convert hapmap format ldfiles to hdf5 in java

This is a simple java program in a maven project using jhdf5. You still need to install the HDF5 system library.
to re-compile
   mvn clean package
   
to build a jar file.
set java.library.path to location of the native HDF5 library and run the compiled jar

to execute: java -Djava.library.path=<path_to_hdf5lib> -jar target/ld2hdf5-0.2.1-jar-with-dependencies.jar <make-rootdir with 1 2 .. X subdirs> <outfile.h5>
