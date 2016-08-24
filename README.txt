Due to file size limitations of github repositories the example conditioned
model data files are not included. To obtain the necessary data files go to,

https://cloudstor.aarnet.edu.au/plus/index.php/s/j9CkRGRoVu4PiJk

then download and store the BET.zip and YFT.zip files in the OMconditioning
folder of the MSE-IO-BET-YFT github source tree.

To run the MSE-IO-BET-YFT example models first requires that the conditioned
model data is recreated in the correct hierarchy. To extract the data run the
setup-examples.bat (windows) or setup-examples.sh (linux / Mac) script in the
OMconditioning folder. To run the example models run a 64 bit version of R,
change directory to the root folder for the project, then run the MSE's with
the command,

source('RScripts/YFTandBET_demo_ script.R', print.eval=TRUE)

Note that to run on Mac OS X or Linux it is necessary to build the 64 bit Mseom
shared library (libMseom.so) and copy it to the lib sub-directory. A 64 bit
Windows dll is only provided with this distribution at this point in time.

To build the shared library on OS X, Linux, Cygwin or mingWin you must have
automake and libtool installed. You can build the Mseom shared library for
OS X or Linux by downloading the ADT project from https://github.com/pjumppanen/ADT ,
then from a command prompt, change directory to the root folder for the project
and execute the following shell commands.

chmod 755 ./*
./init_config.sh
./configure --disable-tapenade --with-r-lib-path={path to the R shared library import lib} --with-r-include-path={path to the R include files}

After configure script has run successfully build the library using the
shell commands,

make clean
make

then copy the libMseom.so file from the folder adt/objs/sample/Mseom/.libs
to the lib folder under the MSE-IO-BET-YFT project.
