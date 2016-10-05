GMES (GIST Maxwell's Equations Solver) is a free Python package
for the electromagnetic simulation using FDTD method. This file 
describes how to compile and set up GMES.

1. Download Latest Version
You can download the latest release at the GMES project page of 
sourceforge:

http://sourceforge.net/projects/gmes

The latest development version is always available from our code 
repository:

http://sourceforge.net/p/gmes/code-git/

2. Compatibility
GMES follows the C++11 standard and is compatible to Python 2.7, 
thus it is OS portable in source code level. However, I usually 
execute GMES only on linux, especially Ubuntu. If you have 
succeeded to compile and execute on your environment, please 
report your setup and environment. Any contributions are welcome.

2.1. GNU/Linux
GMES has been tested on Ubuntu 14.04 LTS Trusty Tahr 64 bit.

3. Preliminaries

3.1. GNU/Linux
GMES uses various Python and non-Python programs. You will need 
the following programs to install and execute GMES. The following 
programs are provided as packages on Ubuntu distributions.

* python
* g++
* swig
* cython
* numpy (Ubuntu package: python-numpy)
* scipy (Ubuntu package: python-scipy)
* matplotlib (Ubuntu package: python-matplotlib)
* MPI for Python (optional) (Ubuntu package: python-mpi4py) 
* pytables (optional) (Ubuntu package: python-tables)

4. Compile and Install
If you like to test the GMES without or before installation, 
follow 'compile for the inplace test' section. Otherwise, follow 
'compile' and 'install' sections in order. 

4.1. For GNU/Linux

4.1.1. compile for the inplace test
$ python setup.py build_ext --inplace

4.1.2. compile
$ python setup.py build_ext --swig-opts="-c++"

4.1.3. install
$ sudo python setup.py install

5. Testing
If you completed the installation or compiled with inplace option,
you can test the proper execution of GMES using the unit tests in 
the 'tests' directory. If you see 'OK' without any error, it's OK.

$ python <test.py>

6. Execution Commands

6.1. GNU/Linux with MPI2
Your MPI2 implementation may requires pre- or post- commands to 
exploit the parallel environment. Refer the manuals of your MPI2 
implementation for these commands. <N> and <GMES SCRIPT> indicate 
the number of parallel process and your script which uses the GMES
package.
 
$ mpiexec -l -n <N> python-mpi <GMES SCRIPT>

Also, in the parallel environment, a MPI-enabled Python 
interpreter is required. Check the MPI for Python user manual at 
http://mpi4py.scipy.org/docs/usrman/.

7. User Manual
The current version of GMES does not have a separate user manual. 
However, you can start from examples in the 'examples' and 'tests'
directories. If you need more detailed documentations on the 
classes and parameters of GMES, please use the 'dir' and 'help'
command in the Python interpreter.

8. Participation
Please report any errors or submit patches (if possible)! Any 
contributions are welcome. If you would like to join the GMES 
development team or contribute to the distribution, please contact
the gmes-user mailing list, details at 
http://sourceforge.net/p/gmes/mailman/.

9. Known Issues

9.1. numpy.inf
Please use a large number instead of numpy.inf since GMES does not
correctly treat the value of numpy.inf as infinity.

9.2. Warnings on Compilation
Cython sets '-Wunused-but-set-variable' compile option and it 
generates the following warnings during compilation.

src/material.c: In function ‘__pyx_pf_4gmes_8material_7Lorentz_10get_pw_material_hz’:
src/material.c:33412:13: warning: variable ‘__pyx_v_coords’ set but not used [-Wunused-but-set-variable]

These warnings are not harmful and you can safely ignore them.

9.3. Python 3 Support
GMES has not tested with Python 3.

 -- Kyungwon Chun <kwchun@gist.ac.kr>
