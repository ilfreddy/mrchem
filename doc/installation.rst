============
Installation
============


-------------------
Build prerequisites
-------------------

The supplied setup script should be able to configure things
correctly, provided all the necessary tools are available.

Stallo
------

Using the Intel tool chain on Stallo::

    $ module load intel/2017a
    $ module load CMake/3.9.1
    $ module load Boost/1.63.0-intel-2017a-Python-2.7.13
    $ module load Eigen/3.3.4

Fram
----

Using the Intel tool chain on Fram::

    $ module load intel/2017a
    $ module load CMake/3.9.1
    $ module load Boost/1.63.0-intel-2017a-Python-2.7.13

Eigen is not available through the module system on Fram, so it must be
installed manually by the user, see the `Eigen3
<http://eigen.tuxfamily.org/index.php?title=Main_Page>`_ home page or the
``.ci/eigen.sh`` source file for details.

-------------------------------
Obtaining and building the code
-------------------------------

The public version of MRChem is available on GitHub::

    $ git clone git@github.com:MRChemSoft/mrchem.git

To build the code with only shared memory OpenMP parallelization::

    $ cd mrchem
    $ ./setup --prefix=<install-dir> --omp <build-dir>
    $ cd <build-dir>
    $ make
    $ make install

With the Intel tool chain on Stallo or Fram you need to specify the compilers
in the setup::

    $ ./setup --prefix=<install-dir> --omp --cxx=icpc <build-dir>

To build the code with hybrid MPI/OpenMP parallelization::

    $ ./setup --prefix=<install-dir> --omp --mpi --cxx=mpiicpc <build-dir>

The MRChem executables will be installed in ``<install-dir>/bin``.


-------
Testing
-------

A test suite is provided to make sure that everything compiled properly. To run
a collection of small unit tests::

    $ cd <build-dir>
    $ ctest -L unit

To run a couple of more involved integration tests::

    $ cd <build-dir>
    $ ctest -L integration

In order to run the integration tests you must first install the ``runtest``
Python package (listed as required in ``Pipfile``).

-------------------
Running the program
-------------------

A Python input parser will be provided along with the mrchem executable::

    $ build/bin/mrchem              // Python input parser
    $ build/bin/mrchem.x            // MRChem executable

The input parser takes a single file argument (default ``mrchem.inp``),
processes the input and calls the main executable. Output is written to stdout
but can be redirected to an output file::

    $ ./mrchem mrchem.inp > mrchem.out &

To run the program in OpenMP parallel use the environment variable
``OMP_NUM_THREADS`` (``unset OMP_NUM_THREADS`` will give you all threads
available, otherwise use ``export OMP_NUM_THREADS N``)::

    $ export OMP_NUM_THREADS 16
    $ ./mrchem mrchem.inp

When you run the program in hybrid MPI/OpenMP parallel, you must run the input
parser manually first with the dryrun ``-D`` option, before launching the main
executable with ``mpirun`` (or equivalent). For 20 threads each on 5 MPI
processes::

    $ ./mrchem -D mrchem.inp
    $ OMP_NUM_THREADS=20  mpirun -np 5 @mrchem.inp >mrchem.out &

