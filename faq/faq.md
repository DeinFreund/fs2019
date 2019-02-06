# Frequently Asked Questions
This is a collection of questions that have been asked on the mailing list and
shall be used a resource of information.

_Please consult this FAQ pool first before you ask a question. If you still want to ask a question, please use piazza_

## Table of Contents
- [Eigen library on Euler](#the-eigen-library-on-the-euler-cluster)
- [Installing Python packages on Euler](#installing-python-packages-on-euler)
- [How to use an interactive node on Euler](#how-to-use-an-interactive-node-on-euler)
- [Daint tutorial](#daint-tutorial)
- [Git tutorial](#git-tutorial)

## The Eigen library on the Euler cluster

The Eigen library is a templated library and does not provide shared object
files but only header files.  In general, if you need additional libraries, you
should use
```
module load <package>
```
You can find the available `packages` by executing
```
module avail
```

In case of the Eigen library, the location of the header files is not added to
your environment when loading the library.  Consequently, the compiler can not
find the headers and will complain at compile time.  You can fix the problem by
specifying the include path manually
```
g++ -I/cluster/apps/eigen/3.2.1/x86_64/gcc_4.8.2/serial/include/eigen3 mycode.cpp
```
If you plan to use the Eigen library more often, you can set the search path of
the `gcc` compiler in your `$HOME/.bashrc` file by adding the line
```
export CPATH=/cluster/apps/eigen/3.2.1/x86_64/gcc_4.8.2/serial/include/eigen3
```
at the end of the file.  This will set the search path for both `gcc` and
`g++`.  But remember, this will always load the Eigen version `3.2.1`.  You
must manually update the path if you want to compile using another version.


## Installing Python packages on Euler

If you need to install a Python package on Euler, you must install it in your
`$HOME` due to write permissions.  For example, you can use the following
command to install the `beautifulsoup` package with python2.7:
```
python -m pip install --user bs4
```
After that you can use `import bs4` in your Python code.


## How to use an interactive node on Euler

Sometimes it is convenient to checkout an interactive node on Euler, rather
than submitting jobs to queue.  Interactive means that you will get back a
shell with the requested resources and you can run your code interactively.
You can request an interactive node on Euler with the following command:
```
bsub -n 24 -W 03:00 -Is /bin/bash
```
As you already know, `bsub` will place your request on the job scheduler queue.
Instead of executing your program, you request a BASH shell in which you will
work once the requested resources have been allocated for you.  The above
example allocates 24 cores for 3 hours.  When you submit this request, you will
have to wait until the job has been allocated for you.  If the cluster is busy,
you might have to wait a while until you receive your interactive shell.

If you want a Haswell node, you can specify it explicitly
```
bsub -n 24 -W 03:00 -R "select[model==XeonE5_2680v3]" -Is /bin/bash
```





## Working with Piz Daint supercomputer

You will be granted access to the main supercomputer in the CSCS (Swiss Supercomputer Center in Lugano): Piz Daint. The credentials (a username and a password) will be distributed to you when time comes.
Piz Daint consists of two parts: a GPU-enabled cluster Cray XC50 and a multicore cluster Cray XC40 (former Piz Dora). You will only have to use the GPU part, that offers 5320 nodes with an Intel R Xeon E5-2690 v3 CPU and an NVIDIAR Tesla P100 GPU.
Detailed information about the machines and how to use them can be found on the [CSCS web site](https://www.cscs.ch/computers/dismissed/piz-daint-piz-dora/m) as well as the [user portal](https://user.cscs.ch/)







