TRISURF NG
==========
Modified by Yoav based on Rajkumar's cluster version of Samo and Miha's trisurf-ng

0. Diff
--------------
* vicsek interaction: setting vicsek_model=1, the force now sums normals on the 'connected cluster up to vicsek_radius' with weight vicsek_strength
* adhesion works as it did, does not conflict: turn off vicsek for normal adhesion (vicsek_model=0)
* added random_seed option to the tape: default (0) to epoch
* commandline --tape-options now works for adhesion and vicsek parameter too and is saved on the .vtu \<tape\>: use in the form
> $trisurf --tape-options vicsek_model=0,adhesion_model=2,nshell=5,random_seed=9
* gaussian curvature with angle sum formula
* anisotropy (in progress)
* more vertex data is outputted to the .vtu (type, force, normal, curvatures...)


1. Installation
--------------

To compile the program, user must have ``automake``, and ``gcc`` tools installed on the computer.

Required C libraries are:
* libconfuse
* libgsl
* libxml2
* zlib

and python libraries for running tsmgr:
* tabulate
* psutil
* configobj

On Debian based systems, install prerequisities by typing the following command in the command line:

``sudo apt-get install libconfuse-dev libgsl0-dev libxml2-dev zlib1g-dev automake gcc python3-psutil python3 python3-pip``
``sudo pip3 install tabulate configobj``

Move to the project root directory and compile with:

``./configure``

``make``

``make install``

If you are experiencing difficulties due to different automake versions, proceed with the longer procedure:

``aclocal``

``autoheader``

``automake -ac``

``autoconf``

``autoreconf -if``

``./configure``

``make``

``make install``


This procedure can be done automatically by calling the build.sh script.
(build.sh is modified not to make install)

2. Use
------

Prepare tape file, storing the definition for the simulation. You can use the sample tape file in the ``src/`` directory as a template for your simulation.

Run simulations with ``cluster-trisurf/src/trisurf --force-from-tape`` for initial run, or ``cluster-trisurf/src/trisurf`` for continuing aborted simulations.

======== LIBRARY VERSION ================

THis line seemed to fixed everything:
libtoolize --force


======= STATICALLY LINKING ==============
./configure --enable-static LDFLAGS=-static

