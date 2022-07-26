TRISURF NG
==========
Modified by Yoav based on Rajkumar's cluster version of Samo and Miha's trisurf-ng followed further developments with Samo, Mitja, Raj, and Shubhadeep

### 0. Diff
--------------
- vicsek interaction: setting vicsek_model=1, the force direction is average of normals in connected cluster, up to vicsek_radius, with weight vicsek_strength
- adhesion from Raj: step, parabolic, spherical and cylindrical
- added random_seed option to the tape: default (0) to this moment (linux epoch)
- gaussian curvature with angle sum formula $\frac{1}{A}\left( 2\pi-\sum\theta_i\right)$
- anisotropy (in progress): director $\hat{d}$ and shape operator energy $E(\hat{S})$
- option to have angle limits between triangles (to prevent spikiness) - optional to tape (default to -1 i.e. no effect, to use add e.g. "min_dihedral_angle_cosine=0.5")
- more vertex data is individualized as well as outputted to the .vtu 
    - type
    - force
    - normal
    - bending modulii
    - curvatures...
- added -flto=auto flag to configure.ac (link time optimization)
- thrown away multiprocessing options
- did some changes to spherical harmonics (some sort of refactoring out of ts_vtx?)
- did some changes to polymers (some sort of refactoring but can't even remember what)
- lots of unaccounted changed to bonds in order to make sure vertices and bonds and triangles remain ordered.  
    - vtx->neigh = [0, 1, 2, 3, 4, 5]  
    - vtx->bond = [{0,vtx}, {1,vtx}, {2,vtx}, {3,vtx}, {4,vtx}, {5,vtx}]  
    - vtx->tristar = [{0,1,vtx}, {1,2,vtx}, {2,3,vtx}, {3,4,vtx}, {4,5,vtx}, {5,0,vtx}]  
- lots of tabs and spaces changes due to viewing and working from VScode
- commandline --tape-options (-c) now works for most things and is saved on the .vtu \<tape\> section: use in the form
> ```$trisurf --tape-options adhesion_model=2,nshell=5,random_seed=9,iterations=10,#added_this_comment_at_end```

**Things to work on:**
- shape operator energy still has problems : we think it's the area, the bending constants in energy, and factors of 2.
    - do we want to add circumcenter to triangle and update with normal? triangle->circx
- spikiness still happens:  
    - probably because $\kappa$ is incorrect
- Add shear force
- Cross products into functions?
- constvol and constarea are almost certainly broken somehow, since they involve "side steps" which have not been updated
- initial distribution need more work for adding more types (partway there)  
- ghost (unmoving) and edge vertices are not working yet: energy doesn't work, bondflip is iffy, no initialization  
    - general rework of bonding and vertex order: to make a vtx-vtx bond there are several different functions which do different parts.
- making less huge vtu  
    - export less debug data or make it optional
    - however the heck compression works (in tape?)
    - somehow discard <trisurf> tag and use the vtk tags (like connectivity) instead: something with parseTrisurfTag
- keeping the same random seed sequence after restorations (currently it's resetting each time to the tape value/current time). might be okay?  
- in addition to lto, add pgo (profile guided optimization) to the compilation? 
    - I would guess it could be very effective, since the simulation is fairly constrained and we can run several very short but realistic examples
    - No clue how to mix it with autotools (makefile.am and configure.am)
- probably much much more  


### 1. Installation
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

### 2. Use
------

Prepare tape file, storing the definition for the simulation. You can use the sample tape file in the ``src/`` directory as a template for your simulation.

Run simulations with ``cluster-trisurf/src/trisurf --force-from-tape`` for initial run, or ``cluster-trisurf/src/trisurf`` for continuing aborted simulations.

======== LIBRARY VERSION ================

THis line seemed to fixed everything:
libtoolize --force


======= STATICALLY LINKING ==============
./configure --enable-static LDFLAGS=-static

