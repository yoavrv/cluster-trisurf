TRISURF NG
==========
A Monte-Carlo Triangulated-Surface simulator for membrane dynamics.  

At the core, this simulates of the bending energy due to curvature of a 2D membrane, based on the Helfrich Hamiltonian  

$$ \iint \frac{\kappa}{2}\left(C_1+C_2-C_0\right)^2 dA $$

Where the membrane is described by a graph of vertices/nodes $i\in V$, bonds/edges $\left\langle i,j\right\rangle \in E$, and triangles $\left\langle i,j,k\right\rangle\in T$.




### 1. Installation
--------------

To compile the program, user must have ``automake``, and ``gcc`` tools installed on the computer.

Required C libraries are:
* libconfuse
* libgsl
* libxml2
* zlib


On Debian based systems, install prerequisities by typing the following command in the command line:

``sudo apt-get install libconfuse-dev libgsl0-dev libxml2-dev zlib1g-dev automake gcc``

Move to the project root directory and compile with:

``./configure``  
``make``  
``make install``  

If you are experiencing difficulties due to different automake versions, proceed with the longer procedure:
 
``aclocal``  
``autoheader``  
``automake -ac``  
``autoconf``  
``./configure``  
``make``  
``make install``  

This procedure can be done automatically by calling the build.sh script.  
build_no_install.sh is modified not to install, only making `src/trisurf` in the project directory.

Hint: To install different versions on the cluster  
``./configure --prefix=~/apps --program_prefix=special_``  
``make install`` will install a special_trisurf in ~/apps.  
Check that it is correct with  
``make install --dry-run``

### 2. Use
------

Prepare tape file, storing the definition for the simulation. You can use the sample tape file in the ``src/`` directory as a template for your simulation.

Run simulations with ``trisurf --force-from-tape`` for initial run, or ``cluster-trisurf/src/trisurf`` for continuing stopped simulations from binary ``dump.bin``.  
Run ``trisurf --restore-from-vtk some_file.vtu`` to continue simulation from a VTU file.  
see ``trisurf --help`` for more information.

======== LIBRARY VERSION ================  

This line seemed to fixed everything:  
``libtoolize --force``


======= STATICALLY LINKING ==============  
(Not vetted on the new version)  

``./configure --enable-static LDFLAGS=-static``


### 3. Difference from previous versions of Samo
--------------------------------------------------
- vicsek interaction: setting force_model=16, the force direction is average of normals in connected cluster, up to vicsek_radius, with weight vicsek_strength
- Adhesion from Raj-kumar Sadhu: step, parabolic energy from surface of plane, spherical or cylindrical geometry.
- Added random_seed option to the tape fro reproducability: default (0) to this moment (linux epoch)
- Gaussian curvature with angle sum formula $\frac{1}{A}\left( 2\pi-\sum\theta_i\right)$. Note gaussian modulus is negative!
- Anisotropy: vertices have director $\hat{d}$ and shape operator energy $E(\hat{S})$, as well as spontaneous deviator $d_0$
- Dihedral angle limits between triangles to prevent "spikiness" pathologys (in tape  "min_dihedral_angle_cosine=0.1". value of -1 means no effect)
- More vertex data is individualized per-vertex and outputted to the .vtu 
    - type
    - force
    - normal
    - bending modulii
    - curvatures
- Added -flto=auto flag to configure.ac (link time optimization)
- thrown away multiprocessing options
- Refactror spherical harmonic out of ts_vtx
- Refactors to polymers to disentangle from other parts of the code (not sure if it still works!)
- Major changes to bonds in order to make sure vertices and bonds and triangles remain ordered.  
    - vtx->neigh = [0, 1, 2, 3, 4, 5]  
    - vtx->bond = [{0,vtx}, {1,vtx}, {2,vtx}, {3,vtx}, {4,vtx}, {5,vtx}]  
    - vtx->tristar = [{0,1,vtx}, {1,2,vtx}, {2,3,vtx}, {3,4,vtx}, {4,5,vtx}, {5,0,vtx}]  
- Commandline --tape-options (-c) now works for most things and is saved on the .vtu \<tape\> section: use in the form
> ```$trisurf --tape-options adhesion_model=2,nshell=5,random_seed=9,iterations=10,#added_this_comment_at_end```

**Things to work on:**
- Figure out all the factors of 2 and -1 for consistency with theory.
- Shape operator is only approximately the same as the Helfrich energy and Gaussian curvature
- Add shear force
- Cross products into functions?
- constvol and constarea are broken, since they involve "side steps" which have not been updated
- Initial distribution need more work for adding more types (partway there)  
- ghost (unmoving) and edge vertices are not working yet: energy doesn't work, bondflip is iffy, no initialization  
    - general rework of bonding and vertex order: to make a vtx-vtx bond there are several different functions which do different parts.
- Making smaller vtu  
    - export less debug data or make it optional
    - Re-enable compression (in tape?)
    - Somehow discard <trisurf> tag and use the vtk tags (like connectivity) instead: Probably modify parseTrisurfTag
- ? Keep the same random seed sequence after restorations (currently it's resetting each time to the tape value/current time). 
- ? in addition to lto, add pgo (profile guided optimization) to the compilation
    - I would guess it could be very effective, since the simulation is fairly constrained and we can run several very short but realistic examples
    - No clue how to mix it with autotools (makefile.am and configure.ac)
 
