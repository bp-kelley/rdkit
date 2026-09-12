# Fragment based Conformation generation

This conformer generated was implemented as a test bed for using different
techniques for rotor-driven approaches.  This is a tried and true approach
for many different conformation generators:

    References:
      https://pubs.acs.org/doi/10.1021/ci100031x (omega)
      https://pmc.ncbi.nlm.nih.gov/articles/PMC2896087/ (frog2)
      https://pubs.acs.org/doi/10.1021/acs.jcim.3c00563 (conforge)


Why another?

  1. Licensing.  BSD is friendly for both academic and commercial licenses.
     CONFORGE is a good piece of kit, but uses the Hamburg Torsion Library
     which appears to be non-free for commercial use.

  2. The RDKit has, in my opinion, the hardest pieces required for rigid
     conformer generation.
     
       * MMFF forcefield,
       * a very good DG implementation
       * and, surprisingly, it's own torsion library in the form of ETKDG.

  3. Adds a playground for search strategies, Systematic, ThompsonSampling, Tree and
     (the default) combination of both.

  4. Has the ability to use non-mmff forcefields such as UFF for non parameterized
     atoms.  (Not implemented yet)

  5. Finally, the original goal was to make fast Synthon 3D searching which requires
     an easily modifiable conformation generating core.

Additions to the RDKit Core:

  ChemReactions/Enumerate/EnumerateSynthons.h/.cpp - this adds quick synthon zipping
    to enumerate synthon libraries or, optionally, combine exact reagents.

  ForceFieldHelpers/MMFF - Adds the _TOR varient from:
    J. Wahl, J. Freyss, M. von Korff, T. Sander, J. Cheminform. 2019, 11, 53
    (doi:10.1186/s13321-019-0371-6), which mainly correct the overestimated
    aromatic C-N amide rotation barrier.

    note: this doesn't appear to help, so we don't technically need it.
    
Core Components
===============

Conformer Generation
--------------------

Conformer generation has the following components:

  1. A fraglib that generates initial conformations of fragments.
  2. A set of torsion angles to sample.
  3. A joiner that forms the initial conformations from the fragments.
  4. A rotor driving core generating a conformation given a set or torsion angles.
  5. A search strategy to find the minimum score (MMFF energy)


### 1. Fraglib

The fragment library uses distance geometry to generate the initial set of fragments.
The current strategy employs ETKDG + MMFF minimiation.  Each fragment is assigned
a fragment class, i.e. SmallRing, LargeRing, Rigid etc and computed with varying
parameters to find optimal conformations:

| fragment   | min     | max     | maxConfs | eWindow | rmsd | conf count  | confs/rotor |
| class      | sampled | sampled | output   |         |      | operator    |             |
|------------|--------:|--------:|---------:|--------:|-----:|-------------|------------:|
| Rigid      | 30      | 30      | 100      | 8       | 0.1  | MAXIMUM     | 0           |
| SmallRing  | 30      | 800     | 800      | 8       | 0.1  | MULTIPLY    | 6           |
| LargeRing  | 100     | 1600    | 1000     | 24      | 0.1  | EXPONENTIAL | 10          |
| Acyclic    | 5       | 50      | **1**    | 8       | 0.1  | EXPONENTIAL | 10          |
| Exhaustive | 30      | 1000    | 100      | 12      | 0.1  | MAXIMUM     | 0           |
| Fast       | 10      | 100     | 20       | 8       | 0.1  | MAXIMUM     | 0           |

Of note, the operator of Maximum computes rotor*confs/rotor for the desired number of samples.
The operator of expoential computes 2^confs/rotor to search floppier rings.  For
minimization, exit vectors (dummy atoms) are replaced with carbons to generate realistic
geometries.

### 2. A set of torsions to sample

This package supports multiple torsion angle generators.

  * ETKDG torsion sampling.
  * Uniform sampling, thorough but quite slow
  * Smarts style library files like the Hamburg Torsion Library (not supplied for licensing reasons).
  * Smirnoff FF sampling - requires the OpenFF SMIRNOFF parameter files

The default is ETKDG torsion sampling with a uniform backup for missing torsions.

One novelty is symmetric torsions are reduced via graph theory (canonicalization)
not hand curated smarts patterns so it is hopefully more robust.
See Samplers/TorsionSymmetries.h for details


### 3. Joining Fragments

Each fragment is embedded with an "exit vector" so the joiner just zips these together
preserving the exit vector and ensuring that the bond length is sane.  For each join,
the preferred torsion angle is chosen from the torsion sampler.  The output is
a rotor driving context for the setup.  This has the appropriate scoring functions
for the molecule prebuilt.

### 4. Rotor Driving

The rotor driver is designed to keep the largest fragment fixed and sorts for rotors
from inside out.  This means that earlier rotors will move more atoms and have
larger changes in conformation.

The rotor driver is given a FragmentJoinerContext that can set rotors and
return the current score. Since each fragment has it's own internal energy computed,
this means that only contributions across fragments needs to be computed.  The
default scorer also clamps VdW distances to avoid unnecessary computations.


### 5. 

Now that we have rotor driving and scoring, the problem becomes a search problem.

| method     | Description |
|------------|-------------|
| Tree       | sets one rotor at a time, widest-moving first, keeping only the best N partial conformers at each step |
| Systematic | builds up combinations of fragment confs and angles, keeping everything under the energy window rather than a fixed count |
| Thompson   | bandit sampling over joint (fragment conf, rotor angle) arms, energy as reward, budget scaling with rotor and conf count |
| Merge      | runs systematic and thompson and merges the ensembles: they solve largely different molecules, so this measures how much of their union survives ranking |
| Auto       | either tree or thompson switching to systematic above N rotors |


### Benchmark

There are two timings, no fraglib (cold) and precomputed (warm.  Methods shown
are tree and thompson sampling.

### Platinum (2,753)

| mode | RigidRotorSearchMode | %<1Å | %<2Å | cold ms/mol | warm ms/mol | confs |
|---|---|---:|---:|---:|---:|---:|
| tree | Tree | 73.1% | 95.9% | 251 | 96 | 34 |
| auto (tree) | Auto + autoBudget=false | 74.1% | 96.3% | 328 | 177 | 37 |
| systematic | Systematic | 82.2% | 98.0% | 470 | 341 | 72 |
| thompson | Thompson | 89.0% | 99.2% | 356 | 202 | 109 |
| auto (ts) | Auto + autoBudget=true | 89.3% | 99.2% | 413 | 185 | 108 |
| merged | Merged | 89.7% | 99.4% | 701 | 547 | 118 |

### PDBbind (3,344)

| mode | RigidRotorSearchMode | %<1Å | %<2Å | cold ms/mol | warm ms/mol | confs |
|---|---|---:|---:|---:|---:|---:|
| tree | Tree | 72.1% | 95.8% | 273 | 113 | 46 |
| auto (tree) | Auto + autoBudget=false | 73.8% | 96.5% | 445 | 326 | 53 |
| systematic | Systematic | 81.2% | 97.8% | 615 | 513 | 80 |
| thompson | Thompson | 85.1% | 98.7% | 434 | 268 | 114 |
| auto (ts) | Auto + autoBudget=true | 85.2% | 98.6% | 548 | 342 | 112 |
| merged | Merged | 86.5% | 98.9% | 908 | 794 | 122 |

Code organization
==================
This is a brief overview of the organization of the conformer generation:

 * Confgen/FragmentConfGen.h(.cpp)h - core class for creating conformations.

   Usage:
     FragConfGen gen;
     auto result = gen.build(mol):;

   The FragConfGen uses a thread-safe fragment library that stores the
   ETKDG embedded fragments for the builds.  When threading, always
   share the fraglib.

   Fragment libraries can be pregenerated, we don't currently ship one
   with the RDKit for size issues.  see genFragLib for options.

 * Fraglib/Fraglib.h(.cpp)

  Fragment library.  This holds the embedded fragments used when
  rigidly rotating.

 * Joiner/FragmentJoiner.h(.cpp)

  This generates the initial rigid fragments with appropriate bond lengths
  and returns the FragmentJoinerContext that the searches will use
  when rotor driving.

 * Search/RigidRotorSearch.h(.cpp)

  This is the base class for rotor driving.  It takes a FragmentJoinerContext
   and searches it.


3D Synthon Searching
--------------------

The initial reason for developing this was to enable synthon searching.  The original
implementation was simply driving junction bonds of the joined rotors, but the
question was asked, "was this good enough" which drove the development of the
full rotor driving algorithm and, to be fair, spiraled a bit out of hand.


 * SynthonSearch/EnumerateSynthon3D.h(.cpp) - a variant of ChemReactions/Enumearte
   that generates conformations on demand for a synthon library.

   This can operate in two modes:
       Coarse:  embed as much as the synthon as possible and only rotate one rotor
                per synthon.
       Full:    perform the full conformer search

   It has a few extra features such as quickly finding the number of atoms
    in a product or the number of atoms added when using a particular synthon
    
 * SynthonSearch/SynthonSearch3D.h

   The actual search.  The search is staightforward, run through every
   synthon, grabbing a random few products and rank the synthons
   based on score.

   Sort the best synthons and use them to refine a final result.
   Using multiple trajectories increases the probability of finding
   a result.

Along the way, there were some interesting issues to tackle:

  1. What happens when combined synthons form new rings?
  2. Is a greddy algorithm enough?  Are synthons additive in shape/color
     space?
  3. What is a good result when you can't or don't want to search the
     whole space?

BENCHMARKS: (need to put here)

Note:

 This all came to pass because Matt Stahl asked "Why not use Thompson sampling
 as a search strategy?" which proves, yet again, he's smarter than I am.

 Claude was used for writing the benchmarks and tests, reviewing for bugs and
  the (quite often) large scale refactoring of code.
 Claude also ported the QCP RMSD code from: https://github.com/pandegroup/IRMSD

 
 I had naively assumed I could slap together a framework
 for doing rigid rotation techniques quickly.  This was pretty far from what happened.
 The majority of the time went into making an efficient and decomposable core to
 become a playground for implementation of search algorithms.  Scoring functions
 were validated, rotor driving was optimized and so on.

 Similar to what Pat Walters has written, this proved immensely beneficial for A/B testing
 and writing regression tests (which I'm not terribly great at)
 where Claude was invaluable.  On the core code, understanding the complexities of
 cheminformatics and general organization, not so much and was abandonded early in
 development.

 The thompson sampling and other search code is pretty bare bones but shows real
 promise for enhancements as I learn more about what priors help.

 The next steps are:
 
    1. add invertible nitrogens
    2. more search optimizations
    3. macrocycles, perhaps we need electrostatics back one here
        for fragment generation, or maybe templates or a new force field.
    4. More 3D synthon analysis