# Omni Tool

Program to generate structures in lammps data format

Supported structures:

- Sphere
- TennisBall
- SphereJanus
- Globular_Sphere
- SpherePatch
- OblateSpheroid

- Icosahedron
- Dodecahedron
- Pentamer

- Chain
- Monomer


./omni_tool gen_glob
Example file to generate globular_sphere:
$ cat gen_glob
Particle_type: globular_sphere
Output_type: lammps_full # other keywords: xyz pdb lammps_full
Num_of_beads: 50
Scale: 5
Number_of_ligands: 200
Mol_tag: 1
Populate: 40 random
Sim_box: 0 50 0 50 0 50
ff_lj: 1 1.0 2.0 3.0 # type epsilon, sigma, cutoff
ff_cos2: 1 1.0 2.2 1.0 # type epsilon start_distance range
ff_lj: 2 1.0 1 3.0
ff_cos2: 2 1.0 0.7 1.0
ff_lj: 3 1.0 1 3.0
ff_cos2: 3 1.0 0.7 1.0
Patch: 0 0 1 0 0 0.75 2 # equation of plane: a*(x-cx) + b*(y-cy) ... > 0; a b c cx cy cz type_lig
Patch: 0 0 -1 0 0 -0.75 3
end

