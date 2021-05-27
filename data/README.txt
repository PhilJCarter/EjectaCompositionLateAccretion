PKDGRAV ss files:

The ss files (both initial conditions and outputs) contain properties for each 
particle at a particular time. These files can be read using the included python 
code. Contains particle mass (solar masses), particle radius (au), particle 
cartesian coordinates (au), particle cartesian velocity components 
(au/year/(2*pi)), particle ID number, and particle color/type.

The file header stores the time, the number of particles, and a file type flag.
                
--------------------------------------------------------------------------------

PKDGRAV origin_bins files:

The origin_bins files give the mass weighted provenance data for each particle 
in the corresponding ss output file. The format is as follows:

Number of particles (N)
N lines giving fraction of mass of each particle that originated in each of 10 
or 25 0.1 au wide bins
                
--------------------------------------------------------------------------------

PKDGRAV core_origin files:

The core_origin files have the same format as origin_bins files and give the 
provenance data for each particles core weighted by core mass.
                
--------------------------------------------------------------------------------

PKDGRAV cfrac files:

The cfrac files give the core mass fraction for each particle 
in the corresponding ss output file. The format is as follows:

Number of particles (N)
N lines giving core mass fraction of each particle
                
--------------------------------------------------------------------------------

PKDGRAV iord files:

The iord files give the particle ID number for each particle in the 
corresponding ss output file. The format is as follows:

Number of particles (N)
N lines giving ID for each particle
                
--------------------------------------------------------------------------------

PKDGRAV ss.par files:

The par files give the parameters used for the simulations.
                
