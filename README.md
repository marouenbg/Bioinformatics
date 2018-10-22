HeliStatistiX is a python v2.7 tool that allows to compute tilt angles of helices and packing angles of helices in interactions, mainly.

Contact: marouen.b.guebila@gmail.com

## Input

1. Put the .pdb files downloaded from OPM database in the same folder as the script of HeliStatistiX

2. Write a file named config containing 3 lines:

  - PDB codes of the proteins analyzed separated by ',' (comma)
	E.g., 1dxr,1e12,1eys,

  - Subunits of each protein separated by ',' and by ';' to separate proteins
	E.g., L,M;A;L,M;

  -Beginning and end of each helix (these information can be taken from OPM database). The start and end of each helix are separated by ','
Helices are separated by ':', Subunits are separated by '/' ,Proteins are separated by ';'
	E.g., 33,54:85,106:116,138:171,192:228,251/53,74:111,133:143,165:199,220:263,284;27,50:63,82:106,123:134,153:159,180:199,221:227,250;33,54:92,114:124,146:179,200:237,258/55,75:112,133:143,165:197,220:264,283;3,24:38,57:70,87:98,117:122,141:163,180:190,211/24,41:60,81;

## Ouput 

2 files :

-Tilt_angles.csv

-Packing_angles.csv

Tilt_angles.csv contains :

-The tilt angle,

-The helix length in aminoacids number,

-The helix length in Angstrom,

-The mean hydrophobicity of the helix,

-The hydrophobic moment of Eisenberg,

-The number of neighbor helices,

-The number of kinks,

-The greatest angle of kink (relative to the helix axis),

-The composition of aminoacids of the first tour 

-The composition of aminoacids of the last tour

Packing_angles.csv contains :

-The packing angle between a couple of helices in interaction,

-The helix orientation (A/P),

-Number of contacts between helices (based on Chothia definition),

-The minimum distance between the centers of each helix,

-The length of the conenxion loop,

-The number of kinks totalized by both helices

-The hydrophobic moment of Eisenberg of each helix


