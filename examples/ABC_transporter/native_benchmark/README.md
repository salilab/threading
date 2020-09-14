# Threading model of the ABC transporter using the native structure

The scripts in this folder create StructureElements from the PDB for
the ABC transporter (6RAG).  Default is to take every contiguous set
of sequence longer than four residues with the same SS designation
(from ../data/6rag.stride).  This results in 51 SEs covering 821 of 
the 1180 residues in the structure.

## Running script usage
`IMP_BUILD/setup_environment.sh python build_from_native_model.py`

This starts from an "empty" model and builds in one SE at each step
into its "correct" assignment.  


`IMP_BUILD/setup_environment.sh python sample_from_native_model.py`

This starts from an empty model and uses MonteCarlo sampling to place
SEs into the ABC transporter sequence. The current implementation
has an adaptive temperature function.

See the script header for a number of parameters that can be played with



