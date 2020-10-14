import IMP
import IMP.core
import IMP.atom
import IMP.threading
import IMP.algebra
import IMP.pmi
import IMP.pmi.tools
import random
import numpy
import math
import operator
import os
from copy import deepcopy

class Loop():
    '''
    A Loop object stores a set of contiguous residues corresponding to an unmodeled length of sequence
    
    Loops also contain pointers to StructureElement object IDs immediately N-terminal and C-terminal to this loop.
    These IDs are key values for a StructureElement dictionary.

    If this loop starts at the N-terminus, then self.se_before_id = -1 and similarly is it is a C-terminal loop
    '''
    def __init__(self, chain_id="A", start=0, length=0, se_before_id=-1, se_after_id=-1):
        self.start = start
        self.chain_id = chain_id
        self.length = length
        self.se_before_id = se_before_id
        self.se_after_id = se_after_id
        
    def get_loop_residues(self):
        return list(range(self.start, self.start+self.length))

    def get_loop_chain(self):
        return self.chain_id

    def generate_loop_id(self):
        # The loop ID is a unique identifier that contains all of the information about this loop:
        # chain, start residue, length, SE before and SE after.
        return self.chain_id+"_"+str(self.start)+"_"+str(self.length)+"_"+str(self.se_before_id)+"_"+str(self.se_after_id)

    def is_res_in_loop(self, residue):
        # Check if the given residue [as a tuple: (chain_id, resnum)] is contained in this loop
        chain_id = residue[0]
        if resnum in self.get_loop_residues() and chain_id == self.chain_id:
            return True
        else:
            return False

class SSEThread(IMP.ModelObject):
    def __init__(self, sequences={}, max_res_dist=4.0, se_clash_dist=4.0):
        '''
        @param sequences :: dictionary of sequences {chain_id:FASTA} or list of tuples [(chain_id, FASTA)]
        @param max_res_dist :: The maximum distance-per-residue (in angstroms) for computing loop lengths
        @param se_clash_dist :: Minimum distance that two SEs can be to each other (in angstroms)
        '''
        self.model = IMP.Model()
        self._setup_SSE_system()
        self.max_res_dist = max_res_dist
        self.se_clash_dist = se_clash_dist

        # if sequences is a list, make it a dictionary
        if isinstance(sequences, list):
            self.sequences = {}
            for chain in sequences:
                self.sequences[chain[0]]=chain[1]

        if len(sequences.keys()) > 0:
            self.build_sequence_chains(sequences)
        
        self.system_xyzs = [IMP.core.XYZ(p) for p in IMP.atom.Selection(self.sequence_hierarchy).get_selected_particles()]
    
    def _setup_SSE_system(self):
        self.root_hier = IMP.atom.Hierarchy.setup_particle(IMP.Particle(self.model))

        self.structure_chain = IMP.atom.Chain.setup_particle(IMP.Particle(self.model), "X")
        self.sequence_hierarchy = IMP.atom.Hierarchy.setup_particle(IMP.Particle(self.model))

        self.root_hier.add_child(self.structure_chain)
        self.root_hier.add_child(self.sequence_hierarchy)

    def build_sequence_chains(self, sequences):
        # Given a list of (chain_id, fasta sequence) tuples, build the sequence chains

        for chain_id in sequences.keys():#range(len(sequences)):
            fasta = sequences[chain_id]
            new_seq_chain = IMP.atom.Chain.setup_particle(IMP.Particle(self.model), str(chain_id))
            self.sequence_hierarchy.add_child(new_seq_chain)

            sse_res_id = 1
            for i in range(len(fasta)):
                pr = IMP.Particle(self.model)
                # Setup residue/mass/xyzr
                res = IMP.atom.Residue.setup_particle(pr,
                                            IMP.pmi.tools.get_residue_type_from_one_letter_code(fasta[i]),
                                            sse_res_id)
                #res_particles.append(res.get_particle())
                IMP.atom.Mass.setup_particle(res.get_particle(), IMP.atom.get_mass(res.get_residue_type()))
                IMP.core.XYZR.setup_particle(res.get_particle())

                # Initialize coordinates to 0,0,0
                IMP.core.XYZ(pr).set_coordinates((0,0,0))

                # We use the "are optimized" flag to indicate structured or not
                # Here, we begin with an unstructured residue
                IMP.core.XYZ(pr).set_coordinates_are_optimized(False)

                new_seq_chain.add_child(res)
                sse_res_id+=1

        return self.sequence_hierarchy

    def get_sequence_chain_lengths(self):
        chain_lengths = {}
        
        for c in self.sequence_hierarchy.get_children():
            chain_id = IMP.atom.Chain(c).get_id()
            chain_lengths[chain_id] = len(c.get_children())

        return chain_lengths

    def extract_structure_elements_from_smotif_pdbs(self, smotif_pdbs, stride=None):
        '''
        Given a set of PDBs, create two structure elements for each
        '''
        m2 = IMP.Model()
        self.structure_elements = {}
        se_ix = 0
        for p in smotif_pdbs:
            h = IMP.atom.read_pdb(p, m2, IMP.atom.CAlphaPDBSelector())

            # First, get CA atoms
            cas = IMP.atom.Selection(h, atom_type=IMP.atom.AT_CA).get_selected_particles()

            coords0 = []
            coords1 = []
            first = True

            residm1 = 0

            # First number of residues is in filename
            # A_1_26_43_44_56_0.pdb
            fname = p.split("/")[-1]
            nres0 = int(fname.split("_")[3])-int(fname.split("_")[2]) + 1
            #nres1 = int(p.split("_")[5])-int(p.split("_")[4]) + 1

            # Collect coordinates
            for c in range(nres0):
                coords0.append(IMP.core.XYZ(cas[c]).get_coordinates())
            for d in range(c+1, len(cas)):
                coords1.append(IMP.core.XYZ(cas[d]).get_coordinates())

            if stride is not None:
                n1 = int(fname.split("_")[2])
                n2 = int(fname.split("_")[3])+1
                n3 = int(fname.split("_")[4])
                n4 = int(fname.split("_")[5])+1
                stride0 = [stride[r] for r in range(n1, n2)]
                stride1 = [stride[r] for r in range(n3, n4)]
                #print(stride0, stride1)
                self.structure_elements[se_ix] = self.setup_structure_element(coords0, "H", 0, 1, 0, stride0)
                self.structure_elements[se_ix+1] = self.setup_structure_element(coords1, "H", 0, 1, 0, stride1)
            else:
                self.structure_elements[se_ix] = self.setup_structure_element(coords0, "H", 0, 1, 0)
                self.structure_elements[se_ix+1] = self.setup_structure_element(coords1, "H", 0, 1, 0)

            se_ix+=2

        #print("Keys", self.structure_elements.keys(), self.structure_elements.get_chain())
        return self.structure_elements

    def setup_structure_element(self, coordinates, sec_struct, start_residue=0, polarity=1, offset=0, chain_id="A", stride=None):

        se_pi = IMP.Particle(self.model)
        se_hier = IMP.atom.Hierarchy.setup_particle(se_pi)
        self.structure_chain.add_child(se_hier)

        # Add coordinates to hierarchy
        for c in range(len(coordinates)):
            coord = coordinates[c]
            np = IMP.Particle(self.model)
            hp = IMP.atom.Hierarchy.setup_particle(np)
            xyz = IMP.core.XYZR.setup_particle(np)
            xyz.set_coordinates(coord)
            xyz.set_radius(1.0)                     # radius doesn't really matter atm.
            IMP.atom.Mass.setup_particle(np, 1.0)   # Mass doesn't matter either
            se_hier.add_child(hp)

            if stride is not None:
                IMP.atom.SecondaryStructureResidue.setup_particle(np, stride[c][0][0], stride[c][0][1], stride[c][0][2])

        se = IMP.threading.StructureElement.setup_particle(self.model, 
                                                        se_pi.get_index(),
                                                        start_residue,
                                                        polarity,
                                                        len(coordinates),
                                                        offset,
                                                        chain_id)

        if sec_struct=="H":
            IMP.atom.SecondaryStructureResidue.setup_particle(se_pi, 1, 0, 0)
        elif sec_struct=="S":
            IMP.atom.SecondaryStructureResidue.setup_particle(se_pi, 0, 1, 0)
        elif sec_struct=="C":
            IMP.atom.SecondaryStructureResidue.setup_particle(se_pi, 0, 0, 1)
        else:
            raise Exception("Secondary structure designation must be H, S or C.", sec_struct, "was provided")
        return se

    def construct_se_min_loop_table(self):
        # Now build the hash table
        sids = list(self.structure_elements.keys())
        self.se_min_loop_table = numpy.ones((len(sids), len(sids)))

        for sid0 in range(len(sids)):
            si0 = self.structure_elements[sids[sid0]]
            for sid1 in range(sid0, len(sids)):
                si1 = self.structure_elements[sids[sid1]]

                if self.are_structure_elements_clashing(si0, si1):
                    val = -1
                else:
                    val = self.get_minimum_spanning_residues(si0.get_coordinates()[-1], si1.get_coordinates()[0])
                    self.se_min_loop_table[sid0][sid1] = int(val)
                    
                    val = self.get_minimum_spanning_residues(si1.get_coordinates()[-1], si0.get_coordinates()[0])
                    self.se_min_loop_table[sid1][sid0] = int(val)

        return self.se_min_loop_table 

    def get_minimum_spanning_residues(self, xyz1, xyz2):
        return math.ceil(IMP.algebra.get_distance(xyz1, xyz2) / self.max_res_dist)

    def are_structure_elements_clashing(self, se0, se1):
        dists = []
        for c0 in se0.get_coordinates():
            for c1 in se1.get_coordinates():
                dist = IMP.algebra.get_distance(c0, c1)
                dists.append(dist)
                if dist < self.se_clash_dist:
                    return True
        return False

    def get_built_residues(self):
        # Implementing a lookup table would be faster.
        built_res = []

        for xyz in self.system_xyzs:
            if xyz.get_coordinates_are_optimized():
                built_res.append(xyz.get_particle())
        #for p in IMP.atom.Selection(self.sequence_hierarchy).get_selected_particles():
        #    if IMP.core.XYZ(p).get_coordinates_are_optimized():
        #        #res = IMP.atom.Residue(self.model, p.get_index())
        #        built_res.append(p)
        return built_res

    def get_built_structure_element_ids(self, sort=False):
        built_ids = []
        built_sr = []

        # Ensure that we have a start res list
        try:
            ln = len(self.start_res_list)
        except:
            self.get_start_res_list()

        # Use numpy.where instead
        for s in self.structure_elements.keys():
            sr = self.start_res_list[s]#int(self.structure_elements[s].get_start_res())
            if sr[1] != 0: 
                built_ids.append(s)

        if sort:
            return self.sort_seids(built_ids)
        else:
            return built_ids

    def get_start_res_list(self):
        # Return a list of start residues as (chain_id, resnum)
        
        self.start_res_list=[]
        for se in range(len(self.structure_elements.keys())):
            self.start_res_list.append( (self.structure_elements[se].get_chain(), 
                                        int(self.structure_elements[se].get_start_res()))
                                        )
        return self.start_res_list

    def get_all_loops(self):
        return self.get_loops(self.get_built_structure_element_ids())

    def sort_seids(self, seids):
        # Returns a list of SEIDs
        if seids is None or len(seids)==0:
            return []

        ses = []
        
        seq_chains = [IMP.atom.Chain(c).get_id() for c in self.sequence_hierarchy.get_children()]

        sorted_seids = {}
        for c in seq_chains:
            chain_seids = []
            for s in seids:
                se = self.structure_elements[s]
                if se.get_chain() == c:
                    chain_seids.append((s, se.get_start_res()))

            s_chain_seids = sorted(chain_seids, key=operator.itemgetter(0))
            sorted_seids[c] = [x[0] for x in s_chain_seids]

        return sorted_seids

    def get_loops(self, seids, set_as=True):
        # Given a set of StructureElement IDs, return a dictionary of Loop objects

        # First, get a dictionary of SEIDs sorted by chain and residue number
        sorted_chain_ses_ids = self.sort_seids(seids)

        loops = {}

        for chain_id in sorted_chain_ses_ids:
            sorted_ses_ids = sorted_chain_ses_ids[chain_id]
          
            # The first loop in a chain always starts at 1 and has no SE N-terminal to it
            new_loop_start = 1
            new_se_before_id = -1
          
            for s in range(len(sorted_ses_ids)):
              
                seid = sorted_ses_ids[s]
                start_res = int(self.structure_elements[seid].get_start_res())
                start = new_loop_start
                length = start_res - new_loop_start - 1
                new_loop = Loop(chain_id=chain_id, 
                                start=start,
                                length=length,
                                se_before_id=new_se_before_id,
                                se_after_id=seid)
                loop_id = new_loop.generate_loop_id()
                
                #Generate loop ID and add to dictionary
                loops[loop_id] = new_loop

                # New loop starts at the end of this StructureElement
                new_loop_start = start_res + self.structure_elements[seid].get_length()                

            # Now do the last loop
            length = self.get_sequence_chain_lengths()[chain_id]-new_loop_start
            
            new_loop = Loop(start=new_loop_start,
                        length=length,
                        se_before_id=new_se_before_id,
                        se_after_id=-1)
            loop_id = new_loop.generate_loop_id()
            loops[loop_id] = new_loop
            new_se_before_id=-1

        # Set this loop dictionary as the system's
        if set_as:
            self.loops=loops

        return loops

    def is_se_allowed_in_model(self, seid):
        built_seids = self.get_built_structure_element_ids()

        # If structure element is built, then it is allowed
        if seid in built_seids:
            return True

        # If any built SE clashes (table=-1), then it is not allowed
        for b in built_seids:
            if self.se_min_loop_table[seid][b]==-1:
                return False
        return True

    def get_available_start_residues(self, seid, same_loop=False):
        # Given a structure element, find the start residues available to it
        available_residues = []

        se = self.structure_elements[seid]
        # If not allowed due to clash, then return nothing
        if not self.is_se_allowed_in_model(seid):
            return []

        # If allowed, then look at all loops
        se_len = int(se.get_length())

        # If same loop, we want to only return residues that are within the loop vacated by this SE
        if same_loop:
            before_loop, after_loop = self.get_loop_ids_bracketing_seid(seid)
            all_loop_resis = list(range(self.loops[before_loop].start, self.loops[after_loop].start + self.loops[after_loop].length))
            br = 0
            ar = 0
            if self.loops[before_loop].se_before_id != -1:
                br = int(self.se_min_loop_table[self.loops[before_loop].se_before_id][seid])
            if self.loops[after_loop].se_after_id != -1:
                ar = int(self.se_min_loop_table[seid][self.loops[after_loop].se_after_id])
            
            #print("Same loop", structure_element_id, se.get_start_res(), "||", before_loop, after_loop, "||", len(all_loop_resis), ar, se_len, all_loop_resis[br:len(all_loop_resis)-ar-se_len])
            if len(all_loop_resis) - ar - se_len - br + 1 <=0:
                print("No resis.  Weird!!")
                print("Same loop", seid, se.get_start_res(), se_len, "||", before_loop, after_loop, "||", len(all_loop_resis), br, ar, all_loop_resis[br:len(all_loop_resis)-ar-se_len+1])
                exit()
                return []
            else:
                for resnum in all_loop_resis[br:len(all_loop_resis)-ar-se_len+1]:
                    available_residues.append((se.get_chain(), resnum))
                return available_residues       
        else:
            loops = self.loops

        for lk in loops.keys():

            loop = loops[lk]

            if loop.length < 0:
                raise Exception("Loop length is less than zero", lk, len(loops.keys()))
            
            # If loop is not large enough to hold SE, then continue
            if loop.length < se_len:
                continue

            # if loop starts at N-terminus, then there is no loop before offset
            if loop.se_before_id==-1:
                br=0
            # Otherwise, see how many residues this SE needs to be from the SE before
            else:
                br = int(self.se_min_loop_table[loop.se_before_id][seid])
            # Now do the same for the C-terminus
            if loop.se_after_id==-1:
                ar=0
            else:
                ar = int(self.se_min_loop_table[seid][loop.se_after_id])

            # If, after subtracting the offsets, the loop is now too small, continue without adding residues
            if loop.length - br - ar - se_len <= 0:
                continue

            else:
                resis = loop.get_loop_residues()

                for resnum in resis[br:len(resis)-ar-se_len-1]:
                    available_residues.append(loop.chain_id, resnum)
        
        return available_residues

    def get_loop_ids_bracketing_seid(self, seid):
        # From the loops dictionary, return the loop before and loop after the given SEID
        # If there is no loop before or after (i.e., it's at a terminus), return None for that element

        before_loop_id = None
        after_loop_id = None

        for lk in self.loops.keys():
            loop = self.loops[lk]
            if loop.se_before_id == seid:
                after_loop_id = lk
            if loop.se_after_id == seid:
                before_loop_id = lk
        return before_loop_id, after_loop_id

    def compare_seids_and_loop_table(self):
        # Self-consistency check that loops and structure elements are the same
        # Currently used as a debugger and prints out all errors.

        # We check things until we find an inconsistency.
        good = True 

        # First, check the before and after SEIDs for the loops
        for lk in self.loops.keys():
            l = self.loops[lk]

            # If not an N-terminal loop
            if l.se_before_id!=-1:
                se = self.structure_elements[l.se_before_id]
                
                # Check that chains are the same
                if se.get_chain() != l.chain_id:
                    print("Loop and SE have different chains!!", seid, lk)
                    good=False

                # Check that the loop starts one residue after the N-terminal SE
                if se.get_resindex_list()[-1]+1 != l.start:
                    print(se.get_resindex_list(), se.get_start_res(), se.get_length())
                    print("--XX-- SEID", l.se_before_id, se.get_start_res(), se.get_resindex_list()[-1], ":: LOOP after", l.start, l.start+l.length, l.get_loop_residues)
                    good=False

            # If not a C_terminal loop
            if l.se_after_id!=-1:
                se = self.structure_elements[l.se_after_id]

                # Check that chains are the same
                if se.get_chain() != l.chain_id:
                    print("Loop and SE have different chains!!", seid, lk)
                    good=False

                # Check that the C-terminal SE starts one residue after the loop ends
                if se.get_start_res() != l.start + l.length:
                    print(se.get_resindex_list(), se.get_start_res(), se.get_length())
                    print("--XX-- LOOP", l.start, l.start+l.length, l.get_loop_residues(),"SEID After", l.se_after_id, se.get_start_res(), se.get_resindex_list()[-1],)
                    good=False
        
        # TODO:: Can't currently do this second check.  Need to get all residues.
        '''
        # Second, check that there are no overlapping SEs
        all_residues = list(range(1,self.get_sequence_lengths()[l.chain_id]))

        for i in self.get_built_structure_element_ids():
            se = self.structure_elements[i]
            for r in se.get_resindex_list():
                if r not in all_residues:
                    print("RESID", r, "in multiple SEs")
                    for ix in self.get_built_structure_element_ids():
                        if r in self.structure_elements[ix].get_resindex_list():
                            print("  ", ix, self.structure_elements[ix].get_resindex_list())
                    good=False
                else:
                    all_residues.remove(r)
        '''

        # Third, check that each SEID is only used once as a beginning and end
        before_ids = []
        after_ids = []
        for lk in self.loops:
            l = self.loops[lk]
            if l.se_before_id != -1:
              if l.se_before_id in before_ids:
                print("****", l.se_before_id, "SEID is before two loops")
                print("SEID resis:", self.structure_elements[l.se_before_id].get_resindex_list())
                print("This loop:", l.get_loop_residues())
                for lkx in self.loops:
                    if self.loops[lkx].se_before_id==l.se_before_id:
                        print("Other loop:", self.loops[lkx].get_loop_residues())
                good=False
              else:
                before_ids.append(l.se_before_id)
            if l.se_after_id !=-1:
              if l.se_after_id in after_ids:
                print("****", l.se_after_id, "SEID is after two loops")
                print("SEID resis:", self.structure_elements[l.se_after_id].get_resindex_list())
                print("This loop:", l.get_loop_residues())
                for lkx in self.loops:
                    if self.loops[lkx].se_after_id==l.se_after_id:
                        print("Other loop:", self.loops[lkx].get_loop_residues())
                good=False
              else:
                after_ids.append(l.se_after_id)
        
        return good

    def remove_se_from_loop_table(self, seid):
        # Combine the loops immediately N- and C-terminal of this SE

        # First, get the two loops that we will combine        
        before_loop_id, after_loop_id = self.get_loop_ids_bracketing_seid(seid)
      
        chain_id = self.sequence_elements[seid].get_chain()
        # If no loops are found, this SEID was not built into the model
        # Just return.
        if before_loop_id is None and after_loop_id is None:
            return

        # Compute parameters for the new loop
        if before_loop_id is not None:
            before_loop = self.loops[before_loop_id]
            start = before_loop.start
            se_before_id = before_loop.se_before_id
        else:
            # This SE starts at the N-terminus of a chain
            # Start of the new loop will now be the start residue of the SEID
            # And there is no SE before this loop now (since it starts at N-terminus
            after_loop = self.loops[after_loop_id]
            start = 1
            se_before_id = -1

        if after_loop_id is not None:
            after_loop = self.loops[after_loop_id]
            length = after_loop.start + after_loop.length - start
            se_after_id = after_loop.se_after_id
        else:
            # This SE ends at the C-terminus.
            length = self.get_sequence_chain_lengths()[chain_id] - start
            se_after_id = -1
       
        if length < 0:
            print("rem HELP!", length, start)
        new_loop = Loop(chain_id=chain_id,
                        start=start,
                        length=length,
                        se_before_id=se_before_id,
                        se_after_id=se_after_id)

        new_loop_id = new_loop.generate_loop_id()

        # Remove the two loops we connected
        if before_loop_id is not None:
            del self.loops[before_loop_id]
        if after_loop_id is not None:
            del self.loops[after_loop_id]
        
        #print("Deleted", seid, self.structure_elements[seid].get_start_res(), self.structure_elements[seid].get_length(),
        #        "| del loops:", before_loop_id, after_loop_id, "| Add loops:", new_loop_id)

        # And add our new loop to the dictionary
        self.loops[new_loop_id] = new_loop

    def add_se_to_loop_table(self, seid, new_sr):
        # Given an SEID and a new start residue, modify the loop table to account for this.

        chain_id = new_sr[0]
        resnum = new_sr[1]

        # ensure no overlap:
        # TODO - reimplement for multi-chain. Add chain_id to "get_built_residues".
        '''
        built = [IMP.atom.Residue(r).get_index() for r in self.get_built_residues(chain_id=chain_id)] 
        for r in range(new_sr, new_sr+self.structure_elements[seid].get_length()):
            if r in built:
                print(built)
                print("What you doing!?!?!", r, "is already built", self.structure_elements[seid].get_resindex_list())
                exit()
        '''
        # If the new SR is zero, do not add anything
        if new_sr == 0:
            return

        this_loop = None
        
        # Find the loop containing new_sr
        for lk in self.loops.keys():
            if self.loops[lk].is_res_in_loop(chain_id, resnum):
                this_loop = self.loops[lk]
                break
        
        if this_loop is None:
            self.compare_seids_and_loop_table()
            raise Exception("No loop found for SEID", seid, "in chain", chain_id, "starting at", resnum, [(self.loops[s].start, self.loops[s].start+self.loops[s].length-1) for s in self.loops.keys()])
        
        extra_string=""
        # Mode is the number of new loops we make
        mode = None
        # if the new start res is the start of this loop, then we only create one new loop
        '''
        if new_sr == this_loop.start:
            start = new_sr+self.structure_elements[seid].get_length()
            new_loop = Loop(start=start,
                            length=this_loop.start+this_loop.length-start,
                            se_before_id=seid,
                            se_after_id=this_loop.se_after_id)
            self.loops[new_loop.generate_loop_id()] = new_loop
            mode = 0
        
        # If the new SE goes all the way to the end of this loop, then again, only create one loop
        elif new_sr == this_loop.start+this_loop.length-self.structure_elements[seid].get_length():
            start = this_loop.start
            length = new_sr-start
            se_before_id = this_loop.se_before_id
            se_after_id = seid
            new_loop = Loop(start=start,
                           length=length,
                           se_before_id=se_before_id,
                           se_after_id=se_after_id)
            self.loops[new_loop.generate_loop_id()] = new_loop
            se_before_id = this_loop.se_before_id
            mode = 1
        
        # Otherwise, we need to create two loops out of one.
        '''
        if 0==0:
            # Make the first loop
            start = this_loop.start
            length = resnum - start

            if length < 0:
                print("HELP!", length, start, this_loop)

            se_before_id = this_loop.se_before_id
            se_after_id = seid
            new_loop = Loop(start=start,
                           length=length,
                           se_before_id=se_before_id,
                           se_after_id=se_after_id)
            self.loops[new_loop.generate_loop_id()] = new_loop
            extra_string = new_loop.generate_loop_id()

            # Second loop starts from the end of the SEID
            start = new_sr + self.structure_elements[seid].get_length() 
            
            # Length is the end of the previous loop - start residue.
            length = this_loop.start+this_loop.length-start
            
            if length < 0:
                print("HELP!", length, this_loop.start, this_loop.length, "|", chain_id, start, resnum, self.structure_elements[seid].get_length(), "|", this_loop.generate_loop_id(), extra_string)
            
            se_before_id = seid
            se_after_id = this_loop.se_after_id

            new_loop = Loop(chain_id=chain_id,
                           start=start,
                           length=length,
                           se_before_id=se_before_id,
                           se_after_id=se_after_id)
            self.loops[new_loop.generate_loop_id()] = new_loop
            mode = 2
        
        if mode == 2:
            snd_loop = self.loops[extra_string]
        
        else:
            i=0

        del self.loops[lk]


    def modify_loop_table(self, seid, new_sr):
        # Given a new start residue (as a (chain_id, resnum)), modify the loop table for this seid

        # First, remove the seid loop 
        self.remove_se_from_loop_table(seid)
        
        # Second, add the new loops to the table
        self.add_se_to_loop_table(seid, new_sr)

    def create_se_connectivity_restraints(self, structure_elements=None, function=None, n_sds=2):
        # Make N-1 connectivity restraints.

        if structure_elements is None:
            structure_elements = self.structure_elements
        # All restraints will not be used at all times.  
        # Any SECR restraints not used will be assigned dummy particles
        # which indicates the restraint should be evaluated as zero

        # make dummy particles for when we want to zero out SECR restraints
        self.secr_dummy_p0 = IMP.Particle(self.model)
        self.secr_dummy_p1 = IMP.Particle(self.model)
        
        if function is None:
            function = IMP.core.HarmonicUpperBound(0, 0.2) 

        self.secrs = []
        for i in range(len(structure_elements)-1):
            secr = IMP.threading.StructureElementConnectivityRestraint(self.model,
                                            function,
                                            self.secr_dummy_p0.get_index(), 
                                            self.secr_dummy_p1.get_index(),
                                            n_sds)

            self.secrs.append(secr)
        return self.secrs

    def update_SECR_restraints(self):
        try:
            n_secr = len(self.secrs)
        #if there are no secr restraints, then just exit
        except:
            return

        secr_id = 0
    
        # Find the connected SEs. Easiest through loop table
        for lk in self.loops.keys():
            if self.loops[lk].se_before_id == -1 or self.loops[lk].se_after_id == -1:
                continue
            else:
                # Those loops with SEs before and after indicate a connection
                se0 = self.structure_elements[self.loops[lk].se_before_id]
                se1 = self.structure_elements[self.loops[lk].se_after_id]

                # Assign these seids to the restraint
                self.secrs[secr_id].assign_particles(se0.get_particle_index(), se1.get_particle_index())

                secr_id+=1
        
        # Place the dummy particles in the remaining SECRs
        for i in range(secr_id, n_secr):
            self.secrs[secr_id].assign_particles(self.secr_dummy_p0.get_index(), self.secr_dummy_p1.get_index())

    def update_system(self):
        # Place in here all items that should be updated each step. 
        self.update_SECR_restraints()

    def add_SS_to_sequence(self, psipred, chain):
        self.sequence_ss = {}
        for res in psipred.keys():
            sys_part = IMP.atom.Selection(self.sequence_hierarchy, residue_index=res, chain_id=chain).get_selected_particle_indexes()[0]
            IMP.atom.SecondaryStructureResidue.setup_particle(self.model, sys_part, psipred[res][0], psipred[res][1], psipred[res][2])
         
            self.sequence_ss[res] = psipred[res]

    def add_SS_to_ses(self, secstruct):
        for seid in secstruct.keys():
            se_res = self.structure_elements[seid].get_children()
            for res in range(len(se_res)):
                ss = secstruct[seid][res]
                IMP.atom.SecondaryStructureResidue.setup_particle(self.model, se_res[res].get_particle_index(), 
                        ss[0], ss[1], ss[2])

    def print_model(self):
        # Print a string with the loops and structure elements by sequence:

        #First, order the loops by SR:

        sorted_loop_keys = []
        loop_ids = self.loops.keys()
        lkd = {}
        for lk in loop_ids:
            sr = int(lk.split("_")[0])
            lkd[sr] = lk

        for i in sorted(lkd.keys()):
            sorted_loop_keys.append(lkd[i])

        outstring = ""
        last_after_id = -1
        for lk in sorted_loop_keys:
            outstring+=":"+str(self.loops[lk].start)+"-loop-"+str(self.loops[lk].start+self.loops[lk].length-1)+":"
            if self.loops[lk].se_after_id != -1:
                se = self.structure_elements[self.loops[lk].se_after_id]
                outstring+=":"+str(se.get_start_res())+"##STR##"+str(se.get_start_res()+se.get_length()-1)+":"
        
        print(outstring)

    def get_all_model_ss_propensities(self, unbuilt_as_coil=True):
        # Returns a dictionary of tuples
        self.model_ss = {}
        built_ses = self.get_built_structure_element_ids()
        all_resids = set(self.sequence_ss.keys())
        for seid in built_ses:
            se = self.structure_elements[seid]
            se_res = IMP.atom.Hierarchy(se).get_children()
            resids = se.get_resindex_list()
            # Remove these built resids
            for r in resids:
                if r not in all_resids:
                    print("Residue", r, " not in residue table")
                    for seid in built_ses:
                        se = self.structure_elements[seid]
                        print(seid, "||", se.get_resindex_list())
                    raise Exception("Residue", r, " not in residue table")
                all_resids.remove(r)

            for r in resids:
                self.model_ss[r] = IMP.atom.SecondaryStructureResidue(self.model, 
                                    se_res[r-se.get_start_res()].get_particle_index()).get_all_probabilities()
        # All unbuilt residues are modeled as coil or None
        for r in all_resids:
            if unbuilt_as_coil:
                self.model_ss[r]=(0,0,1)
            else:
                self.model_ss[r]=(0,0,0)

        return self.model_ss

class CompletenessRestraint(IMP.Restraint):
    '''
    A restraint on the number of residues in the model that are built

    Evaluated as the number of unbuilt residues time a constant (slope)

    Unbuilt residues are determined by X = 0.0, as is used in SSEThread
    '''
    def __init__(self, root_hierarchy, slope=1.0):
        self.model = root_hierarchy.get_model()
        self.xyzs = [IMP.core.XYZ(p) for p in IMP.atom.Selection(root_hierarchy, resolution=1).get_selected_particles()]
        self.slope = slope
        IMP.Restraint.__init__(self, self.model, "MockRestraint %1%")

    def get_number_of_built_residues(self):
        built=0
        for xyz in self.xyzs:
            if xyz.get_coordinate(0) != 0.0:
                built+=1

        return built

    def unprotected_evaluate(self, da):
        return ((len(self.xyzs)-self.get_number_of_built_residues())*self.slope)

    def do_get_inputs(self):
        return []

    def do_show(self, fh):
        fh.write('CompletenessRestraint')


class SecondaryStructureParsimonyRestraint(IMP.Restraint):
    def __init__(self, system, psipred, evaluate_unbuilt=True):
        # psipred is a list of tuples: (chain_id, psipred dictionary)
        # psipred_dictionary is formatted as: { resnum:(frac_helix, frac_sheet, frac_coil), ... }
        self.system = system

        for c in psipred.keys():
            pp = psipred[c]
            self.system.add_SS_to_sequence(pp, c)
        self.evaluate_unbuilt=evaluate_unbuilt
        IMP.Restraint.__init__(self, self.system.model, "SecondaryStructureParsimonyRestraint")

    def unprotected_evaluate(self, da):
        score = 0
        # First compute structured residues, as determined by built structure elements
        model_ss = self.system.get_all_model_ss_propensities(unbuilt_as_coil=self.evaluate_unbuilt)
        
        for res in self.system.sequence_ss:
            resscore=0
            for i in range(3):
                resscore+=self.system.sequence_ss[res][i]*model_ss[res][i]
            if resscore < 0.1:
                resscore = 0.1
            score+=-1*math.log(resscore)

        return score

    def do_get_inputs(self):
        return []
    
    def do_show(self, fh):
        fh.write('SecondaryStructureParsimonyRestraint')



class SSEAdditionMover(IMP.threading.SSEThreadMover):
    def __init__(self, system, structure_elements):
        '''
        A move that proposes deletion of a random SE from the model
        @param system :: the SSEThread object
        @param structure_elements :: A dictionary of structure elements with key - SEID
        '''
        se_pis = [structure_elements[seid].get_particle_index() for seid in structure_elements.keys()]
        
        IMP.threading.SSEThreadMover.__init__(self, 
                system.model, 
                se_pis,
                system.sequence_hierarchy.get_particle_index())
        
        self.structure_elements = structure_elements
        self.system = system

    def get_se(self, seid):
        return self.structure_elements[self.seid]

    def get_random_structure_element(self):
        seids = self.structure_elements.keys()
        ri = numpy.randint(0,len(seids))
        self.structure_elements[seids[ri]]

    def do_propose(self, new_start_res=None, seid=None):
        '''
        Propose a new start residue for a StructureElement.

        If no start residue (new_start_res) is given, choose one at random from the system

        If no StructureElement (seid) is given, choose one at random

        '''
        # Get a random structure element that is not built
        if seid is None:
            rand_srs = list(self.structure_elements.keys())
            random.shuffle(rand_srs)

            for seid in rand_srs:
                se = self.structure_elements[seid]
                sr = se.get_start_res()
                if sr == 0:
                    break

                # If none are zero, we can't add one.  Return nothing.
                if seid == rand_srs[-1]:
                    return IMP.core.MonteCarloMoverResult([se.get_particle_index()], 0.0)

        self.seid = seid # Store the seid in a persistent variable
        se = self.structure_elements[self.seid]

        # Pick a random available start residue and change the SE start_res to that
        if new_start_res is None:
            available_start_res = self.system.get_available_start_residues(self.seid)
            if len(available_start_res)==0:
                self.mod = False
                return IMP.core.MonteCarloMoverResult([se.get_particle_index()], 0.0)

            new_start_res = random.choice(available_start_res)
        
        self.mod = True
        # Modify the loops table and start res list for the system
        self.system.add_se_to_loop_table(self.seid, new_start_res)
        self.system.start_res_list[self.seid] = new_start_res
        
        # Change the start residue and chain keys on the SE
        se.set_start_res_key(new_start_res[1])
        se.set_chain_key(new_start_res[0])
        self.new_start_res = new_start_res

        # Update the model
        se_pix = self.structure_elements[self.seid].get_particle_index()
        self.transform_coordinates(se_pix)
        
        self.system.update_system()
        self.old_start_res = 0
        return IMP.core.MonteCarloMoverResult([se.get_particle_index()], 1.0)

    def do_reject(self):
        
        if not self.mod:
            return
        se_pix = self.structure_elements[self.seid].get_particle_index()
        
        # Zero out the coordinates we added
        self.zero_coordinates(se_pix)
        
        # Reset the loops table and SR list
        self.system.start_res_list[self.seid]=0
        self.system.remove_se_from_loop_table(self.seid)
        
        # Change the start_res for this SEID
        self.structure_elements[self.seid].set_start_res_key(0)
        
        self.system.update_system()
        
class SSEDeletionMover(IMP.threading.SSEThreadMover):
    '''
    A move that proposes deletion of a random SE from the model
    '''
    def __init__(self, system, structure_elements):
        se_pis = [structure_elements[seid].get_particle_index() for seid in structure_elements.keys()]
        IMP.threading.SSEThreadMover.__init__(self, 
                system.model, 
                se_pis,
                system.sequence_hierarchy.get_particle_index())
        self.structure_elements = structure_elements
        self.system = system

    def get_se(self, seid):
        return self.structure_elements[seid]

    def do_propose(self):
        
        # Get a random structure element that is built
        rand_srs = list(self.structure_elements.keys())
        random.shuffle(rand_srs)

        for seid in rand_srs:
            se = self.structure_elements[seid]
            sr = se.get_start_res()

            # Find the first SE that is not built
            if sr != 0:
                break

            # If none are unbuilt, we can't add one.  Retnurn nothing.
            if seid == rand_srs[-1]:
                return IMP.core.MonteCarloMoverResult([se.get_particle_index()], 0.0)

        self.seid = seid
        self.old_start_res = se.get_start_res()
        self.old_chain = se.get_chain() # Don't need to use this since we don't change the chain
        se_pix = se.get_particle_index()
        
        # Change sequence coordinates to zero and set start_res_key of SE to zero
        self.zero_coordinates(se_pix)
        se.set_start_res_key(0)
        self.new_start_res = 0

        self.system.start_res_list[self.seid] = (self.old_chain, 0)

        # Modify the loops table for the system
        self.system.remove_se_from_loop_table(self.seid)
        self.system.update_system()

        return IMP.core.MonteCarloMoverResult([se.get_particle_index()], 1.0)

    def do_reject(self):
        se = self.get_se(self.seid)
        se_pix = se.get_particle_index() 

        # Reset the loops table, reset start_Res_key and transform the coordinates back.
        self.system.add_se_to_loop_table(self.seid, (self.old_chain, self.old_start_res))
        self.system.start_res_list[self.seid] = (self.old_chain, self.old_start_res)
        se.set_start_res_key(self.old_start_res)
        
        # transform the coordinates back to the model that we deleted
        self.transform_coordinates(se_pix)
        self.system.update_system()

class SSEShiftMover(IMP.threading.SSEThreadMover):
    '''
    A mover that shifts a random SSE start residue within the same loop
    '''
    def __init__(self, system, structure_elements, max_disp=None):
        se_pis = [structure_elements[seid].get_particle_index() for seid in structure_elements.keys()]
        
        IMP.threading.SSEThreadMover.__init__(self,
                system.model,
                se_pis,
                system.sequence_hierarchy.get_particle_index())
        
        self.structure_elements = structure_elements
        self.system = system
        self.max_disp = max_disp

    def get_se(self, seid):
        return self.structure_elements[seid]

    def set_max_disp(self, max_disp):
        self.max_disp = max_disp

    def do_propose(self):
        # Find random built SE
        rand_srs = list(self.system.get_built_structure_element_ids())
        random.shuffle(rand_srs)
        
        for seid in rand_srs:
            se = self.structure_elements[seid]

            # Find a built SE
            if sr != 0:
                break

            # If none are built, we can't shift.  Retnurn nothing.
            if seid == rand_srs[-1]:
                self.seid=-1
                return IMP.core.MonteCarloMoverResult([se.get_particle_index()], 0.0)
        
        self.seid = seid
        self.old_start_res = se.get_start_res()
        self.old_chain = se.get_chain()
        se_pix = se.get_particle_index()
        
        # find start residues within this loop (if true)
        available_start_res = self.system.get_available_start_residues(self.seid, same_loop=self.same_loop)

        lasr = len(available_start_res)

        # Parse available start res by those only within N residues in the same chain
        if self.max_disp is not None:
            available_start_res = [i for i in available_start_res if (sr[0]==self.old_chain and abs(i-sr[1])<=self.max_disp)]
       
        if len(available_start_res) == 0:
            print("ShiftMover WARNING: No available residues returned for this loop. Suspicious.", lasr)
            return IMP.core.MonteCarloMoverResult([se.get_particle_index()], 0.0)
        
        # Choose one of these at random
        new_start_res = random.choice(available_start_res)
        self.new_start_res = new_start_res

        # Modify the loops table for the system
        self.system.remove_se_from_loop_table(self.seid)
        
        # Zero the coordinates
        self.zero_coordinates(se_pix)
        
        # Add new loop to table
        self.system.add_se_to_loop_table(self.seid, new_start_res)
        
        # Transform coordinates and update tables
        se.set_start_res_key(new_start_res[1])
        se.set_chain_key(new_start_res[0])
        self.transform_coordinates(se_pix)
        self.system.start_res_list[self.seid] = new_start_res
        self.system.update_system()

        return IMP.core.MonteCarloMoverResult([se.get_particle_index()], 1.0)

    def do_reject(self):
        if self.seid == -1:
            return
        se = self.structure_elements[self.seid]
        se_pix = se.get_particle_index()
        new_sr = se.get_start_res()
        
        # Reset the loops table
        self.system.remove_se_from_loop_table(self.seid)
        self.system.start_res_list[self.seid] = (self.old_chain, self.old_start_res)
        
        # If new start res is not zero, then zero out the coordiantes we added
        if new_sr != 0:
            self.zero_coordinates(se_pix)
        
        self.system.add_se_to_loop_table(self.seid, (self.old_chain, self.old_start_res))

        # Set back the old new start residue
        se.set_start_res_key(self.old_start_res)
        se.set_chain_key(self.old_chain)
        
        # Transform the old coordinates back
        if self.old_start_res != 0:
            self.transform_coordinates(se_pix)
        
        self.system.update_system()

class SEMover(IMP.threading.StructureElementMover):
    # Mover wrapping an IMP.threading.StructureElementMover
    # This controls a single StructureElement
    def __init__(self, system, seid, zero_pct=50):
        IMP.threading.StructureElementMover.__init__(self, 
                system.model, 
                system.structure_elements[seid].get_particle_index(),
                system.sequence_hierarchy.get_particle_index())
        self.seid = seid
        self.system = system
        if zero_pct > 0:
            if zero_pct < 1.0:
                self.zero_pct = zero_pct*100
            elif zero_pct < 100:
                self.zero_pct = zero_pct
            else:
                raise Exception("zero_pct must be between zero and 100:", zero_pct)
        else:
            raise Exception("zero_pct must be between zero and 100:", zero_pct)

    def get_se(self):
        return self.system.structure_elements[self.seid]

    def do_propose(self):
        se = self.get_se()
        self.old_start_res = se.get_start_res()

        if self.old_start_res != 0:
            self.zero_coordinates()
        
        # if SE is not built, try to put it in:
        if self.old_start_res==0:
            start_res = self.system.get_available_start_residues(self.seid)
            if len(start_res)==0:
                return IMP.core.MonteCarloMoverResult([se.get_particle_index()], 0.0)

        # Otherwise, try to move it to another residue OR set it to zero 
        else:
            start_res = self.system.get_available_start_residues(self.seid, exclude_self=True)
            start_res = start_res+[0]*int(len(start_res)*self.zero_pct/100.0) # simple hack to add zeroes

        if len(start_res) == 0:
            return IMP.core.MonteCarloMoverResult([se.get_particle_index()], 0.0)
        
        new_start_res = random.choice(start_res)
        se.set_start_res_key(new_start_res)
        if new_start_res != 0:
            self.transform_coordinates()
        
        self.system.start_res_list[self.seid]=new_start_res

        # Modify the loops table for the system
        self.system.modify_loop_table(self.seid, new_start_res)

        return IMP.core.MonteCarloMoverResult([se.get_particle_index()], 1.0)


    def do_reject(self):
        new_sr = self.get_se().get_start_res()
        
        # Reset the loops table
        self.system.modify_loop_table(self.seid, self.old_start_res)
        
        if new_sr != 0:
            self.zero_coordinates()
        
        self.get_se().set_start_res_key(self.old_start_res)
        
        if self.old_start_res != 0:
            self.transform_coordinates()

def compare_threading_models(pdb_queries, pdb_reference, close_ca_threshold=5.0, exact_match_threshold=6.0):
    # Given two PDB files, return a number of metrics regarding the fit

    # Completeness:  How many residues built in query vs. reference

    # Close:  Number of query CA atoms within 5 angstroms of any reference CA atom. 

    # Sequence:  Number of query CA atoms within 5 angstroms of exact reference CA atom.

    if not hasattr(pdb_queries, "__iter__"):
        pdb_queries = pdb_queries

    m = IMP.Model()

    h = IMP.atom.read_pdb(pdb_reference, m, IMP.atom.CAlphaPDBSelector())

    ref_parts = IMP.atom.Selection(h, resolution=1).get_selected_particles()

    coords = {}
    for rp in ref_parts:
        resid = IMP.atom.Residue(IMP.atom.Hierarchy(rp).get_parent()).get_index()
        coords[resid] = IMP.core.XYZ(rp).get_coordinates()
        
    for pdb in pdb_queries:
        h2 = IMP.atom.read_pdb(pdb, m, IMP.atom.CAlphaPDBSelector())

        query_parts = IMP.atom.Selection(h2, resolution=1).get_selected_particles()
        
        query_coords = {}
        
        diffs = {}
        # Log coordinates and compute differences
        n_seq_close = 0
        n_ca_close = 0
        for pi in query_parts:
            resid = IMP.atom.Residue(IMP.atom.Hierarchy(pi).get_parent()).get_index()
            crds = IMP.core.XYZ(pi).get_coordinates()
            query_coords[resid] = crds
            try: 
                diffs[resid] = IMP.algebra.get_distance(crds, coords[resid]) 
                print(resid, diffs[resid])
                if diffs[resid] < exact_match_threshold:
                    n_seq_close+=1

                for ck in coords.keys():
                    if IMP.algebra.get_distance(crds, coords[ck]) < close_ca_threshold:
                        n_ca_close+=1
                        break
            except:
                continue

        print(len(coords.keys()), len(query_coords.keys()), n_seq_close, n_ca_close)

class MonteCarlo():
    def __init__(self, system, scoring_function, temp=1.0):
        self.scoring_function = scoring_function
        self.system = system
        self.movers = []
        self.current_score = None
        self.temp = temp
        self.reset_metrics()
        self.num_moves_per_step = 10
        self.addition_pct = 40
        self.deletion_pct = 100
        self.shift_pct = 100

    def set_num_moves_per_step(self, num):
        self.num_moves_per_step = num

    def get_number_of_accepted_steps(self):
        return self.number_of_accepted_steps

    def get_number_of_proposed_steps(self):
        return self.number_of_proposed_steps

    def get_number_of_rejected_steps(self):
        return self.number_of_rejected_steps

    def get_number_of_downward_steps(self):
        return self.number_of_downward_steps

    def get_number_of_upward_steps(self):
        return self.number_of_upward_steps

    def get_last_accepted_energy(self):
        return self.current_score

    def set_temp(self, temp):
        self.temp = temp

    def reset_metrics(self):
        self.number_of_downward_steps = 0
        self.number_of_upward_steps = 0
        self.number_of_accepted_steps = 0
        self.number_of_rejected_steps = 0
        self.number_of_proposed_steps = 0

    def add_movers(self, movers):
        if not hasattr(movers, "__iter__"):
            movers = [movers]
        self.movers += movers

    def get_movers(self):
        return self.movers

    def move(self, move_set):
        prob = 1.0
        #print("MOVE SET:", move_set)
        
        # First, move everything
        for m in move_set:
            #print("  ----------------------------------------", m)
            #self.system.print_model()
            #print([(seid, self.system.start_res_list[seid]) for seid in self.system.get_built_structure_element_ids()])
            mcr = self.movers[m].do_propose()
            if not self.system.compare_seids_and_loop_table():
                print("PROBLEM!")
                exit()
            prob*=mcr.get_proposal_ratio()
        return prob

    def get_proposal_move_set(self, num_moves):
        # Hack here to tune step sequences
        # Assume mover 0 is addition, mover 1 is deletion and mover 2 is shift.
        
        if num_moves==None:
            num_moves = self.num_moves_per_step

        # Want additions and deletions based on number of built residues
        
        prop_additions = (1.0*(len(self.system.system_xyzs)-len(self.system.get_built_residues()))/len(self.system.system_xyzs))*self.addition_pct
        prop_deletions = 1.0*len(self.system.get_built_structure_element_ids())*self.deletion_pct/100.0
        prop_shift = 1.0*len(self.system.get_built_structure_element_ids())*self.shift_pct/100.0
       
        #print("PROPS", prop_additions, prop_deletions, prop_shift)
        moves = []
        for i in range(num_moves):
            rnum = numpy.random.rand()*(prop_additions+prop_deletions+prop_shift)
            if rnum < prop_additions:
                moves.append(0)
            elif rnum < prop_additions + prop_deletions:
                moves.append(1)
            else:
                moves.append(2)
        #print("MOVES:", moves, rnum, "|", prop_additions, prop_deletions, prop_shift)
        self.move_set = moves
        return self.move_set

    def do_one_step(self, num_moves=None):
        move_set = self.get_proposal_move_set(num_moves)
        #start_res = deepcopy(self.system.get_start_res_list()) 
        prob = self.move(move_set)
        new_score = self.scoring_function.evaluate(False)
        accept = self.metropolis(new_score, prob)
        #new_start_res = self.system.get_start_res_list()
        
        if accept:
            self.number_of_accepted_steps+=1
            self.current_score = new_score
        else:
            self.number_of_rejected_steps+=1
            self.reject_moves(move_set)
            #rej_start_res = self.system.get_start_res_list()
            #for i in range(len(start_res)):
            #    if start_res[i] != rej_start_res[i]:
            #        print("Move set:", move_set, prob, accept, self.current_score, new_score)
            #        raise Exception("Rejected mvoe did not revert to original model. SE,", i, start_res[i],rej_start_res[i]) 
        return accept, new_score, move_set

    def metropolis(self, new_score, proposal_ratio):
        # If there is no score, make it higher than new_score
        # so we always accept this step
        if self.current_score is None:
            self.current_score = new_score + 1
        
        self.number_of_proposed_steps+=1

        diff = new_score - self.current_score
        if diff < 0:
            self.number_of_downward_steps+=1
            accept = True
        else:
            e = math.exp(-1*diff/self.temp)
            r = numpy.random.rand()
            if e * proposal_ratio > r:
                self.number_of_upward_steps+=1
                accept = True
            else:
                accept = False
        return accept
    
    def reject_moves(self, move_set):

        for m in move_set:
            self.movers[m].do_reject()


def read_config_file(config_file, config_dict):
    # Given a configuration file, returns a dictionary of parameters
    config_path = os.path.realpath(config_file)
    
    f = open(config_file, "r")
    
    for line in f.readlines():    
        # Burn any comment lines or blank lines
        if line.strip()=="" or line[0]=="#":
            continue

        fields = line.split()

        # Add sequences to system
        if fields[0] == "sequence":
            if len(fields[1])!=1:
                raise Exception("Chain ID must be a single letter:", fields[1])
                exit()

            if "sequence" not in config_dict.keys:
                config_dict["sequence"]=[]
            config_dict["sequence"].append((fields[1], fields[2]))

        elif fields[0] == "localization_table":
            path = config_path+"/"+fields[1]
            # should check if file exists
            if not os.path.isfile(path):
                raise Exception("Localization table file", path, "does not exist")
            config_dict["localization_table"] = path

        elif fields[0] == "pdb_file":
            path = config_path+"/"+fields[1]
            # should check if file exists
            if not os.path.isfile(path):
                raise Exception("PDB file", path, "does not exist")

            if "pdb_files" not in config_dict.keys():
                config_dict["pdb_files"]=[]
            config_dict["pdb_files"].append(path)

        elif fields[0] == "crosslink_dataset":
            path = config_path+"/"+fields[1]
            if not os.path.isfile(path):
                raise Exception("Crosslink file", path, "does not exist")

            if "crosslink_dataset" not in config_dict.keys():
                config_dict["crosslink_dataset"]=[]
            config_dict["crosslink_dataset"].append((path, float(fields[2]), float(fields[3]), float(fields[4])))

        else:
            # Everything else is a single parameter/pair
            # convert numbers to floats (will need to int all of the integer parameters in the main code)
            try:
                val = float(fields[1])
            except:
                val = fields[1]
            config_dict[fields[0]] = val

        return config_dict

def is_config_dict_complete(config_dict):
    # Does this configuration dictionary have all of the required elements to run a simulation?

    # Need sequences of some sort
    if "sequence" not in config_dict.keys():
        return False
    if len(config_dict["sequence"]) == 0:
        return False

    if "pdb_file" not in config_dict.keys() or "localization_table" not in config_dict.keys():
        return False

    return True

usage = '''
SSEThread command-line usage:

python SSEThread.py mode [args]

Available modes: setup sample anaysis

# Setup:
python SSEThread.py setup -c config_file -o output_folder

# Running
python SSEThread.py sample -c config_file [-b]
-b : benchmark - run 1000 steps to time a sampling run

# Analysis
python SSEThread.py analysis -o output_folder [-t truth.pdb]
'''

default_configuration_file = os.path.dirname(os.path.realpath(__file__))+"/default_configuration.dat"

# Scripts for setting up and running the system
if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-c", "--config", dest="config_file", type="string", help="Path to configuration file")
    parser.add_option("-b", "--benchmark", action="store_true", dest="benchmark", default=False, help="Run short benchmark for timing Only used with 'sampling' mode")
    parser.add_option("-o", "--output", dest="output_folder", type="string", help="Directory to put/read output for this system.")
    parser.add_option("-t", "--truth", dest="truth", type="string", default=None, help="Reference PDB file for comparing SSEThread results.")
    (options, args) = parser.parse_args(sys.argv[2:])
    
    config_dict = read_configuration_file(default_configuration_file, config_dict = {})

    mode = sys.argv[1].lower()

    if mode == "setup":
        configuration = read_config_file(options.config_file, config_dict=config_dict)
        setup_ssethread(config_dict=configuration, output_directory=options.output_folder)

    elif mode == "sampling":
        configuration = read_config_file(options.config_file, config_dict=config_dict)
        sample_ssethread(config_dict=configuration, output_directory=options.output_folder, benchmark=options.benchmark)

    elif mode == "analysis":
        analyze_ssethread(output_directory=options.output_folder, truth=options.truth)

    else:
        print("Unknown mode")
        print("")
        print(usage)
    
