'''
Build the ABC transporter from Crystal-structure derived StructureElements
'''
import IMP
import IMP.core
import IMP.atom
import IMP.threading
import IMP.threading.SSEThread as SSEThread
import glob
import sys
import os
import time
import math

outputdir = "native_build"

# Mover Parameters
shift_max_displacement = 10  # Maximum displacement in residues for shift moves

# Scoring Parameters
completeness_slope = 5 # Score per-unmodeled residue. Higher values weight completing the sequence more
xl_lower_dist = 20     # Minimum distance for simulated crosslinks (file contains all Lys-Lys pairs)
xl_upper_dist = 35     # Maximum distance for simulated crosslinks
frac_xl_to_use = 100   # Percent of crosslinks within the distance thresholds to use
use_perfect= True      # Set the restraint distance equal to the model distance. If False, need XL_dist
xl_dist = 35           # Restraint distance (like your XL linker length + side chain. Required if use_perfect is False.
                       #   In theory, it should be higher than xl_upper_distance, but there is no check in place for that.

# Functions
chains = { "A" : 0,
           "B" : 1}

def open_minimal_stride_file(filename, nres):
    f= open(filename, "r")

    outdict = {}

    # Initialize everything to coil
    for i in range(1, nres+1):
        outdict[i] = (0.0,0.0,1.0)

    for line in f.readlines():
        fields=line.split(" ")
        ss = fields[1].strip()

        if ss=="H":
            ssp = (1.0,0.0,0.0)
        elif ss=="E":
            ssp = (0.0,1.0,0.0)
        else:
            ssp = (0.0,0.0,1.0)

        outdict[int(fields[0])]=ssp
    
    return outdict

def open_ss2_file(filename):
    f= open(filename, "r")

    outdict = {}
    for line in f.readlines():
        fields = line.split(" ")
        resid = int(fields[0])
        outdict[resid] = (float(fields[4]), float(fields[3]), float(fields[5].strip()))

    return outdict

def setup_conditional_pair_restraint(p1, p2, length, xl_slope, constant):
    r = IMP.threading.LoopPairDistanceRestraint(p1.get_model(), IMP.core.HarmonicUpperBound(length, xl_slope), 
        p1, p2, constant)
    return r

def import_crosslinks(system, infile, xl_dist, perfect=True, upper_threshold=35, low_threshold=0, random_choice=100, slope=0.2, constant=0.2):
    '''
    perfect = True means the restraint is built at the distance in the PDB
    otherwise, the restraint is built at the threshold distance (better simulating a crosslink set)
    AA 107.53908962663361 585 175
    '''
    rests = []
    f=open(infile, "r")
    for l in f.readlines():
        fields=l.split()
        if float(fields[1]) < upper_threshold and float(fields[1])>low_threshold:
            chain1 = fields[0][0]
            chain2 = fields[0][1]

            if chains[chain1]==0:
                resindex1 = int(fields[2])
            else:
                resindex1 = int(fields[2])+602#system.chain_breaks[chains[chain1]-1]
            
            if chains[chain2]==0:
                resindex2 = int(fields[3])
            else:
                resindex2 = int(fields[3])+602#+system.chain_breaks[chains[chain2]-1]
            
            p1 = IMP.atom.Selection(system.sequence_hierarchy, residue_index=resindex1).get_selected_particles()
            p2 = IMP.atom.Selection(system.sequence_hierarchy, residue_index=resindex2).get_selected_particles()

            if perfect:
                xl_dist = float(fields[1])

            cpr = setup_conditional_pair_restraint(p1[0], p2[0], xl_dist, slope, constant)
            rests.append(cpr)
    
    return rests

def build_smotifs_from_native(system, pdb, stride_dict):
    # Given an SSEThread system, PDB file and stride files
    # generate structure elements from the native PDB
    #
    # @param stride :: list of stride files in order
    system.structure_elements = {}
    start_residues = {}
    ss_segment = "X"
    ss_segment_residues = []
    structured_resids = set()
    m = IMP.Model()
    native_hier = IMP.atom.read_pdb(pdb, m, IMP.atom.CAlphaPDBSelector())
    seid = 0
    srs = []
    lens = []
    for res in range(1, len(stride_dict.keys())+1):
        ss_res = stride_dict[res]
        # If this STRIDE designation is not equal to the previous, we move on
        if ss_res != ss_segment:
            # If the length of the SS segment is greater than 4, 
            # then we make this into a structure element
            if len(ss_segment_residues) > 4:
                structured_resids.update(ss_segment_residues)
                if stride_dict[res-1] == (1,0,0):
                    sec_struct="H"
                elif stride_dict[res-1] == (0,1,0):
                    sec_struct="S"
                else:
                    sec_struct="C"

                if sec_struct != "C":
                    srs.append(ss_segment_residues[0])
                    lens.append(len(ss_segment_residues))
                    stride = [stride_dict[r] for r in ss_segment_residues]
                    native_parts = IMP.atom.Selection(native_hier, residue_indexes=ss_segment_residues, atom_type=IMP.atom.AT_CA).get_selected_particles()
                    native_coords = [IMP.core.XYZ(p).get_coordinates() for p in native_parts]
                    se = system.setup_structure_element(native_coords, sec_struct, start_residue=0, polarity=1, offset=0, stride=stride)
                    start_residues[seid] = ss_segment_residues[0]
                    system.structure_elements[seid] = se
                    seid+=1
            ss_segment = ss_res
            ss_segment_residues = [res]
         
        # If this is equal to the previous, we extend the ss_segment chain
        else:
            ss_segment_residues.append(res)
    outstring=""

    # Sometimes, nice to print out the start residues and lengths for these SEs
    for s in srs:
        outstring+=str(s)+", " 
    #print(outstring)
    
    outstring = ""
    for l in lens:
        outstring+=str(l)+", "
    #print(outstring)
    
    return start_residues


sequenceA = "MTEDTYSKAFDRALFARILRYVWPYRLQVVLALLFLLVVTLAAAATPLFFKWAIDLALVPTEPRPLAERFHLLLWISLGFLAVRAVHFAATYGETYLIQWVGQRVLFDLRSDLFAKLMRLHPGFYDRNPVGRLMTRVTSDVDAINQFITGGLVGVIADLFTLVGLLGFMLFLSPKLTLVVLLVAPVLLAVTTWVRLGMRSAYREMRLRLARVNAALQENLSGVETIQLFVKEREREEKFDRLNRDLFRAWVEIIRWFALFFPVVGFLGDFAVASLVYYGGGEVVRGAVSLGLLVAFVDYTRQLFQPLQDLSDKFNLFQGAMASAERIFGVLDTEEELKDPEDPTPIRGFRGEVEFRDVWLAYTPKGVEPTEKDWVLKGVSFRVRPGEKVALVGATGAGKTSVVSLIARFYDPQRGCVFLDGVDVRRYRQEELRRHVGIVLQEPFLFSGTVLDNLRLFDPSVPPERVEEVARFLGAHEFILRLPKGYQTVLGERGAGLSTGEKQLLALVRALLASPDILLILDEATASVDSETEKRLQEALYKAMEGRTSLIIAHRLSTIRHVDRILVFRKGRLVEEGSHEELLAKGGYYAALYRLQFQEAKL"
sequenceB = "MTGRSAAPLLRRLWPYVGRYRWRYLWAVLAGLVSIFFFVLTPYFLRLAVDAVQAGRGFGVYALAIVASAALSGLLSYAMRRLAVVASRQVEYDLRRDLLHHLLTLDRDFYHKHRVGDLMNRLNTDLSAVREMVGPGILMGSRLSFLVLLAFLSMYAVNARLAFYLTLILPGIFLAMRFLLRLIDRRYREAQEVFDRISTLAQEAFSGIRVVKGYALERRMVAWFQDLNRLYVEKSLALARVEGPLHALLGFLMGFAFLTVLWAGGAMVVRGELSVGELVQFNAYLAQLTWPILGLGWVMALYQRGLTSLRRLFELLDEKPAIRDEDPLPLALEDLSGEVRFEGVGLKRDGRWLLRGLTLTIPEGMTLGITGRTGSGKSLLAALVPRLLDPSEGRVYVGGHEARRIPLAVLRKAVGVAPQEPFLFSETILENIAFGLDEVDRERVEWAARLAGIHEEILAFPKGYETVLGERGITLSGGQRQRVALARALAKRPKILILDDALSAVDAETEARILQGLKTVLGKQTTLLISHRTAALRHADWIIVLDGGRIVEEGTHESLLQAGGLYAEMDRLQKEVEA"


datadir = "../data/"
outputdir = "./native_build_output/"

os.makedirs(outputdir, exist_ok=True)

sequences = [sequenceA, sequenceB]

# Setup system
system = SSEThread.SSEThread(sequences)
stride_dict = open_minimal_stride_file(datadir+"6rag.stride", len(sequenceA)+len(sequenceB))
start_residues = build_smotifs_from_native(system, datadir+"../data/6rag_CA.pdb", stride_dict)
psipred_dict = open_ss2_file(datadir+"6rag.ss2")

# We need an addition mover to add the SEs
seam = SSEThread.SSEAdditionMover(system, system.structure_elements)

srl=system.get_start_res_list()
loops = system.get_all_loops()
hash_table = system.construct_se_min_loop_table()

# Restraints
secr = system.create_se_connectivity_restraints()
completeness = SSEThread.CompletenessRestraint(system.sequence_hierarchy, completeness_slope)
sspr = SSEThread.SecondaryStructureParsimonyRestraint(system, psipred_dict, evaluate_unbuilt=True)
cprs = import_crosslinks(system,"potential_crosslinks.txt", perfect=use_perfect, xl_dist=xl_dist, upper_threshold=xl_upper_dist, low_threshold=xl_lower_dist, random_choice=frac_xl_to_use)

restraints = secr + cprs + [completeness] + [sspr]
sf = IMP.core.RestraintsScoringFunction(restraints)
fh = open(outputdir+"native_output.dat", "w")

outstring = "Frame:"
for r in restraints:
    outstring+=" "+r.get_name()

fh.write(outstring+"Total_Score\n")

print("SEID Total_Score Completeness SS Connectivity XL SE_start_res SE_length Total_built_residues")
for seid in system.structure_elements.keys():
    seam.do_propose(start_residues[seid], seid)
    
    outstring = str(seid)+" "
    ts = 0
    for r in restraints:
        sc = r.unprotected_evaluate(None)
        ts += sc
        outstring += " " + str(sc)
    outstring += " " + str(ts)
    fh.write(outstring +"\n")
    IMP.atom.write_pdb(system.sequence_hierarchy, outputdir+"native"+str(seid)+"_fit.pdb")
    secr_s = 0
    cprs_s = 0
    for s in secr:
        secr_s+=s.unprotected_evaluate(None)
    for s in cprs:
        cprs_s+=s.unprotected_evaluate(None)
    print(seid, ts, completeness.unprotected_evaluate(None), sspr.unprotected_evaluate(None), secr_s, cprs_s, system.structure_elements[seid].get_start_res(), system.structure_elements[seid].get_length(), len(system.get_built_residues()))



exit()
# Look at crosslinking restraint satisfaction:

for i in range(len(cprs)):
    r = cprs[i]
    pis = r.get_sequence_residue_particles()
    res1 = IMP.atom.Residue(system.model, pis[0])
    res2 = IMP.atom.Residue(system.model, pis[1])
    xyz1 = IMP.core.XYZ(system.model, pis[0])
    xyz2 = IMP.core.XYZ(system.model, pis[1])

    if xyz1.get_coordinates_are_optimized() and xyz2.get_coordinates_are_optimized():
        dist = IMP.core.get_distance(xyz1, xyz2)
        print(i, res1, res2, dist, xyz1.get_coordinates_are_optimized(), xyz2.get_coordinates_are_optimized(), r.unprotected_evaluate(None))
    else:
        print(i, res1, res2, xyz1.get_coordinates_are_optimized(), xyz2.get_coordinates_are_optimized(), r.unprotected_evaluate(None))
