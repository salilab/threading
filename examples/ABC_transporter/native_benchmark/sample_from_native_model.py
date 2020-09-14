import IMP
import IMP.core
import IMP.atom
import IMP.threading
import IMP.threading.SSEThread as SSEThread
import RMF
import IMP.rmf
import glob
import os
import sys
import time
import cProfile
import math

###########################################
#  Parameters
##########################################


outputdir = "native_sampling"

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

# Adaptive temperature sampling params
tmax = 100             # Maximum MC temperature
tmin = 1               # Minimum MC temperature
T_steps = 100          # Number of temperature steps
annealing_steps = 200  # Number of temperature cycles to run
steps_per_anneal = 500 # Number of steps per cycle
low_accept_range = 0.1 # If MC acceptance rate is below this, raise the temp
hi_accept_range = 0.30 # If MC acceptance rate is above this, lower the temp

# Equilibrium sampling params
t_equil = "last"         # Equilibrium sampling temperature. "last" uses last T from adaptive
n_equil_steps = 10000    # Number of equilibrium sampling steps after adaptive T
writefreq = 100          # Steps between writes in equilibration sampling

#######################################
# Functions and stuff
#######################################


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
                resindex1 = int(fields[2])+system.chain_breaks[chains[chain1]-1]
            
            if chains[chain2]==0:
                resindex2 = int(fields[3])
            else:
                resindex2 = int(fields[3])+system.chain_breaks[chains[chain2]-1]
            
            p1 = IMP.atom.Selection(system.sequence_hierarchy, residue_index=resindex1).get_selected_particles()
            p2 = IMP.atom.Selection(system.sequence_hierarchy, residue_index=resindex2).get_selected_particles()

            if perfect:
                xl_dist = float(fields[1])
            else:
                xl_dist = xl_dist

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

    return start_residues


sequenceA = "MTEDTYSKAFDRALFARILRYVWPYRLQVVLALLFLLVVTLAAAATPLFFKWAIDLALVPTEPRPLAERFHLLLWISLGFLAVRAVHFAATYGETYLIQWVGQRVLFDLRSDLFAKLMRLHPGFYDRNPVGRLMTRVTSDVDAINQFITGGLVGVIADLFTLVGLLGFMLFLSPKLTLVVLLVAPVLLAVTTWVRLGMRSAYREMRLRLARVNAALQENLSGVETIQLFVKEREREEKFDRLNRDLFRAWVEIIRWFALFFPVVGFLGDFAVASLVYYGGGEVVRGAVSLGLLVAFVDYTRQLFQPLQDLSDKFNLFQGAMASAERIFGVLDTEEELKDPEDPTPIRGFRGEVEFRDVWLAYTPKGVEPTEKDWVLKGVSFRVRPGEKVALVGATGAGKTSVVSLIARFYDPQRGCVFLDGVDVRRYRQEELRRHVGIVLQEPFLFSGTVLDNLRLFDPSVPPERVEEVARFLGAHEFILRLPKGYQTVLGERGAGLSTGEKQLLALVRALLASPDILLILDEATASVDSETEKRLQEALYKAMEGRTSLIIAHRLSTIRHVDRILVFRKGRLVEEGSHEELLAKGGYYAALYRLQFQEAKL"
sequenceB = "MTGRSAAPLLRRLWPYVGRYRWRYLWAVLAGLVSIFFFVLTPYFLRLAVDAVQAGRGFGVYALAIVASAALSGLLSYAMRRLAVVASRQVEYDLRRDLLHHLLTLDRDFYHKHRVGDLMNRLNTDLSAVREMVGPGILMGSRLSFLVLLAFLSMYAVNARLAFYLTLILPGIFLAMRFLLRLIDRRYREAQEVFDRISTLAQEAFSGIRVVKGYALERRMVAWFQDLNRLYVEKSLALARVEGPLHALLGFLMGFAFLTVLWAGGAMVVRGELSVGELVQFNAYLAQLTWPILGLGWVMALYQRGLTSLRRLFELLDEKPAIRDEDPLPLALEDLSGEVRFEGVGLKRDGRWLLRGLTLTIPEGMTLGITGRTGSGKSLLAALVPRLLDPSEGRVYVGGHEARRIPLAVLRKAVGVAPQEPFLFSETILENIAFGLDEVDRERVEWAARLAGIHEEILAFPKGYETVLGERGITLSGGQRQRVALARALAKRPKILILDDALSAVDAETEARILQGLKTVLGKQTTLLISHRTAALRHADWIIVLDGGRIVEEGTHESLLQAGGLYAEMDRLQKEVEA"


datadir = "../data/"


###############################
# Here is where the real work begins
###############################


os.makedirs(outputdir, exist_ok=True)

sequences = [sequenceA, sequenceB]
system = SSEThread.SSEThread(sequences)
stride_dict = open_minimal_stride_file(datadir+"6rag.stride", len(sequenceA)+len(sequenceB))
start_residues = build_smotifs_from_native(system, datadir+"6rag_CA.pdb", stride_dict)

loops = system.get_all_loops()
hash_table = system.construct_se_min_loop_table()

all_resis = set(list(range(1,len(sequenceA)+len(sequenceB))))

# Build in the movers
seam = SSEThread.SSEAdditionMover(system, system.structure_elements)
sedm = SSEThread.SSEDeletionMover(system, system.structure_elements)
sesm = SSEThread.SSEShiftMover(system, system.structure_elements)
# Set the shift maximum displacement
sesm.set_max_disp(shift_max_displacement)

# Build the restraints
psipred_dict = open_ss2_file(datadir+"6rag.ss2")

secr = system.create_se_connectivity_restraints()
sspr = SSEThread.SecondaryStructureParsimonyRestraint(system, psipred_dict)
cprs = import_crosslinks(system,"potential_crosslinks.txt", perfect=use_perfect, xl_dist=xl_dist, upper_threshold=xl_upper_dist, low_threshold=xl_lower_dist, random_choice=frac_xl_to_use)
comp = SSEThread.CompletenessRestraint(system.sequence_hierarchy, slope=completeness_slope)

restraints = secr + cprs+ [comp] + [sspr]
sf = IMP.core.RestraintsScoringFunction(restraints)


# Initialize MonteCarlo Object and add movers
mc = SSEThread.MonteCarlo(system, sf)
mc.add_movers([seam, sedm, sesm])
mc.set_num_moves_per_step(1)


# Initialize Output File
fh = open(outputdir+"/output.dat", "w")
outstring = "Frame Temp | "
for r in restraints:
    outstring+=" "+r.get_name()
fh.write(outstring+"Total_Score || start_resis\n")

print("INIT:", sf.evaluate(False), comp.unprotected_evaluate(None), mc.get_number_of_accepted_steps(), mc.get_number_of_proposed_steps())

writefreq = 100

t0 = time.time()

temps = []

deltaT = 0.955

# Initialize the temperature index to zero (highest temp).  I guess we could start elsewhere
Tix = 0
for n in range(annealing_steps):
    Tm = tmax * (deltaT ** Tix)
    mc.set_temp(Tm)
    moves = []
    for i in range(steps_per_anneal):
        #print("#-------------------------")
        result, score, move_set = mc.do_one_step()
        #print("Score, etc...", score, result, system.get_built_structure_element_ids(), len(system.get_built_residues()))
        #print("Startreslist:", system.start_res_list)
        #print("Built resis", [IMP.atom.Residue(r).get_index() for r in system.get_built_residues()])
        #print("MOVES:", move_set, result, score, "SEID:", mc.movers[move_set[0]].seid, mc.movers[move_set[0]].old_start_res, mc.movers[move_set[0]].new_start_res)
        n_built = len(system.get_built_residues())
        n_se_built  = 0
        for seid in system.get_built_structure_element_ids():
            n_se_built += system.structure_elements[seid].get_length()
        
        if n_se_built != n_built:
            print("Uh Oh!  Inconsistency in our model!", n_se_built, n_built)
            print("Moves:", move_set, result, "SEID:", mc.movers[move_set[0]].seid, mc.movers[move_set[0]].old_start_res)
            print("Built_resis:", [IMP.atom.Residue(r).get_index() for r in system.get_built_residues()])
            print("Start_res_list:", system.start_res_list)
            built_seids = []
            srs = []
            for seid in system.structure_elements.keys():
                if system.structure_elements[seid].get_start_res() != 0:
                    built_seids.append(seid)
                    srs.append(system.structure_elements[seid].get_start_res())
            six = sorted(range(len(srs)), key=lambda k: srs[k])
            totres = 0
            for ix in six:
                se = system.structure_elements[built_seids[ix]]
                print(" :",se.get_start_res(), se.get_length(), built_seids[ix])
                totres += se.get_length()

            print(system.get_available_start_residues(mc.movers[move_set[0]].seid))

            print("SE, Built:", totres, n_se_built, n_built)
            print("Exiting.  Use the above info to try and diagnose the problem")
            exit()
        moves+=move_set
    ts = 0

    # Yes, we're computing the restraints again. Need to get this from the restraints themselves.
    for r in restraints:
        sc = r.unprotected_evaluate(None)
        ts += sc
        outstring += " " + str(sc)
    outstring +=" " + str(ts)+" ||"
    for sr in system.get_start_res_list():
        outstring += " "+str(sr)
    fh.write(outstring)
    outstring += "\n"
    IMP.atom.write_pdb(system.sequence_hierarchy, outputdir+"/native_sa"+str(n)+".pdb")
    print("i, Temp:",n,  Tm, "| Score:", mc.get_last_accepted_energy(), score, "| Nres:", len(system.get_built_residues())) 
    print("Moves: %Add %Del %Shift", [1.0*moves.count(i)/steps_per_anneal for i in range(3)])
    print("MC Metrics, %accept #upward #downward:", 1.0*mc.get_number_of_accepted_steps()/steps_per_anneal, mc.get_number_of_upward_steps(), mc.get_number_of_downward_steps())
    print("----------------")
    if 1.0*mc.get_number_of_accepted_steps()/steps_per_anneal > hi_accept_range and Tix != T_steps:      
        Tix+=1
    elif 1.0*mc.get_number_of_accepted_steps()/steps_per_anneal < low_accept_range and Tix != 0:
        Tix-=1
    mc.reset_metrics()


if t_equil != "last":
    mc.set_temp(t_equil)
mc.reset_metrics()

for i in range(n_equil_steps):
    result = mc.do_one_step()
    
    # Write all scores to output
    outstring = str(i)+" "+str(mc.temp)+" | "
    ts = 0
    if i%writefreq == 0:
        for r in restraints:
            sc = r.unprotected_evaluate(None)
            ts += sc
            outstring += " " + str(sc)
        outstring +=" " + str(ts)+" ||"
        for sr in system.get_start_res_list():
            outstring += " "+str(sr)
        outstring += "\n"
        fh.write(outstring)

        print(len(system.get_built_residues()), ts, comp.unprotected_evaluate(None), mc.get_last_accepted_energy(), "|", mc.get_number_of_accepted_steps(), mc.get_number_of_proposed_steps(), mc.get_number_of_upward_steps(), mc.get_number_of_downward_steps()) 

        # Write frame to pdb file
        IMP.atom.write_pdb(system.sequence_hierarchy, outputdir+"/native_step"+str(i)+".pdb")
    
print("Sampler timing")
print("::", time.time()-t0, "for", n_steps, "steps")


