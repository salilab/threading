import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy


def get_best_scoring_models(inp, num_best=10, threshold=None):
    # Hard-code in thresholds
    #if num_best > len(inp):
    #    print "You have fewer models than number of best scoring models!"
    #    print "Returning you all of the models"

    #for sorted_list = sorted(inp, key=lambda x: (x[1]))

    return sorted_list[:num_best]

def plot_scores(inp):

    fig = plt.subplot()
    sorted_list = sorted(inp, key=lambda x: (x[1]))

    scores = [i[1] for i in sorted_list]

    fig.plot(scores)

    plt.show()

def calc_resis_from_keys(se_length, keys):
    # given a set of keys, return a list of se_length numbers
    # that correspond to the mapping og SE to sequence
    # start_res, polarity, length, offset

    outlist = []

    if keys[1]==1:
        for i in range(int(keys[3])):
            outlist.append(0)
        for i in range(int(keys[2])):
            outlist.append(int(keys[0] + i))
        for i in range(int(se_length-keys[2]-keys[3])):
            outlist.append(0)
    else:
        # Start with the end of the list
        for i in range(int(se_length-keys[2]-keys[3])):
            outlist.append(0) 
        # insert the residues at the begining (reverse order)
        for i in range(int(keys[2])):
            outlist.insert(int(keys[0] + i),0)
        # now put the offset Xs in the beginning
        for i in range(int(keys[3])):
            outlist.insert(0,0)
    return outlist 

def calc_se_residue_coverage(res_lists):
    '''
    Given a list of residue numbers for this structural element,
    return a list of the % of times that structural residue is assigned to a
    residue in sequence (is not zero)
    '''
    out_pcts = numpy.zeros(len(res_lists[0]))
    res_list_array = numpy.array(res_lists)

    len_array = res_list_array.shape[0]

    for i in range(len(out_pcts)):
        out_pcts[i]=len(res_list_array[:,i].nonzero()[0]) * 1.0 / len_array

    return out_pcts

def calculate_all_se_coverages(inp):
    coverage = []
    for s in range(len(inp[0][0])):
        se = []
        for i in inp:
            se.append(calc_resis_from_keys(se_lengths[s], i[0][s]))
        coverage.append(list(calc_se_residue_coverage(se)))

    return coverage


def calc_se_residue_histogram(res_lists, length):
    # given a list of list of residue numbers, for this SE,
    # return a numpy array with the frequency of each residue number

    histogram = numpy.zeros(length-1)

    for r in res_lists:
        #h = numpy.array(numpy.histogram(r, bins=range(1, length+1))[0])
        #print type(h), type(histogram), h.shape, histogram.shape
        #rint h, histogram
        #print h,  histogram
        histogram += numpy.histogram(r, bins=range(1, length+1))[0]

    return histogram

def calc_se_residue_2D_histogram(res_lists, seq_length, se_length):
    # given a list of list of residue numbers, for this SE,
    # return a numpy array with the frequency of each residue number

    histogram = numpy.zeros(length-1)

    for r in res_lists:
        #h = numpy.array(numpy.histogram(r, bins=range(1, length+1))[0])
        #print type(h), type(histogram), h.shape, histogram.shape
        #rint h, histogram
        #print h,  histogram
        histogram += numpy.histogram(r, bins=range(1, length+1))[0]

    return histogram

def count_2d_list(l):
    # Given a 2D list, return the total number of elements
    counter = 0
    for i in l:
        counter += len(i)

    return counter

def plot_all_se_coverages(covs):
    nres = count_2d_list(covs)
    fig, axes = plt.subplots(1,len(covs), figsize=(nres/2+(len(covs)-1)*0.1 + 1,3))
    fig.subplots_adjust(wspace=0.1)
    start=1
    axes[0].set_ylabel("Percent Occupied in BSMs")
    for i in range(len(covs)):
        axes[i].set_xlabel("Residue Number")

        x = numpy.arange(1,len(covs[i])+1)
        axes[i].set_xlim(1,len(covs[i]))
        axes[i].plot(x, numpy.array(covs[i]), color = 'blue')#numpy.array(covs[0]))
        axes[i].fill_between(x, 0, numpy.array(covs[i]), facecolor='blue')
        axes[i].set_xticks(numpy.arange(1,len(covs[i])+1,4))
        axes[i].tick_params(top='off')
        if i!=0:
            axes[i].tick_params(left='off', right='off')
            axes[i].set_yticklabels([])

    return fig

def get_list_of_seq_resis_by_timepoint(inp, se_lengths):
    all_tp_all_resis = []
    for i in inp:
        all_resis = []
        # Get the residues covered in this step over all SEs
        for s in range(len(se_lengths)):
            resis = calc_resis_from_keys(se_lengths[s], i[0][s])
            for r in resis:
                if r != 0:
                    all_resis.append(r)

        all_tp_all_resis.append(all_resis)

    return all_tp_all_resis

def plot_sequence_coverage(seq_length, models, se_lengths):
    seq_resis_by_tp = get_list_of_seq_resis_by_timepoint(models, se_lengths)

    seq_resis = [item for sublist in seq_resis_by_tp for item in sublist]

    hist = numpy.histogram(numpy.array(seq_resis), bins=numpy.arange(0.5, seq_length+1.5,1))

    fig, ax = plt.subplots(1,1, figsize=(seq_length*0.5,3))

    ax.set_ylim([0,1])

    ax.plot(range(1,71), hist[0]/(1.0*len(models)))

    return fig, ax

def calculate_threading_matrix_by_se(seq_length, models, se_lengths):
    # Get % of models with that residue each SE.
    n_se = len(se_lengths)
    thread_matrix = numpy.zeros((n_se, seq_length))

    for i in models:
        for s in range(n_se):
            resis = calc_resis_from_keys(se_lengths[s], i[0][s])
            for r in resis:
                if r != 0:
                    thread_matrix[s][r-1]+=1.0
    return thread_matrix / (1.0*len(inp))

def calculate_threading_matrix(seq_length, models, se_lengths):
    # Get % of models with that residue in SE sequence position
    n_se_res = numpy.sum(se_lengths)
    n_se = len(se_lengths)
    thread_matrix = numpy.zeros((n_se_res, seq_length))
    for i in models:
        offset=0
        for s in range(n_se):
            resis = calc_resis_from_keys(se_lengths[s], i[0][s])
            for r in range(se_lengths[s]):
                if resis[r] != 0:
                    thread_matrix[offset+r][resis[r]-1]+=1.0 / len(models)

            offset+=se_lengths[s]

    return thread_matrix

def plot_matrix(matrix, ylabel="Structure Element", xlabel="Residue Number"):
    plt.imshow(matrix, cmap="hot")
    plt.ytics([])

#print len(inp3)

rex_inp = []
import json
import glob
import matplotlib.patches as patches
import math

for outfile in glob.glob("enumerate_multi_*.dat"):
    with open(outfile) as f:
        lst = json.load(f)

        for m in lst['models']:
	    r = m['restraints']
            #rex_inp.append([m['model'], m['score']+(m['restraints']['SecondaryStructureParsimonyRestraint0']-1.73)*500,m['restraints']])
            if r['SecondaryStructureParsimonyRestraint0'] < 1.78 and r['SEConnectivityRestraint0'] + r['SEConnectivityRestraint1'] + r['SEConnectivityRestraint2'] + r['SEConnectivityRestraint3'] < 25 and not math.isnan(m['score']) and m['score'] < 75:
	    	rex_inp.append([m['model'], m['score'] , r])


print "Total parsed models", len(rex_inp)

#se_lengths = [10,25,14]
#se_lengths = [14,10,25]
se_lengths = [14,25,10]
#elements=[(20, 10, 'H'), (113, 25, 'H'), (152, 14, 'H')]
p_sites = [2603,2609,2612,2620,2624]#,2638,2645,2647,2649,2671,2672,2674,2675,2677,2743]



#list(map(lambda x: 100000.0 if math.isnan(x) else x, rex_inp))
#bsms = get_best_scoring_models([[i[0],100000,i[2]] if math.isnan(i[1]) else i for i in rex_inp], 10000)
bsms = rex_inp
f = open("bsms.dat", "w")

for bsm in bsms:
    #print bsm
    f.writelines(str(bsm[0][0][0])+" "+str(bsm[0][1][0])+" "+str(bsm[0][2][0])+" "+str(bsm[1])+" | " + ''.join(str(e)+" " for e in bsm[2].values()) + "\n")

f.close()


#print bsm[2].keys() 


#print bsms[2]
#print bsms[49]
#print bsms



#bsm = [[[[2677.0, 1.0, 14.0, 0.0], [2725.0, 1.0, 25.0, 0.0], [2754.0, 1.0, 10.0, 0.0]], 22.73919386256208]]

#for i in bsm[0][0]:
#    print i, seq[int(i[0])+1:int(i[0])+int(i[2])+1]


seq_length = 4128
#se_lengths = [15,11,12]
tm = calculate_threading_matrix(seq_length, bsms, se_lengths) / len(bsms)
#tm = calculate_threading_matrix(seq_length, bsms, se_lengths)

print numpy.sum(tm[0]), tm[0]

fig, ax=plt.subplots()
ax.imshow(tm, cmap="Greys", interpolation="none")
ax.set_yticklabels([])
ax.set_xlabel("Residue Number")
ax.set_xlim([2550,2800])
#ax.set_ylabel("Structure Element")
#elements=[(12,15,'H'), (34,11,'H'),(47,12,'H')]
'''
correct_elements = [[2677, 14, 1, 0], [2725, 25, 1, 0], [2754, 10, 1, 0]]

y = 0
# For matrix view
for e in correct_elements:
    for i in range(e[0],e[0]+e[1]):
        rect = patches.Rectangle((i-1.5, y-0.5),1,1,linewidth=0.35,edgecolor='r',facecolor='none')
        ax.add_patch(rect)
        print y, i, se_lengths
        y+=1
'''


plt.savefig('pkcs_matrix.png', dpi=900)

'''
for e in elements:
    rect = patches.Rectangle((e[0]-1.5, y-0.5),e[1], 1,linewidth=0.5,edgecolor='r',facecolor='none')
    ax.add_patch(rect)
    print y, se_lengths
    y+=1
plt.show()
'''

#plot_scores(pk_inp)
plt.show()



