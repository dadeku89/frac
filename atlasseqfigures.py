import os, sys, operator, subprocess, numpy, scipy, matplotlib
from numpy import *
from scipy.spatial.distance import pdist, squareform
#forces matplotlib to not use any Xwindows backend                                                
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import cophenet
from scipy.stats.stats import pearsonr
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as dist
from pylab import pcolor, colorbar, xticks, yticks
from scipy.cluster.hierarchy import fcluster
from pylab import *
from collections import Counter
import seaborn
import Mysbatch
from scipy.stats import spearmanr
sys.path.append('/home/dadeku89/scratch/lib_Daa/miso_events')
import miso_analysis
import kmerAnalysis
from miso_analysis import *
import pandas as pd
from os.path import basename
from scipy.stats import binom
from scipy.stats import ks_2samp
from matplotlib import font_manager as fm, rcParams
import matplotlib as mpl

#how I filtered the original tpm table so that at least one fraction has a tpm of zero for each gene 
def filterTpm(tpmTable,d1, d2,  outfile):
    outfile = open(outfile, 'w')
    fiOBltTpm = {}
    gene_sd = {}
    label = []
    d1 = int(d1)
    d2 = int(d2)
    for line in open(tpmTable):
        if line.startswith('#'):
            label = line.strip().split()[d1:d2]
        else:
            i = 0
            rowTot = []
            info = line.strip().split("\t")
            gene_name = info[0]
            tpm = info[d1:d2]
            tpm = map(float,tpm)
            for item in tpm:
                if item > 1:
                    i+=1
                    rowTot.append(i)
            if sum(rowTot) >= 1:
                filtTpm[gene_name] = tpm
    sortedGenes =  sorted(filtTpm.keys())
    print len(sortedGenes)
    outfile.write('#Gene' + '\t')
    for item in label:
        outfile.write(str(item) + '\t')
    outfile.write('\n')
    for item in sortedGenes:
        outfile.write(str(item) + '\t'+  '\t'.join(map(str, filtTpm[item])) + '\t' +  '\n')

#takes in the kallisto tpm table and makes an array of all genes in the analysis and the tpm count of each gene for each fraction with dataMatrix being the expression array and rowHeaders being the list of all genes matched
def tpm_array(tpmTable):
    with open(tpmTable, 'r') as f:
        rowHeaders = [] #this will contain all of the gene information for the tpm table            
        dataMatrix = [] #this will be the matrix of all of the rpkm data for each fraction         
        colHeaders = f.next().strip().split()[1::] #this contains all of the tpm values for each gene for each fraction                                            
        for line in f:
            data = line.strip().split('\t')
            rowHeaders.append(data[0])
            dataMatrix.append([float(x) for x in data[1::]])
    return dataMatrix, rowHeaders

def ms_array3(ms_table):
    msdata = []
    protlist = []
    with open(ms_table, 'r') as f:
        for line in f:
            i = 1 
            if not line.startswith('#'):
                info = line.strip().split('\t')
                desc = info[-8].split(' ')
                test = []
                for item in desc:
                    if item.split('=')[0] == 'GN':
                        test.append(item.split('=')[1])
                        if i == 1:
                            protlist.append(item.split('=')[1])
                        i +=1
                if i == 1:
                    if len(test) == 0:
                        protlist.append(desc[0])
                profile = line.strip().split('\t')[2:20]
                profile = [float(f) for f in profile]
                msdata.append(profile)
    
                
    return msdata, protlist
def genesincluster(clust_f):
    clustertogenes = {}
    for line in open(clust_f, 'r'):
        info = line.strip().split('\t')
        gene = info[0]
        cluster = int(info[1])
        if cluster not in clustertogenes:
            clustertogenes[cluster] = [gene]
        else:
            clustertogenes[cluster].append(gene)
    return clustertogenes

def clustGenerp(tpmTable, clust_f):
    clusterGene = {}
    data, genelist = tpm_array(tpmTable)
    data =  [[j/sum(row[0:17]) for j in row[0:17]] for row in data]
    clustmatchmsprof = {}
    for line in open(clust_f, 'r'):
        info = line.strip().split('\t')
        gene = info[0]
        cluster = int(info[1])
        if cluster not in clusterGene:
            clusterGene[cluster] = [gene]
        else:
            clusterGene[cluster].append(gene)
    for cluster in clusterGene:
        if len(clusterGene[cluster]) >= 20:
            meanclustlist = []
            for gene in clusterGene[cluster]:
                proflist = []
                index = genelist.index(gene)
                profile = data[index]
                prof1 = profile[0:17]
                prof1 = [f/sum(prof1) for f in prof1]
                meanclustlist.append(prof1)
            meanclustlist = array(meanclustlist)
            meanclustlist =  mean(meanclustlist, axis = 0)
            meanclustlist =meanclustlist.tolist()
            clustmatchmsprof[cluster] =  meanclustlist
    return clustmatchmsprof

def rnalistrbpmatch(tpmTable):
    data, genelist = tpm_array(tpmTable)
    data =  [[j/sum(row[0:17]) for j in row[0:17]] for row in data]
    return data, genelist

def rbplist(mannrbp_f, rbpdb_f, orthologue_f):
    totalrbplist = []
    humanrbp = []
    for line in open(mannrbp_f,'r'):
            if not line.startswith('Symbol'):
                humanrbp.append(line.strip())
                
    humantomouse = {}
    for line in open(orthologue_f, 'r'):
        info = line.strip().split()
        humangene = info[1]
        mousegene = info[3]
        humantomouse[humangene] = mousegene
    mannrbplist = []
    for rbp in humanrbp:
        if rbp in humantomouse:
            mannrbplist.append(humantomouse[rbp])
            totalrbplist.append(humantomouse[rbp])
    rpdblist = []
    for line in open(rbpdb_f,'r'):
        if not line.startswith('Symbol'):
            rbp = line.strip()
            rpdblist.append(rbp)
            totalrbplist.append(rbp)
    return totalrbplist

def generaterbpsinms(mannrbp_f, rbpdb_f, orthologue_f, ms_table):
    msdata, protlist = ms_array3(ms_table)
    totalrbplist = rbplist(mannrbp_f, rbpdb_f, orthologue_f)
    rbpinms = []
        
    for pep in totalrbplist:
        if pep in protlist:
            if pep not in rbpinms:
                rbpinms.append(pep)
    print len(rbpinms)
                
    return rbpinms    


#Figure 4 
def rpcorrelations(ms_table, tpmTable):
    msdata, protlist = ms_array3(ms_table)
    dataMatrix, rowHeaders = tpm_array(tpmTable)
    genelist = []
    notin = []
    rpvector = []
    rpcorrelations = []
    correlationlist = []
    genecorr = []
    for pep in protlist:
        if pep in rowHeaders:
            genelist.append(pep)
        else:
            notin.append(pep)
    for gene in genelist:
        pidx = protlist.index(gene)
        ridx = rowHeaders.index(gene)
        proteinprofile = msdata[pidx][0:17]
        proteinprofile = [f/sum(proteinprofile) for f in proteinprofile]
        rnaprofile = dataMatrix[ridx][0:17]
        rnaprofile = [f/sum(rnaprofile) for f in rnaprofile]
        rpvector.append([proteinprofile, rnaprofile])
        corr = corrcoef(proteinprofile, rnaprofile)[0][1]
        correlationlist.append(corr)
        genecorr.append([gene, corr])
    return genelist, rpvector, correlationlist, genecorr

def testrp(ms_table, tpmTable, prot):
    msdata, protlist = ms_array3(ms_table)
    dataMatrix, rowHeaders = tpm_array(tpmTable)
    for pep in protlist:
        if pep == prot:
            pidx = protlist.index(pep)
            proteinprofile = msdata[pidx][0:17]
            proteinprofile = [f/sum(proteinprofile) for f in proteinprofile]
            print 'Peptide:',pep
    for rna in rowHeaders:
        ridx = rowHeaders.index(rna)
        rnaprofile = dataMatrix[ridx][0:17]
        [f/sum(rnaprofile) for f in rnaprofile]
        corr = corrcoef(proteinprofile, rnaprofile)[0][1]
        
        if corr >0.85:
            #print 'Associated RNA',rna, corr
            print rna
    fraction = range(0,17)
    proteinprofile = [f*17 for f in proteinprofile]
    plt.plot(fraction, proteinprofile)
    #plt.ylim(0,3)
    plt.xlabel('Fraction')
    plt.ylabel('Relative Peptide Abundance')
    plt.savefig('test.pdf')
#Figure 3B
def corrhistRP(ms_table, tpmTable):
    genelist, rpvector, correlationlist, genecorr  = rpcorrelations(ms_table, tpmTable)
    bins =linspace(-1, 1, 15)
    plt.figure()
    ax = plt.subplot(1,1,1)
    
    plt.hist(correlationlist, bins, alpha = 0.5, align = 'mid', color = '#000000', rasterized = True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    print len(rpvector)
    plt.savefig('test.pdf')

def predictionTogene(signal_f):
    geneTopredict =  {}
    for line in open(signal_f, 'r'):
        if not line.startswith('#'):
            info = line.strip().split('\t')
            gene = info[1]
            prediction = info[4]
            geneTopredict[gene] = prediction

    return geneTopredict

#next ned to gt the list of genes that are associated with an er cluster and then determine among the genes that are associated with an er cluster how many of those genes are signal sequence containing for correlated vs anticorrelation, plot this result and compute the pvalue using the rank sums test to determine if the difference is statistically significant
def signalcorrrp(signal_f, ms_table, tpmTable):
    genelist, rpvector, correlationlist, genecorr  = rpcorrelations(ms_table, tpmTable)
    geneTopredict = predictionTogene(signal_f)
    plt.figure(figsize = (5,7))
    t = range(80,100,5)
    b = range(5,25,5)[::-1]
    for x in range(len(t)):
        i = 0
        j = 0
        m = 0 
        n = 0
        toppercent = percentile(correlationlist, t[x])
        bottompercent = percentile(correlationlist, b[x])
        for gene in genecorr:
            if gene[1]>= toppercent:
                if gene[0] in geneTopredict:
                    j +=1
                    if geneTopredict[gene[0]].split(' ')[1] == 'Signal':
                        i +=1
            if gene[1]<=bottompercent:
                if gene[0] in geneTopredict:
                    m +=1
                    if geneTopredict[gene[0]].split(' ')[1] == 'Signal':
                        n +=1
        topsignal =  float(i)/float(j)
        bottomsignal = float(n)/float(m)
        plt.bar(str(t[x]), topsignal, color = '#808080')#, edgecolor = '#000000')
        plt.bar(str(b[x]), bottomsignal,  color = '#454444')
    plt.xlabel('Bin')
    plt.ylabel('Percent Signal containing')
    plt.savefig('test.pdf')
def ergenes(er_f, clust_f, signal_f, ms_table, tpmTable):
    clustertogenes = genesincluster(clust_f)
    geneTopredict = predictionTogene(signal_f)
    genelist, rpvector, correlationlist, genecorr  = rpcorrelations(ms_table, tpmTable)
    sortedgenecorr = sorted(genecorr, key = itemgetter(1), reverse = True)
    erclusters = []
    for line in open(er_f, 'r'):
        cluster, organelle = line.strip().split('\t')
        if organelle == 'er':
            erclusters.append(int(cluster))
    ergenedict = {}
    for cluster in erclusters:
        for gene in clustertogenes[cluster]:
            
            if gene in geneTopredict:
                prediction = geneTopredict[gene]
                ergenedict[gene] = prediction
    for gene in sortedgenecorr:
        try:
            print ergenedict[gene[0]], gene[1]
        except:
            pass
            
def signalpcorrrp(ms_table, tpmTable, signal_f, clust_f):
    genelist, rpvector, correlationlist, genecorr  = rpcorrelations(ms_table, tpmTable)
    geneTopredict = predictionTogene(signal_f)
    clustertogenes = genesincluster(clust_f)
    ctcorr = []
    ctanticorr = []
    i = 0 
    j = 0
    c = 0
    a = 0
    for item in genecorr:
        if item[1]>0.5:
            gene = item[0]
            if gene in geneTopredict:
                #print item[0], geneTopredict[item[0]]
                i +=1
                print geneTopredict[item[0]].split(' ' )[1]
                if geneTopredict[item[0]].split(' ')[1] == 'Signal':
                    c+=1
        if item[1]<-0.825:
            gene = item[0]
            if gene in geneTopredict:
                #print item[0], geneTopredict[item[0]]
                j +=1
                if geneTopredict[item[0]].split(' ')[1] == 'Signal':
                    a+=1
    print i , j
    print c, a
            
        

def makegofilescorr(ms_table, tpmTable):
    genelist, rpvector, correlationlist, genecorr = rpcorrelations(ms_table, tpmTable)
    sortedcorr = sorted(genecorr, key = itemgetter(1), reverse = True)
    print 'Correlated:'
    for item in sortedcorr: 
        print item[1]
    sortedanticorr = sorted(genecorr, key = itemgetter(1)) 
    print 'Anticorrelated:' 
    for item in sortedanticorr:
        print item[1]
#figure 4c and figure 4d
def godotplot(go_f):
    compartmentlist = []
    enrichmentlist = []
    numberofgeneslist = []
    pvaluelist = []
    for line in open(go_f, 'r'):
        if not line.startswith("GO "):
            info = line.strip().split('\t')
            compartment = info[1]
            pvalue = float(info[2])
            fdrqvalue = info[3]
            enrichment = float(info[4])
            numberofgenes = float(info[6])
            compartmentlist.append(compartment)
            enrichmentlist.append(enrichment)
            pvaluelist.append(-log10(pvalue))
            numberofgeneslist.append(numberofgenes)
    cmap = matplotlib.cm.get_cmap('plasma')
    zipped = zip(compartmentlist, enrichmentlist, numberofgeneslist, pvaluelist)
    zipsort = sorted(zipped, key = operator.itemgetter(1), reverse = True)
    goanalysis = zipsort[:10]
    i = range(0,len(goanalysis))
    i = reversed(i)
    plt.figure(figsize =(2, 4))
    ax = plt.subplot(1,1,1)
    slist = sorted([x[2] for x in goanalysis], reverse=True)
    slistall =  slist[::4]
    for item in i:
        plt.scatter(zipsort[item][1], goanalysis[item][0], s = goanalysis[item][2], c = goanalysis[item][3], cmap = cmap, vmin = 3, vmax = 8)
        print zipsort[item][1], goanalysis[item][0], goanalysis[item][2]
    for item in reversed(range(0,len(slistall))):
        ax.scatter([], [], s = slistall[item], color = '#000000', label = slistall[item])
    maxenrich = max([x[1] for x in goanalysis])*.04
    cb = plt.colorbar(orientation = 'horizontal')
    cb.ax.tick_params(labelsize=10)
    plt.xticks(fontsize = 6)
    plt.yticks(fontsize = 6)
    plt.xlabel('Enrichment', fontsize = 6)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.xlim(1,3)
    plt.savefig('test.pdf', bbox_inches="tight")

#Figure 5A and Figure 5B





def filterrnas(tpmTable):
    filtrna = []
    for line in open(tpmTable, 'r'):
        if not line.startswith('#'):
            info = line.strip().split('\t')
            rna = info[0]
            tpmlist = info[1::]
            tpmlist = [float(f) for f in tpmlist]
            if min(tpmlist) >= 1:
                filtrna.append(rna)
    return filtrna

#dictionary of the rna from the hek293t data that points from rna to the expression count 
def rnatohekexp(hektpmtable, orthologue_f):
    humantomouse = {}
    for line in open(orthologue_f, 'r'):
        info = line.strip().split()
        humangene = info[1]
        mousegene = info[3]
        humantomouse[humangene] = mousegene
    rnatoexp = {}
    for line in open(hektpmtable, 'r'):
        if not line.startswith('#'):
            info = line.strip().split('\t')
            gene = info[0]
            if gene in humantomouse:
                mousegene = humantomouse[gene]
                tpm = info[1:3]
                tpm = [float(f) for f in tpm]
                meantpm = mean(tpm)
                rnatoexp[mousegene] = meantpm
    return rnatoexp

#creates dictionary from the human postar.csv downloaded file the dictionary is the rna pointing to the number of clip binding sites for each rna for which there is binding site data in the postar table
def clipvalidatesrange(postarrna_f, orthologue_f):
    with open(postarrna_f, 'r') as f:
        next(f)
        rnatobs = []
        for line in f:
            rna = line.strip().split(',')[0]
            rna = rna.split('"')[1]
            bindingsites = line.strip().split(',')[4].split('"')[1]
            bindingsites = int(bindingsites)
            rnatobs.append([rna, bindingsites])
    humantomouse = {}
    for line in open(orthologue_f, 'r'):
        info = line.strip().split()
        humangene = info[1]
        mousegene = info[3]
        humantomouse[humangene] = mousegene
    rnatobsdict = {}
    filtlist = []                                                                                  
    for rna in rnatobs:
        if rna[0] in humantomouse:                                                                 
            filtlist.append(humantomouse[rna[0]])                                                  
            key = humantomouse[rna[0]]
            rnatobsdict[key] = rna[1]
    return rnatobsdict

#plots Figure 6a
def plotproteinrna(tpmTable, ms_table, rbp, rna):
    data, genelist = tpm_array(tpmTable)
    msdata, protlist = ms_array3(ms_table)
    fraction = range(0,17)
    plt.figure()
    idx = genelist.index(rna)
    profile = data[idx]
    profile = [f/sum(profile[0:17]) for f in profile[0:17]]
    idx2 = protlist.index(rbp)
    proteinprof = msdata[idx2][0:17]
    fig,ax1 = plt.subplots()
    ax2 = ax1.twinx()
    profile = [f*17 for f in profile]
    ax1.plot(fraction, profile)
    proteinprof = [f/sum(proteinprof[0:17]) for f in proteinprof[0:17]]
    proteinprof = [f*17 for f in proteinprof[0:17]]
    ax2.plot(fraction, proteinprof, linestyle='dashed', dashes=(5, 5), color = '#FF0200')
    print corrcoef(profile, proteinprof)
    ax1.spines['top'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax1.yaxis.set_ticks(np.arange(0.5, 5.5, 1))
    ax2.yaxis.set_ticks(np.arange(0.5, 3, .5))
    plt.savefig('test.pdf')

#makes the bins of rnas grouped by correlation to rbp that are used in figured 6
def definebinsfigure6(tpmTable, clust_file, uniprot_table, ms_table, rbp):
    msdata, protlist = ms_array3(ms_table)
    idx = protlist.index(rbp)
    pepprofile = msdata[idx]
    pepprofile = [f/sum(pepprofile[0:17]) for f in pepprofile[0:17]]
    data, rowHeaders = tpm_array(tpmTable)
    filtrna = filterrnas(tpmTable)
    matchlenrnalist = []
    for rna in filtrna:
        idx = rowHeaders.index(rna)
        tpm = data[idx]
        a =  tpm[0:17]
        a = [float(f) for f in a]
        a = [x/sum(a) for x in a]
        matchlenrnalist.append(a)
    group1 = []
    group2 = []
    group3 = []
    group4 = []
    group5 = []
    correlationspread = []
    i = 0
    for profile in matchlenrnalist:
        correlation = corrcoef(profile, pepprofile)[0][1]
        correlationspread.append(correlation)
        rna = filtrna[i]
        if correlation >0.85:
            group1.append(rna)
        if correlation > 0.5:
            if correlation < 0.85:
                group2.append(rna)
        if correlation > 0:
            if correlation < 0.5:
                group3.append(rna)
        if correlation >-0.5:
            if correlation < 0:
                group4.append(rna)
        if correlation > -1:
            if correlation < -0.5:
                group5.append(rna)
        i +=1
    return group1, group2, group3, group4, group5
def bindatawrapper():
    indir = '/home/dadeku89/scratch/lib_Daa/RBP_analysis/postarfiles'
    files = [f for f in os.listdir(indir) if f.endswith('.csv')]
    for f in files:
        rbpname = f.split('.')[0]
        cmd = 'module load python \n'
        cmd = cmd + 'python script.py atlasseqfigures.bindata /home/dadeku89/scratch/lib_Daa/tpm_tables/tpm_table_3_24_filtered_al1.txt /home/dadeku89/scratch/lib_Daa/mass_spec/21965987_ms.txt~ %s %s APEXanalysis/ensemblhumantomouseorthologues.txt tpm_tables/hek293totalrnaseqtpmtable.txt /home/dadeku89/scratch/lib_Daa/RBP_analysis/postarplots/'%(rbpname, os.path.join(indir, f))
        scriptOptions = {'jobname': 'bindataplot', 'walltime': '00:30:00', 'memory': '5gb'}
        print cmd
        Mysbatch.launchJob(cmd, scriptOptions, verbose = True)
                                                                                                                                                                                                                                                                                                                                                            
        
             
    

def bindata(tpmTable, ms_table, rbp, postarrna_f, orthologue_f, hektpmtable, outdir):
    rnatobsdict= clipvalidatesrange(postarrna_f, orthologue_f)
    rnatoexp = rnatohekexp(hektpmtable, orthologue_f)
    msdata, protlist = ms_array3(ms_table)
    filtrna = filterrnas(tpmTable)
    idx = protlist.index(rbp)
    pepprofile = msdata[idx]
    pepprofile = [f/sum(pepprofile[0:17]) for f in pepprofile[0:17]]
    data, rowHeaders = tpm_array(tpmTable)
    matchlenrnalist = []
    for rna in filtrna:
        idx = rowHeaders.index(rna)
        tpm = data[idx]
        a =  tpm[0:17]
        a = [float(f) for f in a]
        a = [x/sum(a) for x in a]
        matchlenrnalist.append(a)
    
    i = 0
    correlationspread = []
    for profile in matchlenrnalist:
        correlation = corrcoef(profile, pepprofile)[0][1]
        rna = rowHeaders[i]
        correlationspread.append([correlation, rna])
        i +=1
    sortcorr = sorted(correlationspread, key = itemgetter(0))
    corrlist = []
    rnalist = []
    for x in sortcorr:
        corrlist.append(x[0])
        rnalist.append(x[1])
    n = 713
    groupmeans = []
    grouptornalist = {}
    for item in  [corrlist[i:i + n] for i in range(0, len(corrlist), n)]:
        groupmeans.append(mean(item))
    keyidx = 0
    for item in  [rnalist[i:i + n] for i in range(0, len(rnalist), n)]:
        #print item, len(item)#, mean(item)
        grouptornalist[groupmeans[keyidx]] = item
        keyidx +=1
    #rnalist = []
    corrdist = []
    for rna in rnatobsdict:
        if rna in rnalist:
            corrdist.append(corrlist[rnalist.index(rna)])
    print corrdist
    bins = linspace(-1,1,25)
    plt.hist(corrdist, bins)
    plt.xlabel('Correlation')
    plt.savefig(os.path.join(outdir,rbp+'corrdistcliprna.pdf'))
    """if rnatobsdict[rna] > 1:
            rnalist.append(rna)
    for group in grouptornalist:
        i = 0
        expout = 0 
        for rna in rnalist:
            if rna in grouptornalist[group]:
                i +=1
        explist = []
        for rna in grouptornalist[group]:
            if rna in rnatoexp:
                explist.append(rnatoexp[rna])
            else:
                expout +=1
        meanexp = mean(explist)
        print float(i)/log(meanexp), group , expout"""
        
def figure6fwrapper():
    indir = '/home/dadeku89/scratch/lib_Daa/RBP_analysis/postarfiles'
    outdir = '/home/dadeku89/scratch/lib_Daa/RBP_analysis/postarplots'
    files = [f for f in os.listdir(indir) if f.endswith('.csv')]
    
    for f in files:
        rbpname = f.split('.')[0]
        cmd = 'module load python \n'
        cmd = cmd + 'python script.py atlasseqfigures.figure6f /home/dadeku89/scratch/lib_Daa/tpm_tables/tpm_table_3_24_filtered_al1.txt /home/dadeku89/scratch/lib_Daa/cluster_files/cluster_3_24_t_0.055_i_1.txt /home/dadeku89/scratch/lib_Daa/mm10_protacc.txt /home/dadeku89/scratch/lib_Daa/mass_spec/21965987_ms.txt~ %s %s /home/dadeku89/scratch/lib_Daa/APEXanalysis/ensemblhumantomouseorthologues.txt /home/dadeku89/scratch/lib_Daa/tpm_tables/hek293totalrnaseqtpmtable.txt %s'%(rbpname, os.path.join(indir,f), outdir)
        scriptOptions = {'jobname': '6f', 'walltime': '00:30:00', 'memory': '5gb'}
        print cmd
        Mysbatch.launchJob(cmd, scriptOptions, verbose = True)
        
#plot cdf for 6f, outputs pvalues #delete outdir from functino to change how the file is saved in the future
def figure6f(tpmTable, clust_file, uniprot_table, ms_table, rbp, postarrna_f, orthologue_f, hektpmtable, outdir):
    groups = definebinsfigure6(tpmTable, clust_file, uniprot_table, ms_table, rbp)
    grouplist = ['0.85 to 1.00', '0.50 to 0.85', '0.00 to 0.50', '-0.50 to 0.00', '-1.00 to -0.50']
    rnatobsdict= clipvalidatesrange(postarrna_f, orthologue_f)
    rnatoexp = rnatohekexp(hektpmtable, orthologue_f)
    #bins = linspace(0,2, 2000) paper paramters
    bins = linspace(0,10,2000)
    topbs = []
    topct = []
    plt.figure()
    bsdict = {}
    ax = plt.subplot(1,1,1)
    cmap = matplotlib.cm.get_cmap('viridis')
    clist = []
    c = 0.1
    i = 0
    for group in groups: 
        groupbs = [] 
        for rna in group:
            if rna in rnatobsdict: 
                bsnum = rnatobsdict[rna]
                if rna in rnatoexp:
                    exp = rnatoexp[rna]
                    if bsnum/exp != inf:
                        #groupbs.append(float(bsnum/log(exp)))
                        groupbs.append(float(bsnum))
            else:
                groupbs.append(0.0)
        plt.hist(groupbs, bins, histtype='step', cumulative = True, normed = True, color = cmap(c),linewidth=2, label = grouplist[i])
        if i ==1:
            print groupbs, bins
            print percentile(groupbs, 95)
        bsdict[grouplist[i]] = groupbs
        clist.append(c)
        c +=0.2
        i +=1
    plt.xlabel('Number of Binding Sites normalized by expression from Hek293T rnaseq')       
    plt.ylabel('Normalized Number of RNAs for each group') 
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    fix_hist_step_vertical_line_at_end(ax)
    #plt.ylim(.7,1)
    plt.ylim(0.8,1)# this is the parameters for the figure 6f figure
    #plt.xticks(np.arange(0, 2, 0.5))
    #plt.yticks(np.arange(0.7, 1.0, 0.1))
    #plt.yticks(np.arange(0.7, 1.0, 0.05))
    plt.legend(loc = 'lower right')                                                               
    for a in grouplist:
        pvaluelist = []
        for b in grouplist:
            #print ranksums(bsdict[a], bsdict[b]), a, b, len(bsdict[a]), len(bsdict[b])
            #print ks_2samp(bsdict[a], bsdict[b]), a, b
            pvaluelist.append(-log10(ranksums(bsdict[a], bsdict[b])[1]))
    #plt.savefig('test.pdf')
    plt.savefig(os.path.join(outdir, rbp+'unnorm.pdf'))
    
    return grouplist, bsdict, clist

#make the dictionary of bin/group pointing to clip binding site list for each rna of the group and also ouputs the colorindex list
def figure6fdict(tpmTable, clust_file, uniprot_table, ms_table, rbp, postarrna_f, orthologue_f, hektpmtable):
    groups = definebinsfigure6(tpmTable, clust_file, uniprot_table, ms_table, rbp)
    grouplist = ['0.85 to 1.00', '0.50 to 0.85', '0.00 to 0.50', '-0.50 to 0.00', '-1.00 to -0.50']
    rnatobsdict= clipvalidatesrange(postarrna_f, orthologue_f)
    rnatoexp = rnatohekexp(hektpmtable, orthologue_f)
    bins = linspace(0,2, 100)
    topbs = []
    topct = []
    bsdict = {}
    cmap = matplotlib.cm.get_cmap('viridis')
    clist = []
    c = 0.1
    i = 0
    for group in groups:
        groupbs = []
        for rna in group:
            if rna in rnatobsdict:
                bsnum = rnatobsdict[rna]
                if rna in rnatoexp:
                    exp = rnatoexp[rna]
                    if bsnum/exp != inf:
                        groupbs.append(float(bsnum/log(exp)))
            else:
                groupbs.append(0.0)
        bsdict[grouplist[i]] = groupbs
        clist.append(c)
        c +=0.2
        i +=1
    return grouplist, bsdict, clist

#plots the bargraph of the 95th percentile for the bins, will be in the inset of figure 6f
def bar95percentile(tpmTable, clust_file, uniprot_table, ms_table, rbp, postarrna_f, orthologue_f, hektpmtable,):
    grouplist, bsdict, clist = figure6fdict(tpmTable, clust_file, uniprot_table, ms_table, rbp, postarrna_f, orthologue_f, hektpmtable)
    plt.figure()
    ax = plt.subplot(1,1,1)
    cmap = matplotlib.cm.get_cmap('viridis')
    c = 0
    for group in reversed(grouplist):
        upperpercentile = percentile(bsdict[group], 95)
        cindex = grouplist.index(group)
        plt.bar(group, upperpercentile, color = cmap(clist[cindex]))
        print group, upperpercentile
        #print ranksums(bsdict[group], bsdict[grouplist[-2]])                                       
        #print  ranksums(percentile(bsdict[group], 95), percentile(bsdict[grouplist[-1]], 95))      
        c +=1
    plt.xlabel('Group')
    plt.ylabel('95th percentile of binding sites (normalized by expression)')
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(True)
    ax.spines['right'].set_visible(False)
    #plt.ylim(0,0.6)
    plt.ylim(0,0.75)
    #plt.yticks(np.arange(0, 0.75, 0.25))
    plt.savefig('test.pdf')

#plots the cdf of the expression from the hek293t data to see the expression bias of each of the bins used for the clip validation analysis
def expressionclip(tpmTable, clust_file, uniprot_table, ms_table, rbp, postarrna_f, orthologue_f, hektpmtable):
    groups = definebinsfigure6(tpmTable, clust_file, uniprot_table, ms_table, rbp)
    grouplist = ['0.85 to 1.00', '0.50 to 0.85', '0.00 to 0.50', '-0.50 to 0.00', '-1.00 to -0.50']
    rnatoexp = rnatohekexp(hektpmtable, orthologue_f)
    i = 0
    bins = linspace(0,600,250)
    #plt.figure(figsize = (7,7))                                                                    
    cmap = matplotlib.cm.get_cmap('viridis')
    c = 0.1
    ax = plt.subplot(1,1,1)
    for group in groups:
        explist = []
        for rna in group:
            if rna in rnatoexp:
                exp = rnatoexp[rna]
                explist.append(exp)
        medexp = median(explist)
        plt.hist(explist, bins, histtype='step', cumulative = True, normed = True, color = cmap(c),linewidth=2, label = grouplist[i])
        print max(explist)
        i +=1
        c +=0.2
    plt.xlabel('Group')
    plt.ylabel('TPM counts')
    fix_hist_step_vertical_line_at_end(ax)
    plt.legend(loc = 'lower right')
    plt.savefig('test.pdf')

def rbptocluster(tpmTable, clust_f, mannrbp_f, rbpdb_f, orthologue_f, ms_table):
    clustmatchmsprof = clustGenerp(tpmTable, clust_f)
    rbpinms = generaterbpsinms(mannrbp_f, rbpdb_f, orthologue_f, ms_table)
    msdata, protlist = ms_array3(ms_table)
    rbptoclustercorr = {}
    rbptoclustercorrtop5 = {}
    sortrbpsbycorr = []
    i = 0
    nodups = []
    for rbp in rbpinms:
        if rbp not in nodups:
            nodups.append(rbp)
            idx = protlist.index(rbp)
            proteinprofile = msdata[idx][0:17]
            proteinprofile = [f/sum(proteinprofile) for f in proteinprofile]
            corrlist = []
            for cluster in clustmatchmsprof:
                corr = corrcoef(proteinprofile, clustmatchmsprof[cluster])[0][1]
                corrlist.append([corr, cluster])
            sortedcorr = sorted(corrlist, key = itemgetter(0),reverse = True)
            rbptoclustercorrtop5[rbp] = sortedcorr[0:5]
            rbptoclustercorr[rbp] = sortedcorr[0]
            if sortedcorr[0][0] > 0.8:
                i +=1

            sortrbpsbycorr.append([rbp, sortedcorr[0][0]])
    sortedtable = sorted(sortrbpsbycorr, key = itemgetter(1), reverse = True)
    rbpssortedtable = []
    for item in sortedtable:
        rbpssortedtable.append(item[0])
    print i
    return rbptoclustercorr, rbptoclustercorrtop5, rbpssortedtable

def rbptotoprnas(tpmTable,mannrbp_f, rbpdb_f, orthologue_f, ms_table, rnatorbp_f):
    data, genelist = rnalistrbpmatch(tpmTable)
    rbpinms = generaterbpsinms(mannrbp_f, rbpdb_f, orthologue_f, ms_table)
    msdata, protlist = ms_array3(ms_table)
    out = open(rnatorbp_f, 'w')
    out.write('#RBP' + '\t' + 'RNA' + '\t' + 'Pearson Correlation Coefficient' + '\n')
    for rbp in sorted(rbpinms):
        pidx = protlist.index(rbp)
        proteinprofile = msdata[pidx][0:17]
        corrlist = []
        i = 0
        for rna in data:
            corr = corrcoef(rna, proteinprofile)[0][1]
            if corr >= 0.8:
                corrlist.append([corr, rbp, genelist[i]])
            i +=1
        sortedcorr = sorted(corrlist, key = itemgetter(0), reverse = True)
        if len(sortedcorr) > 1:
            for item in sortedcorr:
                out.write(item[1] + '\t' + item[2] + '\t' + str(item[0]) + '\n') 
    out.close()

def countrbps(rnatorbp_f):
    rbpcount = []
    for line in open(rnatorbp_f, 'r'):
        if not line.startswith('#'):
            rbp = line.strip().split('\t')[0]
            if rbp not in rbpcount:
                rbpcount.append(rbp)
    print len(rbpcount)

def rnalistfrompostarcsv(postarrna_f, orthologue_f, hektpmtable):
    rnatobsdict= clipvalidatesrange(postarrna_f, orthologue_f)
    rnatoexp = rnatohekexp(hektpmtable, orthologue_f)
    rnalist = []
    for line in open(postarrna_f, 'r'):
        if not line.startswith('<U+FEFF>'):
            info = line.strip().split(',')
            hrna = info[0].strip('"')
            first = list(hrna)[0]
            second = list(hrna)[1::]
            rna = first+''.join(second).lower()
            rnatype = info[2].strip('"')

            if rnatype == 'protein_coding':
                rnalist.append(rna)
    for rna in rnalist:
            if rna in rnatobsdict:
                bsnum = rnatobsdict[rna]
                if rna in rnatoexp:
                    exp = rnatoexp[rna]
                    if bsnum/exp != inf:
                        #groupbs.append(float(bsnum/log(exp)))
                        print rna, float(bsnum/log(exp))
            #else:
            #    groupbs.append(0.0)

def tableS7(tpmTable, clust_f, mannrbp_f, rbpdb_f, orthologue_f, ms_table, topcluster_out, top5cluster_out):
    rbptoclustercorr, rbptoclustercorrtop5, rbpssortedtable = rbptocluster(tpmTable, clust_f, mannrbp_f, rbpdb_f, orthologue_f, ms_table)
    out1 = open(topcluster_out, 'w')
    out1.write('#RBP' + '\t' + 'Cluster' + '\t' + 'Pearson Correlation Coefficient' + '\n')
    for rbp in rbpssortedtable:
        out1.write(rbp  + '\t' + str(rbptoclustercorr[rbp][1]) + '\t' + str(rbptoclustercorr[rbp][0]) + '\n')

    out1.close()
    
    out2 = open(top5cluster_out, 'w')
    out2.write('#RBP' + '\t' + 'Cluster' + '\t' + 'Pearson Correlation Coefficient' + '\n')
    for rbp in  rbpssortedtable:
        for clust in rbptoclustercorrtop5[rbp]:                                                    
            out2.write(rbp + '\t' + str(clust[1]) + '\t' + str(clust[0]) + '\n')
    out2.close()
    
def rbplisthumantomouse(tpmTable, clust_f, mannrbp_f, rbpdb_f, orthologue_f, ms_table):
    rbptoclustercorr, rbptoclustercorrtop5, rbpssortedtable = rbptocluster(tpmTable, clust_f, mannrbp_f, rbpdb_f, orthologue_f, ms_table)
    allrbps = []
    for rbp in rbpssortedtable:
        allrbps.append(rbp)
    return allrbps

def lowerdbrbp(clipdbmouse_f):
    clipdbmouselist = []
    for line in open(clipdbmouse_f, 'r'):
        rbp = line.strip()
        first = list(rbp)[0] 
        second = list(rbp)[1::]
        clipdbmouselist.append(first+''.join(second).lower())
    
    return clipdbmouselist

def findrbps(tpmTable, clust_f, mannrbp_f, rbpdb_f, orthologue_f, ms_table, clipdbmouse_f):
    clipdbmouselist =  lowerdbrbp(clipdbmouse_f)
    allrbps = rbplisthumantomouse(tpmTable, clust_f, mannrbp_f, rbpdb_f, orthologue_f, ms_table)
    for item in allrbps:
        if item in clipdbmouselist:
            print item
    
def matchrbpstoortho(tpmTable, clust_f, mannrbp_f, rbpdb_f, orthologue_f, ms_table,  clipdbphumanfile):
    allrbps = rbplisthumantomouse(tpmTable, clust_f, mannrbp_f, rbpdb_f, orthologue_f, ms_table)
    clipdbphumanlist = []
    notin = []
    humanlist = []
    for line in open(clipdbphumanfile, 'r'):
        clipdbphumanlist.append(line.strip())
    humantomouse = {}
    for line in open(orthologue_f, 'r'):
        info = line.strip().split()
        humangene = info[1]
        mousegene = info[3]
        humantomouse[humangene] = mousegene
    for item in clipdbphumanlist:
        if item in humantomouse:
            humanlist.append(humantomouse[item])
        else:
            notin.append(item)
    for item in allrbps:
        if item in humanlist:
            print item
    for item in notin:
        print item
    
def rbpclustercorrelationheatmap(mannrbp_f, rbpdb_f, orthologue_f, tpmTable, clust_f, ms_table):
    fpath = os.path.join(rcParams["datapath"], "/home/dadeku89/.local/lib/python2.7/site-packages/matplotlib/mpl-data/fonts/ttf/helvetica.ttf")
    prop = fm.FontProperties(fname=fpath)
    fname = os.path.split(fpath)[1]
    msdata, protlist = ms_array3(ms_table)
    data, rowHeaders = tpm_array(tpmTable)
    #totalrbplist = rbplist(mannrbp_f, rbpdb_f, orthologue_f)
    clusterGene = clustGenerp(tpmTable, clust_f)
    rbpclustercorr = []
    totalrbplist = generaterbpsinms(mannrbp_f, rbpdb_f, orthologue_f, ms_table)
    clusterlist = []
    rbplabels = []
    orderedclustlist = sorted(clusterGene.keys())
    for rbp in protlist:
        if rbp in totalrbplist:
            if rbp not in rbplabels:
                rbplabels.append(rbp)
                clustercorr = []
                idx = protlist.index(rbp)
                pepprofile = msdata[idx]
                pepprofile = [f/sum(pepprofile[0:17]) for f in pepprofile[0:17]]
                for cluster in orderedclustlist:
                    if cluster not in clusterlist:
                        clusterlist.append(cluster)
                    clustercorr.append(corrcoef(pepprofile, clusterGene[cluster])[0][1])
                rbpclustercorr.append(clustercorr)
    rbpclustercorr = array(rbpclustercorr)
    print rbpclustercorr.shape
    plt.figure(figsize = (20,18))
    G = gridspec.GridSpec(3, 1, height_ratios=[1, 2, 3])
    Y = rbpclustercorr#.T
    Y = pdist(Y)
    Z = sch.linkage(Y, method ='average', metric ='correlation')
    ax = plt.subplot(G[0])
    sch.set_link_color_palette(['#000000'])
    R =sch.dendrogram(Z, no_labels=True, labels = None, above_threshold_color = '#000000')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='y', left = False)
    ax.axes.get_yaxis().set_ticks([])
    ax = plt.subplots_adjust(bottom = 0.01)
    idx = R['leaves']
    orderedmatrix = []
    orderedrbplist = []
    for i in idx:
        orderedmatrix.append(rbpclustercorr[i,:])
        orderedrbplist.append(rbplabels[i])
    orderedmatrix = array(orderedmatrix)
    orderedmatrix = orderedmatrix.T
    print orderedmatrix.shape
    
        
    ax2 = plt.subplot(G[1::])
    ax2.pcolor(orderedmatrix, vmin = -1, vmax = 1, cmap = 'Blues')
    plt.xticks(arange(0.5, len(rbplabels)), orderedrbplist, rotation = 90)
    plt.yticks(arange(0.5, len(clusterlist)), clusterlist)
    ax2.set_xticklabels(rbplabels, fontProperties=prop, rotation = 90, size = 11)
    ax2.set_yticklabels(clusterlist, fontProperties=prop, size = 11)
    ax2.xaxis.tick_top()
    #plt.colorbar()
    plt.savefig('test.pdf')

def colorbarplot(mannrbp_f, rbpdb_f, orthologue_f, tpmTable, clust_f, ms_table):
    fpath = os.path.join(rcParams["datapath"], "/home/dadeku89/.local/lib/python2.7/site-packages/matplotlib/mpl-data/fonts/ttf/helvetica.ttf")
    prop = fm.FontProperties(fname=fpath)
    fname = os.path.split(fpath)[1]
    msdata, protlist = ms_array3(ms_table)
    data, rowHeaders = tpm_array(tpmTable)
    clusterGene = clustGenerp(tpmTable, clust_f)
    rbpclustercorr = []
    totalrbplist = generaterbpsinms(mannrbp_f, rbpdb_f, orthologue_f, ms_table)
    clusterlist = []
    rbplabels = []
    orderedclustlist = sorted(clusterGene.keys())
    for rbp in protlist:
        if rbp in totalrbplist:
            if rbp not in rbplabels:
                rbplabels.append(rbp)
                clustercorr = []
                idx = protlist.index(rbp)
                pepprofile = msdata[idx]
                pepprofile = [f/sum(pepprofile[0:17]) for f in pepprofile[0:17]]
                for cluster in orderedclustlist:
                    if cluster not in clusterlist:
                        clusterlist.append(cluster)
                    clustercorr.append(corrcoef(pepprofile, clusterGene[cluster])[0][1])
                rbpclustercorr.append(clustercorr)
    rbpclustercorr = array(rbpclustercorr)
    #plt.figure(figsize = (20,18))
    Y = rbpclustercorr#.T                                                                                            
    Y = pdist(Y)
    Z = sch.linkage(Y, method ='average', metric ='correlation')
    
    R =sch.dendrogram(Z, no_labels=True, labels = None, above_threshold_color='#AAAAAA')
    """ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='y', left = False)
    ax.axes.get_yaxis().set_ticks([])
    ax = plt.subplots_adjust(bottom = 0.01)"""
    idx = R['leaves']
    orderedmatrix = []
    orderedrbplist = []
    for i in idx:
        orderedmatrix.append(rbpclustercorr[i,:])
        orderedrbplist.append(rbplabels[i])
    orderedmatrix = array(orderedmatrix)
    orderedmatrix = orderedmatrix.T
    plt.figure()
    plt.pcolor(orderedmatrix, vmin = -1, vmax = 1, cmap = 'Blues')
    """plt.xticks(arange(0.5, len(rbplabels)), orderedrbplist, rotation = 90)
    plt.yticks(arange(0.5, len(clusterlist)), clusterlist)
    ax2.set_xticklabels(rbplabels, fontProperties=prop, rotation = 90, size = 11)
    ax2.set_yticklabels(clusterlist, fontProperties=prop, size = 11)
    ax2.xaxis.tick_top()"""
    #plt.figure()
    plt.colorbar()                                                                                                   
    plt.savefig('test.pdf')

def clusterfractions(tpmTable):
    data = tpm_array(tpmTable)[0]
    data =[[x + .1 for x in row] for row in  data]
    data = log(data)
    data = array(data)
    fractions = range(3,25)
    print fractions
    corr = corrcoef(data.T)
    corr = array(corr)
    print corr.shape
    Y = pdist(corr)
    Z = sch.linkage(Y, method ='average', metric ='correlation')
    R =sch.dendrogram(Z,orientation= 'right', no_labels=True)
    idx = R['leaves']
    orderedmatrix = []
    for i in idx:
        print fractions[i]
    #don't include in github
        
    
def geoparse(infile, outfile):
    out = open(outfile, 'w')
    for line in open(infile, 'r'):
        info =line.strip().split()
        out.write(info[0] + '\t' + info[1] + '\n')

def maketrackDbfiles(bamdir, urlbase, outfile, gradient):
    files = [f for f in os.listdir(bamdir) if f.endswith('.bam')]
    out = open(outfile, 'w')
    grad1exclud = ['20','22','23','24']
    for f in files:
        if gradient == '2':
            fracgradnum = f.split('.')[0].split('Aligned')[0]
            gradnum = fracgradnum.split('_')[1]
            fracnum = fracgradnum.split('_')[0]
            track = 'fraction'+ fracgradnum
            url = urlbase + f
            shortlabel = 'Frac_' + fracnum + ' Grad_'+ gradnum
            longlabel = 'ATLAS-Seq Fraction '+ fracnum + ' Gradient '+ gradnum
            typelabel = 'bam'
            out.write('track' + ' ' + fracgradnum + '\n' + 'bigDataUrl' + ' ' + url  +'\n' + 'shortLabel' + ' '+ shortlabel +'\n' + 'longLabel' + ' '+ longlabel +'\n' + 'type' + ' '+ typelabel + '\n'+ '\n'+ '\n')
        if gradient == '1':
            fracgradnum = f.split('.')[0].split('Aligned')[0]
            gradnum = fracgradnum.split('_')[1]
            fracnum = fracgradnum.split('_')[0]
            if fracnum not in grad1exclud:
                track = 'fraction'+ fracgradnum
            url = urlbase + f
            shortlabel = 'Frac_' + fracnum + ' Grad_'+ gradnum
            longlabel = 'ATLAS-Seq Fraction '+ fracnum + ' Gradient '+ gradnum
            typelabel = 'bam'
            out.write('track' + ' ' + fracgradnum + '\n' + 'bigDataUrl' + ' ' + url  +'\n' + 'shortLabel' + ' '+ shortlabel +'\n' + 'longLabel' + ' '+ longlabel +'\n' + 'type' + ' '+ typelabel + '\n'+ '\n'+ '\n')
            
        
