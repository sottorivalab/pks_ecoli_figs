

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from collections import Counter
from functools import reduce
import seaborn as sns
from itertools import groupby
import os
import glob
from statannotations.Annotator import Annotator
import scipy.stats as stats
from statannot import add_stat_annotation
from matplotlib import gridspec
import statsmodels.api as sm

"""
mypath='/Users/bchen/Desktop/Projects/'
folder1='Nomral'
floder2='/11EPICC_Ecoli_files/'
new='/Users/bchen/Desktop/Projects/PKSFiles/Data/'
"""

######
#####Fig1E Rplot
library(ggplot2)
ggplot(testdf, aes(x = BMI_GP, y = Average_Value, colour = age_gp)) + geom_errorbar(aes(ymax = CI_UL, ymin = CI_LL), position = "dodge") + geom_point(position = position_dodge(0.9))


#df<-read.csv(paste0("/Users/bchen/Desktop/Projects/Normal/dNdSresult/summaryAll2.csv"))
df<-read.csv(paste0("/Users/bchen/Desktop/Projects/PKSFiles/Data/dNdSresult.csv"))

lwd_pt <- .pt
theme_set(theme_bw(base_size = 15, base_line_size = 3/lwd_pt, base_rect_size = 3/lwd_pt))

myplot1<-ggplot(df, aes(x=Genes, y=mle, colour=Genes)) + 
  geom_point(size=3) + 
  geom_errorbar(aes(ymin = cilow, ymax = cihigh),size =1.5,width=0.5)+
  facet_grid(type~name)+scale_y_continuous(trans='log10',limits = c(0.05,20))+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  labs(title="dN/dS of Normals",y="dN/dS")

ggsave("/Users/bchen/Desktop/Projects/PKSFiles/Figures/Fig1E.pdf", plot = myplot1)



###Fig1F python

ll=pd.read_csv("/Users/bchen/Desktop/Projects/PKSFiles/Data/mutationLoadGroup.csv",index_col='Unnamed: 0')


dlist=list(ll[ll['dlist']>0]['dlist'])
alist=list(ll[ll['alist']>0]['alist'])
clist=list(ll[ll['clist']>0]['clist'])
titles=['Distant Normals','Adjacent Normals','EPICC Cancer']



fig = plt.figure(figsize =(5, 5))
ax = fig.add_subplot()

bp = ax.boxplot(lll, patch_artist = True,notch ='False')

colors = [ 'orange','cyan','red']
 
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
    
ax.set_yscale('log')
ax.set_ylim([0,200000])
ax.set_xticklabels(titles)
plt.ylabel('Number of Mutations in log10 scale')

# Significance bars
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)

    # Colour of the median lines
    plt.setp(bp['medians'], color='k')

    # Check for statistical significance
    significant_combinations = []
    # Check from the outside pairs of boxes inwards
    combinations = [(1, 3), (2, 3), (1, 2)]
    for c in combinations:
        data1 = lll[c[0] - 1]
        data2 = lll[c[1] - 1]
        # Significance
        U, p = stats.mannwhitneyu(data1, data2, alternative='two-sided')
        if p < 0.05:
            significant_combinations.append([c, p])

# Get the y-axis limits
bottom, top = ax.get_ylim()
y_range = top - bottom
for i, significant_combination in enumerate(significant_combinations):
    # Columns corresponding to the datasets of interest
    x1 = significant_combination[0][0]
    x2 = significant_combination[0][1]
    # What level is this bar among the bars above the plot?
    level = len(significant_combinations) - i
    # Plot the bar
    bar_height = (y_range * 0.2 * level) + top*0.2
    bar_tips = bar_height - (y_range*0.2 * 0.2)
    plt.plot(
        [x1, x1, x2, x2],
        [bar_tips, bar_height, bar_height, bar_tips], lw=1, c='grey'
    )
    # Significance level
    p = significant_combination[1]
    if p < 0.001:
        sig_symbol = '***'
    elif p < 0.01:
        sig_symbol = '**'
    elif p < 0.05:
        sig_symbol = '*'
    text_height = bar_height + (y_range * 0.01)
    plt.text((x1 + x2) * 0.5, text_height, sig_symbol, ha='center', va='bottom', c='grey')
    bottom, top = ax.get_ylim()
    yrange = top - bottom
    ax.set_ylim(bottom - 0.002 * yrange, top)

ax.set_ylim([0,210000])
plt.savefig("/Users/bchen/Desktop/Projects/PKSFiles/Figures/Figure1F.pdf")
plt.show()




#####
#####FIg.2 A and B
dfall1=pd.read_csv("/Users/bchen/Desktop/Projects/PKSFiles/Data/SignatureSPSall1.csv")
dfAll2=pd.read_csv("/Users/bchen/Desktop/Projects/PKSFiles/Data/SignatureSPSall2.csv")


fig=plt.figure(figsize=(10,5))
gs = gridspec.GridSpec(1, 2, width_ratios=[4, 1]) 
ax = plt.subplot(gs[0])
ax = sns.barplot(x='pid', y='Percentage',hue='type',data=dfall1, palette=['darkred','pink',"#00CDCD","#FF7F00"],hue_order=['Cancer clonal','Cancer subclonal', 'adjacent normal', 'distant normal'], errwidth=0.8)
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, ncol=1, loc='upper center', bbox_to_anchor=(0.85, 1), frameon=True)
plt.xticks(rotation=90)
for xtick in ax.get_xticklabels():
    if xtick.get_text() in MSI: 
        xtick.set_color("black")
    else:
        xtick.set_color("black")
plt.ylim([0,0.5])
plt.title("proportion of SPS7 in EPICC samples")
ax.text(-0.03, 1.2, "A", transform=ax.transAxes,fontsize=16, fontweight='bold', va='top', ha='right')
 
ax = plt.subplot(gs[1])
ax = sns.boxplot(x='source', y='Percentage',data=dfAll2, palette=['darkred','pink',"#00CDCD","#FF7F00",'grey'],order=['Cancer clonal','Cancer subclonal', 'EPICC Adjacent Normal', 'EPICC Distant Normal','Healthy Normal (Lee-Six et al.)'])
pairs=[('Cancer clonal', 'Healthy Normal (Lee-Six et al.)'),('EPICC Distant Normal', 'Healthy Normal (Lee-Six et al.)')]
annotator = Annotator(ax, pairs,x='source', y='Percentage',data=dfAll2,order=['Cancer clonal','Cancer subclonal', 'EPICC Adjacent Normal', 'EPICC Distant Normal','Healthy Normal (Lee-Six et al.)'])
annotator.configure(test="Mann-Whitney")
annotator.apply_and_annotate()
plt.ylabel("proportion of SPS7 in samples")
plt.xlabel("Samples")
plt.xticks(rotation=90)
plt.title("")
plt.tight_layout()
ax.text(-0.13, 1.2, "B", transform=ax.transAxes,fontsize=16, fontweight='bold', va='top', ha='right')

plt.tight_layout()
plt.savefig('/Users/bchen/Desktop/Projects/PKSFiles/Figures/Figure2AB.pdf')
# 


#####Fig.2 C and D
dfdel=pd.read_csv("/Users/bchen/Desktop/Projects/PKSFiles/Data/proportion_of_short_tdel_allsample.csv")

fig=plt.figure(figsize=(10,5))
gs = gridspec.GridSpec(1, 2, width_ratios=[4, 1]) 
ax = plt.subplot(gs[0])
ax = sns.barplot(x=dfdel['pid'],y='proportion of short T-del',hue='group',data=dfdel, palette=['darkred','pink',"#00CDCD","#FF7F00"],hue_order=['EPICC Cancer clonal','EPICC Cancer subclonal', 'EPICC Adjacent Normal', 'EPICC Distant Normal'], errwidth=0.8)
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, ncol=1, loc='upper center', bbox_to_anchor=(0.4, 1), frameon=True,fontsize=8)
plt.xticks(rotation=90)
for xtick in ax.get_xticklabels():
    if xtick.get_text() in MSI: 
        xtick.set_color("black")
    else:
        xtick.set_color("black")
plt.ylim([0,1])
plt.title("proportion of T-del in EPICC samples")
ax.text(-0.03, 1.2, "C", transform=ax.transAxes,fontsize=16, fontweight='bold', va='top', ha='right')
 
ax = plt.subplot(gs[1])
ax =sns.boxplot(x='group',y='proportion of short T-del',data=dfdel,palette=['darkred','pink',"#00CDCD","#FF7F00"],hue_order=['EPICC Cancer clonal','EPICC Cancer subclonal', 'EPICC Adjacent Normal', 'EPICC Distant Normal'])

pairs=[('EPICC Cancer clonal', 'EPICC Adjacent Normal'),('EPICC Cancer clonal', 'EPICC Distant Normal')]
# 
annotator = Annotator(ax, pairs,x='group',y='proportion of short T-del',data=dfdel)
annotator.configure(test="Mann-Whitney")
annotator.apply_and_annotate()

plt.ylabel("proportion of  T-del  in samples")
plt.xlabel("Samples")
plt.xticks(rotation=90)
plt.title("")
plt.tight_layout()
ax.text(-0.13, 1.2, "D", transform=ax.transAxes,fontsize=16, fontweight='bold', va='top', ha='right')
plt.tight_layout()
plt.savefig('/Users/bchen/Desktop/Projects/PKSFiles/Figures/Figure2CD.pdf')
plt.show()



########
########Figure 3
df4=pd.read_csv("/Users/bchen/Desktop/Projects/PKSFiles/Data/EPICCpksReads.csv",index_col='Unnamed: 0')
pids=['C516','C518','C519','C524','C527',
      'C528','C530','C531','C536','C537',
      'C538','C539','C542','C543','C544',
      'C547','C548','C549','C550','C551',
      'C552','C555','C560','C561','C562']


with PdfPages('/Users/bchen/Desktop/Projects/PKSFiles/Figures/Figure3.pdf') as pdf:
    fig=plt.figure(figsize=(25,18))
    axi=1
    for pid in pids:               

        ax=fig.add_subplot(5,6,axi)
        axi+=1
        if axi%6==0:
            axi+=1
        else:
            pass 
        plt.title("%s"%pid,fontsize=18)
        dfi=df4[[x.split("_")[0]==pid for x in df4.index]]
        dfi=dfi[["Z1" not in x for x in dfi.index]] 
        clrs = ["red" if x=='Yes' else "grey" for x in list(dfi['pks+island'])]

        dfiorder=[]
        for i in range(len(dfi)):
            if "_".join(dfi.index[i].split("_")[1:3]) in Pnormal[pid]:
                if dfi.index[i].split("_")[1][0]=='E':
                    dfiorder.append(1)
                else:
                    dfiorder.append(2)
            else:   
                if dfi.index[i].split("_")[1][0]=='A':
                    dfiorder.append(3)
                elif dfi.index[i].split("_")[1][0]=='B':
                    dfiorder.append(4)
                elif dfi.index[i].split("_")[1][0]=='C':
                    dfiorder.append(5)
                elif dfi.index[i].split("_")[1][0]=='D':
                    dfiorder.append(6)               
                elif dfi.index[i].split("_")[1][0]=='F':
                    dfiorder.append(7)            
                elif dfi.index[i].split("_")[1][0]=='G':
                    dfiorder.append(8)         
                elif dfi.index[i].split("_")[1][0]=='H':
                    dfiorder.append(9)
                
        dfi['order']=dfiorder
        ax=sns.barplot(x=dfi.index,y='AeadsCount',data=dfi,palette=clrs,order=dfi['order'].sort_values().index)
        

        
        if pid in Pnormal.keys():
            labels = [item.get_text() for item in ax.get_xticklabels()]
            labels = ["_".join(x.split("_")[1:]) for x in labels] 


            ax.set_xticklabels(labels)
        else:
            labels = [item.get_text() for item in ax.get_xticklabels()]
            labels = ["_".join(x.split("_")[1:]) for x in labels]
            ax.set_xticklabels(labels)
    
        for xtick in ax.get_xticklabels():

            if pid in Pnormal.keys():

                if "_".join(xtick.get_text().split("_")[0:2]) in Pnormal[pid]:
                    if "_".join(xtick.get_text().split("_")[0:2]).startswith("E"):
                        xtick.set_color("#FF7F00")
                    else:
                        xtick.set_color("#00CDCD")
                elif xtick.get_text().split("_")[0].startswith("A"):
                    xtick.set_color("#E41A1C")
                elif xtick.get_text().split("_")[0].startswith("B"):
                    xtick.set_color("#377EB8")
                elif xtick.get_text().split("_")[0].startswith("C"):
                    xtick.set_color("#4DAF4A")      
                elif xtick.get_text().split("_")[0].startswith("D"):
                    xtick.set_color("#984EA3")
                else:
                    xtick.set_color("black")
                    
            elif xtick.get_text().split("_")[0].startswith("E"):
                xtick.set_color("#FF7F00")
            elif xtick.get_text().split("_")[0].startswith("A"):

                xtick.set_color("#E41A1C")
            elif xtick.get_text().split("_")[0].startswith("B"):

                xtick.set_color("#377EB8")
            elif xtick.get_text().split("_")[0].startswith("C"):

                xtick.set_color("#4DAF4A")      
            elif xtick.get_text().split("_")[0].startswith("D"):

                xtick.set_color("#984EA3")
            else:

                xtick.set_color("black")
            
        plt.xticks(rotation=90)
        plt.ylabel("Read counts of E.coli")


    plt.tight_layout()    
    pdf.savefig()
    plt.close()
 

#####
#####Fig. 4
#%%

##only SPS7
#group for driver gene
pidlist=[['C519','C524','C527','C528','C530','C549','C531','C537','C543','C544','C547','C560','C550','C551']]

for pids in pidlist:
    dflist=[]
    for pid in pids:
        df=pd.read_csv("/Users/bchen/Desktop/Projects/PKSFiles/Data/drivercmgSBS7/%s_driver_cmg_InSBS7context.csv"%pid,index_col='Unnamed: 0')
        df=df[df['region']!='intron_variant']
        df=df[df['region']!='intron_variant&non_coding_transcript_variant']
        df=df[df['region']!='intron_variant&NMD_transcript_variant']
        df=df[df['region']!='synonymous_variant']
        df=df[df['region']!='non_coding_transcript_exon_variant']

        df=pd.merge(df,dgl,how='inner',left_on='driver_genes',right_on='driver_genes')   
  
        df.index=df['driver_genes']+"_"+df['sample']
        df=df.groupby('driver_genes').mean()
        dfx=df[['SBS1|Patient', 
           'SBS2|Patient',  'SBS3|Patient', 
           'SBS4|Patient', 'SBS5|Patient', 
           'SBS6|Patient', 'SBS7|Patient']]
        dfx.columns=['SBS1','SBS2','SBS3','SBS4','SBS5','SBS6','SBS7']
        contribution=[]
        for i in range(len(dfx)):
            aaa=int(dfx.iloc[i].sort_values(ascending=False).index[0][-1])
            contribution.append(-aaa)
        dfx['contribution']=contribution
        dfx.sort_values(by=['contribution','SBS1','SBS2','SBS3','SBS4','SBS5','SBS6','SBS7'],ascending=False,inplace=True)
        dfx.drop('contribution',axis=1,inplace=True)
        print(len(dfx))
        months = list(dfx.index)
        seasons = [pid]*len(dfx)
        tuples = list(zip(months, seasons))
        index = pd.MultiIndex.from_tuples(tuples, names=['first', 'second'])     
        dfx.index=index
        dflist.append(dfx)
        
    dfalldriver= reduce(lambda  left,right: pd.concat([left,right]), dflist)
    
for pids in pidlist:
    dflist=[]
    for pid in pids:
        df=pd.read_csv("/Users/bchen/Desktop/Projects/PKSFiles/Data/drivercmgSBS7/%s_driver_cmg_InSBS7context.csv"%pid,index_col='Unnamed: 0')
        df=df[df['region']!='intron_variant']
        df=df[df['region']!='intron_variant&non_coding_transcript_variant']
        df=df[df['region']!='intron_variant&NMD_transcript_variant']
        df=df[df['region']!='synonymous_variant']
        df=df[df['region']!='non_coding_transcript_exon_variant']

        df=pd.merge(df,cmgl,how='inner',left_on='driver_genes',right_on='cmgs')

  
        df.index=df['driver_genes']+"_"+df['sample']
        df=df.groupby('driver_genes').mean()
        dfx=df[['SBS1|Patient', 
           'SBS2|Patient',  'SBS3|Patient', 
           'SBS4|Patient', 'SBS5|Patient', 
           'SBS6|Patient', 'SBS7|Patient']]
        dfx.columns=['SBS1','SBS2','SBS3','SBS4','SBS5','SBS6','SBS7']
        contribution=[]
        for i in range(len(dfx)):
            aaa=int(dfx.iloc[i].sort_values(ascending=False).index[0][-1])
            contribution.append(-aaa)
        dfx['contribution']=contribution
        dfx.sort_values(by=['contribution','SBS1','SBS2','SBS3','SBS4','SBS5','SBS6','SBS7'],ascending=False,inplace=True)
        dfx.drop('contribution',axis=1,inplace=True)
        print(len(dfx))
        months = list(dfx.index)
        seasons = [pid]*len(dfx)
        tuples = list(zip(months, seasons))
        index = pd.MultiIndex.from_tuples(tuples, names=['first', 'second'])     
        dfx.index=index
        dflist.append(dfx)
        
    dfallcmg= reduce(lambda  left,right: pd.concat([left,right]), dflist)



fig = plt.figure(figsize=(17,10))
ax = fig.add_subplot(1,2,1)
sns.heatmap(dfalldriver,cmap="coolwarm",cbar_kws={"shrink": 0.5},square=True, cbar=False)    
plt.gcf().set_size_inches(17,10)    
labels = ['' for item in ax.get_yticklabels()]
ax.set_yticklabels(labels)
ax.set_ylabel('')
ax.set_xticklabels(ax.get_xmajorticklabels(), fontsize = 18,rotation =90)
ax.text(-1, 1, "A", transform=ax.transAxes,fontsize=16, fontweight='bold', va='top', ha='right')
 
label_group_bar_table(ax, dfalldriver,fontsizegene=10,fontsizepatient=12)
fig.subplots_adjust(bottom=.1*df.index.nlevels)

ax = fig.add_subplot(1,2,2)
sns.heatmap(dfallcmg,cmap="coolwarm",cbar_kws={"shrink": 0.5},square=True, cbar=False)    
plt.gcf().set_size_inches(17,10)    
labels = ['' for item in ax.get_yticklabels()]
ax.set_yticklabels(labels)
ax.set_ylabel('')
ax.set_xticklabels(ax.get_xmajorticklabels(), fontsize = 18,rotation =90)    
ax.text(-1, 1, "B", transform=ax.transAxes,fontsize=16, fontweight='bold', va='top', ha='right')

label_group_bar_table(ax, dfallcmg,fontsizegene=10,fontsizepatient=12)
fig.subplots_adjust(bottom=.1*df.index.nlevels)

plt.tight_layout()

plt.savefig('/Users/bchen/Desktop/Projects/PKSFiles/Figures/Figure4.pdf')
plt.show()






#######
######Fig. 5
#%%

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from functools import reduce

plt.rc('font', family='Helvetica')
path='/Users/bchen/Desktop/Projects/Normal/mutationInfo/'

df=pd.read_csv("/Users/bchen/Desktop/Projects/PKSFiles/Data/C561.table",sep='\t')
# df=df[df['TYPE']=='SNP']            
#   
df['C561_A1_G8']=[float(max(x.split(",")))/(float(max(y.split(",")))+float(max(x.split(",")))) for x,y in zip(df['EPICC_C561_A1_G8_D1.NV'],df['EPICC_C561_A1_G8_D1.NR'])]
df['C561_A1_G9']=[float(max(x.split(",")))/(float(max(y.split(",")))+float(max(x.split(",")))) for x,y in zip(df['EPICC_C561_A1_G9_D1.NV'],df['EPICC_C561_A1_G9_D1.NR'])]
df['C561_B1_G6']=[float(max(x.split(",")))/(float(max(y.split(",")))+float(max(x.split(",")))) for x,y in zip(df['EPICC_C561_B1_G6_D1.NV'],df['EPICC_C561_B1_G6_D1.NR'])]
df['C561_B1_G8']=[float(max(x.split(",")))/(float(max(y.split(",")))+float(max(x.split(",")))) for x,y in zip(df['EPICC_C561_B1_G8_D1.NV'],df['EPICC_C561_B1_G8_D1.NR'])]
df['C561_B1_G9']=[float(max(x.split(",")))/(float(max(y.split(",")))+float(max(x.split(",")))) for x,y in zip(df['EPICC_C561_B1_G9_D1.NV'],df['EPICC_C561_B1_G9_D1.NR'])]
df['C561_C1_G10']=[float(max(x.split(",")))/(float(max(y.split(",")))+float(max(x.split(",")))) for x,y in zip(df['EPICC_C561_C1_G10_D1.NV'],df['EPICC_C561_C1_G10_D1.NR'])]
df['C561_C1_G1']=[float(max(x.split(",")))/(float(max(y.split(",")))+float(max(x.split(",")))) for x,y in zip(df['EPICC_C561_C1_G1_D1.NV'],df['EPICC_C561_C1_G1_D1.NR'])]
df['C561_D1_G10']=[float(max(x.split(",")))/(float(max(y.split(",")))+float(max(x.split(",")))) for x,y in zip(df['EPICC_C561_D1_G10_D1.NV'],df['EPICC_C561_D1_G10_D1.NR'])]
df['C561_D1_G5']=[float(max(x.split(",")))/(float(max(y.split(",")))+float(max(x.split(",")))) for x,y in zip(df['EPICC_C561_D1_G5_D1.NV'],df['EPICC_C561_D1_G5_D1.NR'])]
df['C561_D1_G6']=[float(max(x.split(",")))/(float(max(y.split(",")))+float(max(x.split(",")))) for x,y in zip(df['EPICC_C561_D1_G6_D1.NV'],df['EPICC_C561_D1_G6_D1.NR'])]
df['C561_E1_B1']=[float(max(x.split(",")))/(float(max(y.split(",")))+float(max(x.split(",")))) for x,y in zip(df['EPICC_C561_E1_B1_D1.NV'],df['EPICC_C561_E1_B1_D1.NR'])]
df['C561_E1_G3']=[float(max(x.split(",")))/(float(max(y.split(",")))+float(max(x.split(",")))) for x,y in zip(df['EPICC_C561_E1_G3_D1.NV'],df['EPICC_C561_E1_G3_D1.NR'])]
df['C561_F1_B1']=[float(max(x.split(",")))/(float(max(y.split(",")))+float(max(x.split(",")))) for x,y in zip(df['EPICC_C561_F1_B1_D1.NV'],df['EPICC_C561_F1_B1_D1.NR'])]
df['C561_F1_B2']=[float(max(x.split(",")))/(float(max(y.split(",")))+float(max(x.split(",")))) for x,y in zip(df['EPICC_C561_F1_B2_D1.NV'],df['EPICC_C561_F1_B2_D1.NR'])]
df['C561_G1_B1']=[float(max(x.split(",")))/(float(max(y.split(",")))+float(max(x.split(",")))) for x,y in zip(df['EPICC_C561_G1_B1_D1.NV'],df['EPICC_C561_G1_B1_D1.NR'])]
df['C561_G1_B2']=[float(max(x.split(",")))/(float(max(y.split(",")))+float(max(x.split(",")))) for x,y in zip(df['EPICC_C561_G1_B2_D1.NV'],df['EPICC_C561_G1_B2_D1.NR'])]
df['C561_H1_B1']=[float(max(x.split(",")))/(float(max(y.split(",")))+float(max(x.split(",")))) for x,y in zip(df['EPICC_C561_H1_B1_D1.NV'],df['EPICC_C561_H1_B1_D1.NR'])]
df['C561_H1_B2']=[float(max(x.split(",")))/(float(max(y.split(",")))+float(max(x.split(",")))) for x,y in zip(df['EPICC_C561_H1_B2_D1.NV'],df['EPICC_C561_H1_B2_D1.NR'])]
df['gene']=[x.split("|")[3] for x in df['CSQ']]
df['mut']=[(x.split("|")[10]).split(":")[1]  if len((x.split("|")[10]).split(":"))==2 else ' ' for x in df['CSQ']]
df561=df[[ 'C561_B1_G9','C561_B1_G6','C561_B1_G8', 'C561_C1_G10', 'C561_C1_G1', 
           'C561_D1_G6', 'C561_D1_G5', 'C561_D1_G10','C561_A1_G9','C561_A1_G8',  
           'C561_G1_B2','C561_G1_B1', 'C561_H1_B2', 'C561_H1_B1','C561_F1_B2','C561_F1_B1',
            'C561_E1_G3']]
# 
# df561[df561>0]=1
df561['mm']=[str(x)+"_"+str(y) for x,y in zip(df['gene'],df['mut'])]
df561=df561[df561[df561[df561.columns[:-1]]>0].count(1)<13]
df561.index=range(len(df561))
sample_names = [
 'APC_c.835-8A>G',
 'APC_c.4464del']


# fig=plt.figure(figsize=(15,20))

g=sns.clustermap(df561[df561.columns[:-1]],figsize=(15,20),cmap="Greens", col_cluster=False,annot_kws={"size": 20})
reorder = g.dendrogram_row.reordered_ind
new_positions = [reorder.index(i) for i in [8609,8612]]
plt.setp(g.ax_heatmap.yaxis.set_ticks(new_positions))
plt.setp(g.ax_heatmap.yaxis.set_ticklabels(sample_names))

hm = g.ax_heatmap.get_position()
plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), fontsize=6)
g.ax_heatmap.set_position([hm.x0, hm.y0+0.1, hm.width*0.5, hm.height])
col = g.ax_col_dendrogram.get_position()
g.ax_col_dendrogram.set_position([col.x0, col.y0, col.width*0.25, col.height*0.5])


g.ax_row_dendrogram.set_visible(False)
for xtick in g.ax_heatmap.axes.get_xticklabels():
    if xtick.get_text().split("_")[1].startswith("E"):
        xtick.set_color("#FF7F00")
    elif xtick.get_text().split("_")[1].startswith("A"):
        xtick.set_color("#E41A1C")
    elif xtick.get_text().split("_")[1].startswith("B"):
        xtick.set_color("#377EB8")
    elif xtick.get_text().split("_")[1].startswith("C"):
        xtick.set_color("#4DAF4A")      
    elif xtick.get_text().split("_")[1].startswith("D"):
        xtick.set_color("#984EA3")
    else:
        xtick.set_color("black")


# plt.tight_layout()
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize = 20)
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = 20)



plt.savefig('/Users/bchen/Desktop/Projects/PKSFiles/Figures/FigureCaseshow.pdf')

plt.show()

#########
#####Extended Data Fig.2
dfall=pd.read_csv("/Users/bchen/Desktop/Projects/PKSFiles/Data/EPICC_CancerNormal_normalized_proportions.csv",index_col='Unnamed: 0')


Pnormal={
'C516':['E1_B1'],
'C518':['E1_B1'],
'C519':['A1_G10','A1_G1','A1_G5','B1_G2','B1_G3','D1_G2','D1_G9','E1_G1','E1_B1'],
'C522':['E1_B1'], 
'C524':['A1_G3','A1_G6','E1_B1'],
'C525':['E1_B2'], 
'C527':['A1_G10','A1_G8','C1_G7','C1_G9','D1_B1','E1_B1'],
'C528':['B1_G6','D1_G6','D1_G8','D1_G9'],
'C530':['B1_G2','B1_G7','C1_G10','C1_G7','E1_G1','E1_B1','E1_G2'],
'C531':['D1_G1','D1_G2','D1_G7'],
'C532':[],
'C536':[],  
'C537':['B1_G3','B1_G4','D1_G7','D1_G10','E1_G3'],
'C538':['E1_B1'],
'C539':['E1_B1'],
'C542':['E1_B1'],
'C543':['A1_G2','A1_G10','C1_G1','C1_G8','D1_G9','D1_G10','E1_B1'],
'C544':['B1_G1','B1_G2','C1_G3','E1_G1','E1_G3','E1_B1'],
'C547':['A1_G4','B1_G2','B1_G3','B1_G4','C1_G7','D1_G6','D1_G8','E1_G1'],
'C548':['E1_B1'],
'C549':['A1_G2','A1_G8','C1_G4','C1_G8'],
'C550':['A1_G2','B1_G2','B1_G3','C1_G3','C1_G4','D1_G2','D1_G5'],
'C551':['D1_G5','E1_B1'],
'C552':['A1_G1','B1_G10','B1_G8','B1_G9','C1_G1','C1_G8','E1_G3'],
'C554':[],  
'C555':['B1_G10','C1_G1','C1_G7','C1_G9'],
'C559':[],  
'C560':['A1_G8','D1_G5','D1_G9','E1_B1'],
'C561':['E1_B1','E1_G3'],
'C562':['B1_G3','B1_G4','C1_G6','D1_G2','E1_G2','E1_B1']
}


pids=['C516','C519','C524','C527','C528',
      'C530','C531','C537','C543','C544',
      'C547','C548','C549','C550','C551',
      'C552','C555','C560','C562']

dfisumlist=[]
fig=plt.figure(figsize=(20,12))
axi=1
for pid in pids:
    ax=fig.add_subplot(4,5,axi)
    axi+=1
    dfi=dfall[[x==pid for x in dfall['Pid']]]
    if len(dfi)>0:
        dfi.drop("Pid",axis=1,inplace=True)
        dfi=dfi[['SPS7', 'SPS6', 'SPS5', 'SPS4', 'SPS3', 'SPS2', 'SPS1']]
        sgroup=[]
        for i in range(len(dfi)):
            if "_".join(dfi.index[i].split("_")[1:3]) in Pnormal[pid]:
                if "_".join(dfi.index[i].split("_")[1:3]).startswith("E"):
                    sgroup.append("3distant normal")
                else:
                    sgroup.append("2adjacent normal")
            elif "_".join(dfi.index[i].split("_")[1:3]).startswith("F"):
                sgroup.append("polyps") 
            elif "_".join(dfi.index[i].split("_")[1:3]).startswith("G"):
                sgroup.append("polyps") 
            elif "_".join(dfi.index[i].split("_")[1:3]).startswith("H"):
                sgroup.append("polyps") 
            else:
                sgroup.append("1cancer")
        dfi['group']=sgroup
        dfi[dfi.columns[:-1]] = dfi[dfi.columns[:-1]].apply(pd.to_numeric)
        dfisum=dfi.groupby('group').mean()
        dfisum.index.name=pid
        
        
        dfisum.plot(kind='bar', stacked=True,ax=ax,color=['brown','yellow','orange','purple','green','blue','red'],alpha=0.8)
        ax.set_title(pid,fontsize=10) 
        labels = [item.get_text() for item in ax.get_xticklabels()]
        plt.xlabel("")
        labels2=[x[1:] for x in labels]
        ax.set_xticklabels(labels2)
        if pid =='C562':
            ax.legend(bbox_to_anchor=(1.04,1), loc="upper left", title="Signature")
        else:
            ax.legend("")
    else:
        pass
    dfisum['pid']=pid
    dfisum.index.name=""
    dfisumlist.append(dfisum)

plt.tight_layout()
plt.savefig('/Users/bchen/Desktop/Projects/PKSFiles/Figures/SuppFigure2.pdf')
plt.show()





#####
#####Extended Data Fig. 3
dfdelall=pd.read_csv("/Users/bchen/Desktop/Projects/PKSFiles/Data/EPICC_all_del.csv",index_col='Unnamed: 0')
dfdeldriver=pd.read_csv("/Users/bchen/Desktop/Projects/PKSFiles/Data/EPICC_driver_del.csv",index_col='Unnamed: 0')          



fig=plt.figure(figsize=(8,4))
ax1=fig.add_subplot(1,2,1)

ax1=sns.boxplot(x='MSIMSS',y='length of T-homopolymer',data=dfdelall,palette='Greys')
ax1.text(-0.1, 1.1, "A", transform=ax1.transAxes,fontsize=16, fontweight='bold', va='top', ha='right')
    
plt.title("short T-del of all genes")
plt.ylabel("length of T-homopolymer")
plt.xlabel("EPICC samples")
pairs = [('MSI', 'MSS')] 
add_stat_annotation(ax1, box_pairs=pairs,x='MSIMSS',y='length of T-homopolymer',data=dfdelall,test='Mann-Whitney', text_format='star', loc='inside', verbose=2)



ax2=fig.add_subplot(1,2,2)
ax2=sns.boxplot(x='MSIMSS',y='rep',data=dfdeldriver,palette='Greys')
plt.title("short T-del of driver genes")
plt.ylabel("length of T-homopolymer")
plt.xlabel("EPICC samples")
plt.tight_layout()
add_stat_annotation(ax2, box_pairs=pairs,x='MSIMSS',y='rep',data=dfdeldriver,test='Mann-Whitney', text_format='star', loc='inside', verbose=2)
ax1.text(-0.1, 1.1, "B", transform=ax2.transAxes,fontsize=16, fontweight='bold', va='top', ha='right')
 
plt.savefig("/Users/bchen/Desktop/Projects/PKSFiles/Figures/ExtendedFigure3.pdf")
plt.show()
print("going on")


#####
#####Extended Data Fig. 4

df4=pd.read_csv("/Users/bchen/Desktop/Projects/PKSFiles/Data/EPICCpksReads.csv",index_col='Unnamed: 0')
pids=['C516','C518','C519','C524','C527',
      'C528','C530','C531','C536','C537',
      'C538','C539','C542','C543','C544',
      'C547','C548','C549','C550','C551',
      'C552','C555','C560','C561','C562']
with PdfPages('/Users/bchen/Desktop/Projects/PKSFiles/Figures/FigureS4.pdf') as pdf:
    fig=plt.figure(figsize=(25,10))
    axi=1
    for pid in pids:
        dfi=df4[[x.split("_")[0]==pid for x in df4.index]]
        dfi=dfi[["Z1" not in x for x in dfi.index]]  
        dfi2=dfi[dfi.columns[2:]].fillna(0)
        
        dfi2[dfi2>0]=1
        if dfi2.sum().sum()>0:
            ax=fig.add_subplot(3,6,axi)
            axi+=1
            plt.title("%s"%pid,fontsize=18)

            dfiorder=[]
            for i in range(len(dfi2)):
                if "_".join(dfi2.index[i].split("_")[1:3]) in Pnormal[pid]:
                    if dfi2.index[i].split("_")[1][0]=='E':
                        dfiorder.append(1)
                    else:
                        dfiorder.append(2)
                else:   
                    if dfi2.index[i].split("_")[1][0]=='A':
                        dfiorder.append(3)
                    elif dfi2.index[i].split("_")[1][0]=='B':
                        dfiorder.append(4)
                    elif dfi2.index[i].split("_")[1][0]=='C':
                        dfiorder.append(5)
                    elif dfi2.index[i].split("_")[1][0]=='D':
                        dfiorder.append(6)               
                    elif dfi2.index[i].split("_")[1][0]=='F':
                        dfiorder.append(7)            
                    elif dfi2.index[i].split("_")[1][0]=='G':
                        dfiorder.append(8)         
                    elif dfi2.index[i].split("_")[1][0]=='H':
                        dfiorder.append(9)
                    else:
                        dfiorder.append(10)
                        
            dfi2['order']=dfiorder
            dfi2.sort_values(by='order',inplace=True)
            dfi2.drop("order",axis=1,inplace=True)
            ax=sns.heatmap(dfi2,cmap="gist_gray_r",linewidths=.5,linecolor='lightgrey',cbar=False,square=True)
            for axis in ['top','bottom','left','right']:
                ax.spines[axis].set_visible(True)
                ax.spines[axis].set_color('black')
                ax.spines[axis].set_visible(True)
                ax.spines[axis].set_color('black')
            plt.xlabel("PKS+ gene")
            
                
            if pid in Pnormal.keys():
            
            
                labels = [item.get_text() for item in ax.get_yticklabels()]
                labels = ["_".join(x.split("_")[1:]) for x in labels]
            

                labels2=[]
                for x in labels:
                    if "_".join(x.split("_")[0:2]) in Pnormal[pid]:
                        if "_".join(x.split("_")[0:2]).startswith("E"):
                            #labels2.append("E_"+x)
                            labels2.append(x)
                        else:
                            #labels2.append("W_"+x)
                            labels2.append(x)
                    else:
                        labels2.append(x)
                

                ax.set_yticklabels(labels2)
        
            else:
                labels = [item.get_text() for item in ax.get_xticklabels()]
                labels = ["_".join(x.split("_")[1:]) for x in labels]
                ax.set_yticklabels(labels)
        

##
            for xtick in ax.get_yticklabels():
    
                if pid in Pnormal.keys():
    
                    if "_".join(xtick.get_text().split("_")[0:2]) in Pnormal[pid]:
                        if "_".join(xtick.get_text().split("_")[0:2]).startswith("E"):
                            xtick.set_color("#FF7F00")
                        else:
                            xtick.set_color("#00CDCD")
                    elif xtick.get_text().split("_")[0].startswith("A"):
                        xtick.set_color("#E41A1C")
                    elif xtick.get_text().split("_")[0].startswith("B"):
                        xtick.set_color("#377EB8")
                    elif xtick.get_text().split("_")[0].startswith("C"):
                        xtick.set_color("#4DAF4A")      
                    elif xtick.get_text().split("_")[0].startswith("D"):
                        xtick.set_color("#984EA3")
                    else:
                        xtick.set_color("black")
                elif xtick.get_text().split("_")[0].startswith("E"):

                    xtick.set_color("#FF7F00")
                elif xtick.get_text().split("_")[0].startswith("A"):

                    xtick.set_color("#E41A1C")
                elif xtick.get_text().split("_")[0].startswith("B"):

                    xtick.set_color("#377EB8")
                elif xtick.get_text().split("_")[0].startswith("C"):

                    xtick.set_color("#4DAF4A")      
                elif xtick.get_text().split("_")[0].startswith("D"):

                    xtick.set_color("#984EA3")
                else:

                    xtick.set_color("black")

##    
    plt.tight_layout()
    pdf.savefig()
    plt.close()


#####
#####Extended Data Fig. 5

dfNN3=pd.read_csv("/Users/bchen/Desktop/Projects/PKSFiles/Data/percentageofpksReadsEPICC.csv")
fig=plt.figure(figsize=(10,6))
ax = sns.barplot(x='pid', y='percentage2',hue='type',data=dfNN3,palette=["#E41A1C","#00CDCD","#FF7F00"],hue_order=['cancer', 'adjacent normal', 'distant normal'])
plt.ylabel("percentage of pks+ ecoli reads")
ax.set_yscale('log')
plt.xticks(rotation=90)
plt.savefig('/Users/bchen/Desktop/Projects/PKSFiles/Figures/FigureS5.pdf')
plt.show()

#####
#####Extended Data Fig.6 

dfTdelPksMSS=pd.read_csv("/Users/bchen/Desktop/Projects/PKSFiles/Data/dfdelpksMSS.csv",index_col='Unnamed: 0')
dfTdelPksMSI=pd.read_csv("/Users/bchen/Desktop/Projects/PKSFiles/Data/dfdelpksMSI.csv",index_col='Unnamed: 0')
p=sns.regplot(x=dfTdelPksMSS['pks'],y=dfTdelPksMSS['t-del'], fit_reg=True, marker="o", color="black", scatter_kws={'s':50})
for line in range(len(dfTdelPksMSS)):
     plt.text(dfTdelPksMSS['pks'][line]+0.001, dfTdelPksMSS['t-del'][line], dfTdelPksMSS.index[line], horizontalalignment='left', size='small', color='black')
sns.scatterplot(x=[0.9],y=[0.6],marker="o", color="black",s=60)
plt.text(0.93, 0.59, "MSS", horizontalalignment='left', size='medium', color='black')

# 
sns.regplot(x=dfTdelPksMSI['pks'],y=dfTdelPksMSI['t-del'], fit_reg=False, marker="o", color="grey", scatter_kws={'s':50})
for line in range(len(dfTdelPksMSI)):
     plt.text(dfTdelPksMSI['pks'][line]+0.001, dfTdelPksMSI['t-del'][line], dfTdelPksMSI.index[line], horizontalalignment='left', size='small', color='black')
sns.scatterplot(x=[0.9],y=[0.57],marker="o", color="grey",s=60)
plt.text(0.93, 0.56, "MSI", horizontalalignment='left', size='medium', color='black')


est = sm.OLS(dfTdelPksMSS['t-del'],dfTdelPksMSS['pks'])
est2 = est.fit()
print("summary()\n",est2.summary())
plt.text(0.005,0.6,"pvalues:%s"%float(est2.pvalues))
plt.text(0.005,0.57,"tvalues:%.2f"%float(est2.tvalues))
plt.text(0.005,0.54,"rsquared:%.2f"%float(est2.rsquared))
plt.text(0.005,0.51,"rsquared_adj:%.2f"%float(est2.rsquared_adj))

plt.ylabel("proportion of short T-dels")
plt.xlabel("proportion of samples with pks+ genes")
plt.savefig('/Users/bchen/Desktop/Projects/PKSFiles/Figures/FigureS6.pdf')
plt.show()




