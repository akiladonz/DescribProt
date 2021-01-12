#!/home/pyenv/versions/py3.7/bin/python
# coding: utf-8

# In[1]:

import os
import numpy as np
import pandas as pd
import csv
import sys

filename = sys.argv[1]


Lines= []


path = "/var/www/servers/DESCRIBEPROT/search_files"
compilename = os.path.join(path, filename)
f = open(compilename)
for line in f:
    data_line = line.rstrip().split('\t')
    Lines.append(data_line)

#print Lines
# In[2]:


ID=Lines[0][0]
Sequence=list(Lines[1][0])
Vsl_Binary=list(map(int,list(Lines[2][0])))
Vsl_Score=list(map(float,Lines[3][0].rstrip().split(',')))
PsiPred_Binary=list(map(int,list(Lines[4][0])))
Helix_Score=list(map(float,Lines[5][0].rstrip().split(',')))
Strand_Score=list(map(float,Lines[6][0].rstrip().split(',')))
Coil_Score=list(map(float,Lines[7][0].rstrip().split(',')))
mmseq_Binary=list(map(int,list(Lines[8][0])))
mmseq_Score=list(map(float,Lines[9][0].rstrip().split(',')))
SignalP_Binary=list(map(int,list(Lines[10][0])))
SignalP_Score=list(map(float,Lines[11][0].rstrip().split(',')))
DFL_Binary=list(map(int,list(Lines[12][0])))
DFL_Score=list(map(float,Lines[13][0].rstrip().split(',')))
ASA_Binary=list(map(int,list(Lines[14][0])))
ASA_Score=list(map(float,Lines[15][0].rstrip().split(',')))
DisoRNA_Binary=list(map(int,list(Lines[16][0])))
DisoRNA_Score=list(map(float,Lines[17][0].rstrip().split(',')))
DrnaRNA_Binary=list(map(int,list(Lines[18][0])))
DrnaRNA_Score=list(map(float,Lines[19][0].rstrip().split(',')))
DisoDNA_Binary=list(map(int,list(Lines[20][0])))
DisoDNA_Score=list(map(float,Lines[21][0].rstrip().split(',')))
DrnaDNA_Binary=list(map(int,list(Lines[22][0])))
DrnaDNA_Score=list(map(float,Lines[23][0].rstrip().split(',')))
DisoProt_Binary=list(map(int,list(Lines[24][0])))
DisoProt_Score=list(map(float,Lines[25][0].rstrip().split(',')))
Scribber_Binary=list(map(int,list(Lines[26][0])))
Scribber_Score=list(map(float,Lines[27][0].rstrip().split(',')))
Morf_Binary=list(map(int,list(Lines[28][0])))
Morf_Score=list(map(float,Lines[29][0].rstrip().split(',')))
#print (np.array(Morf_Score)).dtype

# In[3]:


from operator import itemgetter
from itertools import *
def BinaryCoordiantes(Binary_Prediction):
              #######################################################################################################
              #Getting the positions of predicted Disorder regions
              Disorder_Index= [] 
              b = 0
              for b in range(0, len(Binary_Prediction), 1):
                                                       if Binary_Prediction[b] == 1:Disorder_Index.append(b)
                                                       
              Disorder_Index=np.unique(Disorder_Index)
              ######################################################################################################
              groups = []
              for k, g in groupby(enumerate(Disorder_Index), lambda x: x[0]-x[1]):
                    groups.append(list(map(itemgetter(1), g)))
              Coordinates_x=[]
              Coordinates_y=[]
              for group in groups:
                    Coordinates_x.append(min(group))
                    Coordinates_y.append(max(group))
              if len(Coordinates_x)==0:
                         Coordinates_x.append(len(Binary_Prediction)+100)
                         Coordinates_y.append(len(Binary_Prediction)+100)
              return Coordinates_x,Coordinates_y


# In[4]:


def HelixCoordiantes(Binary_Prediction):
              #######################################################################################################
              #Getting the positions of predicted Disorder regions
              Disorder_Index= [] 
              b = 0
              for b in range(0, len(Binary_Prediction), 1):
                                                       if Binary_Prediction[b] == 0:Disorder_Index.append(b)
                                                       
              Disorder_Index=np.unique(Disorder_Index)
              ######################################################################################################
              groups = []
              for k, g in groupby(enumerate(Disorder_Index), lambda x: x[0]-x[1]):
                    groups.append(list(map(itemgetter(1), g)))
              Coordinates_x=[]
              Coordinates_y=[]
              for group in groups:
                    Coordinates_x.append(min(group))
                    Coordinates_y.append(max(group))
              return Coordinates_x,Coordinates_y


# In[5]:


def CoilCoordiantes(Binary_Prediction):
              #######################################################################################################
              #Getting the positions of predicted Disorder regions
              Disorder_Index= [] 
              b = 0
              for b in range(0, len(Binary_Prediction), 1):
                                                       if Binary_Prediction[b] == 2:Disorder_Index.append(b)
                                                       
              Disorder_Index=np.unique(Disorder_Index)
              ######################################################################################################
              groups = []
              for k, g in groupby(enumerate(Disorder_Index), lambda x: x[0]-x[1]):
                    groups.append(list(map(itemgetter(1), g)))
              Coordinates_x=[]
              Coordinates_y=[]
              for group in groups:
                    Coordinates_x.append(min(group))
                    Coordinates_y.append(max(group))
              return Coordinates_x,Coordinates_y


# In[6]:


def Conv4(Binary_Prediction):
              #######################################################################################################
              #Getting the positions of predicted Disorder regions
              Disorder_Index= [] 
              b = 0
              for b in range(0, len(Binary_Prediction), 1):
                                                       if Binary_Prediction[b] == 3:Disorder_Index.append(b)
                                                       
              Disorder_Index=np.unique(Disorder_Index)
              ######################################################################################################
              groups = []
              for k, g in groupby(enumerate(Disorder_Index), lambda x: x[0]-x[1]):
                    groups.append(list(map(itemgetter(1), g)))
              Coordinates_x=[]
              Coordinates_y=[]
              for group in groups:
                    Coordinates_x.append(min(group))
                    Coordinates_y.append(max(group))
              return Coordinates_x,Coordinates_y


# In[7]:


def Conv5(Binary_Prediction):
              #######################################################################################################
              #Getting the positions of predicted Disorder regions
              Disorder_Index= [] 
              b = 0
              for b in range(0, len(Binary_Prediction), 1):
                                                       if Binary_Prediction[b] == 4:Disorder_Index.append(b)
                                                       
              Disorder_Index=np.unique(Disorder_Index)
              ######################################################################################################
              groups = []
              for k, g in groupby(enumerate(Disorder_Index), lambda x: x[0]-x[1]):
                    groups.append(list(map(itemgetter(1), g)))
              Coordinates_x=[]
              Coordinates_y=[]
              for group in groups:
                    Coordinates_x.append(min(group))
                    Coordinates_y.append(max(group))
              return Coordinates_x,Coordinates_y


# In[8]:


def Conv6(Binary_Prediction):
              #######################################################################################################
              #Getting the positions of predicted Disorder regions
              Disorder_Index= [] 
              b = 0
              for b in range(0, len(Binary_Prediction), 1):
                                                       if Binary_Prediction[b] == 5:Disorder_Index.append(b)
                                                       
              Disorder_Index=np.unique(Disorder_Index)
              ######################################################################################################
              groups = []
              for k, g in groupby(enumerate(Disorder_Index), lambda x: x[0]-x[1]):
                    groups.append(list(map(itemgetter(1), g)))
              Coordinates_x=[]
              Coordinates_y=[]
              for group in groups:
                    Coordinates_x.append(min(group))
                    Coordinates_y.append(max(group))
              return Coordinates_x,Coordinates_y


# In[9]:


def Conv7(Binary_Prediction):
              #######################################################################################################
              #Getting the positions of predicted Disorder regions
              Disorder_Index= [] 
              b = 0
              for b in range(0, len(Binary_Prediction), 1):
                                                       if Binary_Prediction[b] == 6:Disorder_Index.append(b)
                                                       
              Disorder_Index=np.unique(Disorder_Index)
              ######################################################################################################
              groups = []
              for k, g in groupby(enumerate(Disorder_Index), lambda x: x[0]-x[1]):
                    groups.append(list(map(itemgetter(1), g)))
              Coordinates_x=[]
              Coordinates_y=[]
              for group in groups:
                    Coordinates_x.append(min(group))
                    Coordinates_y.append(max(group))
              return Coordinates_x,Coordinates_y


# In[10]:


def Conv8(Binary_Prediction):
              #######################################################################################################
              #Getting the positions of predicted Disorder regions
              Disorder_Index= [] 
              b = 0
              for b in range(0, len(Binary_Prediction), 1):
                                                       if Binary_Prediction[b] == 7:Disorder_Index.append(b)
                                                       
              Disorder_Index=np.unique(Disorder_Index)
              ######################################################################################################
              groups = []
              for k, g in groupby(enumerate(Disorder_Index), lambda x: x[0]-x[1]):
                    groups.append(list(map(itemgetter(1), g)))
              Coordinates_x=[]
              Coordinates_y=[]
              for group in groups:
                    Coordinates_x.append(min(group))
                    Coordinates_y.append(max(group))
              return Coordinates_x,Coordinates_y


# In[11]:


def Conv9(Binary_Prediction):
              #######################################################################################################
              #Getting the positions of predicted Disorder regions
              Disorder_Index= [] 
              b = 0
              for b in range(0, len(Binary_Prediction), 1):
                                                       if Binary_Prediction[b] == 8:Disorder_Index.append(b)
                                                       
              Disorder_Index=np.unique(Disorder_Index)
              ######################################################################################################
              groups = []
              for k, g in groupby(enumerate(Disorder_Index), lambda x: x[0]-x[1]):
                    groups.append(list(map(itemgetter(1), g)))
              Coordinates_x=[]
              Coordinates_y=[]
              for group in groups:
                    Coordinates_x.append(min(group))
                    Coordinates_y.append(max(group))
              return Coordinates_x,Coordinates_y


# In[12]:


def Conv10(Binary_Prediction):
              #######################################################################################################
              #Getting the positions of predicted Disorder regions
              Disorder_Index= [] 
              b = 0
              for b in range(0, len(Binary_Prediction), 1):
                                                       if Binary_Prediction[b] == 9:Disorder_Index.append(b)
                                                       
              Disorder_Index=np.unique(Disorder_Index)
              ######################################################################################################
              groups = []
              for k, g in groupby(enumerate(Disorder_Index), lambda x: x[0]-x[1]):
                    groups.append(list(map(itemgetter(1), g)))
              Coordinates_x=[]
              Coordinates_y=[]
              for group in groups:
                    Coordinates_x.append(min(group))
                    Coordinates_y.append(max(group))
              return Coordinates_x,Coordinates_y


# In[13]:


for b in range(70,len(ASA_Score),1):
              SignalP_Score.append(0.0)
              SignalP_Binary.append(0)
# In[14]:


PsiPred_Score=[]
for b in range(0,len(Helix_Score),1):
                PsiPred_Score.append(max([Helix_Score[b],Strand_Score[b],Coil_Score[b]]))


# In[15]:




# In[16]:


Vsl_Corrdinates_x,Vsl_Corrdinates_y=BinaryCoordiantes(Vsl_Binary)
SignalP_Corrdinates_x,SignalP_Corrdinates_y=BinaryCoordiantes(SignalP_Binary)
DFL_Corrdinates_x,DFL_Corrdinates_y=BinaryCoordiantes(DFL_Binary)
ASA_Corrdinates_x,ASA_Corrdinates_y=BinaryCoordiantes(ASA_Binary)
DisoRNA_Corrdinates_x,DisoRNA_Corrdinates_y=BinaryCoordiantes(DisoRNA_Binary)
DisoDNA_Corrdinates_x,DisoDNA_Corrdinates_y=BinaryCoordiantes(DisoDNA_Binary)
DisoProt_Corrdinates_x,DisoProt_Corrdinates_y=BinaryCoordiantes(DisoProt_Binary)
Morf_Corrdinates_x,Morf_Corrdinates_y=BinaryCoordiantes(Morf_Binary)
DrnaDNA_Corrdinates_x,DrnaDNA_Corrdinates_y=BinaryCoordiantes(DrnaDNA_Binary)
DrnaRNA_Corrdinates_x,DrnaRNA_Corrdinates_y=BinaryCoordiantes(DrnaRNA_Binary)
Scribber_Corrdinates_x,Scribber_Corrdinates_y=BinaryCoordiantes(Scribber_Binary)


# In[17]:


Helix_Corrdinates_x,Helix_Corrdinates_y=HelixCoordiantes(PsiPred_Binary)
Strand_Corrdinates_x,Strand_Corrdinates_y=BinaryCoordiantes(PsiPred_Binary)
Coil_Corrdinates_x,Coil_Corrdinates_y=CoilCoordiantes(PsiPred_Binary)


# In[18]:


Conv1_Corrdinates_x,Conv1_Corrdinates_y=HelixCoordiantes(mmseq_Binary)
Conv2_Corrdinates_x,Conv2_Corrdinates_y=BinaryCoordiantes(mmseq_Binary)
Conv3_Corrdinates_x,Conv3_Corrdinates_y=CoilCoordiantes(mmseq_Binary)
Conv4_Corrdinates_x,Conv4_Corrdinates_y=Conv4(mmseq_Binary)
Conv5_Corrdinates_x,Conv5_Corrdinates_y=Conv5(mmseq_Binary)
Conv6_Corrdinates_x,Conv6_Corrdinates_y=Conv6(mmseq_Binary)
Conv7_Corrdinates_x,Conv7_Corrdinates_y=Conv7(mmseq_Binary)
Conv8_Corrdinates_x,Conv8_Corrdinates_y=Conv8(mmseq_Binary)
Conv9_Corrdinates_x,Conv9_Corrdinates_y=Conv9(mmseq_Binary)
Conv10_Corrdinates_x,Conv10_Corrdinates_y=Conv10(mmseq_Binary)





from sklearn.preprocessing import MinMaxScaler
def Resccaling(Score):
            scaler = MinMaxScaler(feature_range=(0.33,1))
            Score=np.asarray(Score)
            Score=Score.reshape(-1, 1)
            data=[0,1]
            data=np.asarray(data)
            data=data.reshape(-1, 1)
 #           scaler.fit(data)
            Score=scaler.fit_transform(Score)
            Scaled_Score=np.concatenate(Score)
            return Scaled_Score
			
			
def Resccaling1(Score):
            scaler = MinMaxScaler(feature_range=(0,1))
            Score=np.asarray(Score)
            Score=Score.reshape(-1, 1)
            data=[0,1]
            data=np.asarray(data)
            data=data.reshape(-1, 1)
 #           scaler.fit(data)
            Score=scaler.fit_transform(Score)
            Scaled_Score=np.concatenate(Score)
            return Scaled_Score

PsiPred_Score=list(Resccaling(PsiPred_Score))
mmseq_Score=list(Resccaling1(mmseq_Score))
# In[19]:


SS_Code=[]

for b in range(0,len(PsiPred_Binary),1):
    if PsiPred_Binary[b]==0:SS_Code.append('Helix')
    elif PsiPred_Binary[b]==1:SS_Code.append('Strand')
    elif PsiPred_Binary[b]==2:SS_Code.append('Coil')
    else:SS_Code.append('X')
    

x = str(Sequence)
x = x[1:-1]
x = x.replace(",", "")
x = x.replace(" ", "")
x = x.replace("''", "")
x = x.replace("'", "")


#print (np.array(Morf_Score)).dtype

print x
print Conv1_Corrdinates_x
print Conv1_Corrdinates_y
print Conv2_Corrdinates_x
print Conv2_Corrdinates_y
print Conv3_Corrdinates_x
print Conv3_Corrdinates_y
print Conv4_Corrdinates_x
print Conv4_Corrdinates_y
print Conv5_Corrdinates_x
print Conv5_Corrdinates_y
print Conv6_Corrdinates_x
print Conv6_Corrdinates_y
print Conv7_Corrdinates_x
print Conv7_Corrdinates_y
print Conv8_Corrdinates_x
print Conv8_Corrdinates_y
print Conv9_Corrdinates_x
print Conv9_Corrdinates_y
print Conv10_Corrdinates_x
print Conv10_Corrdinates_y
print mmseq_Score
print ASA_Corrdinates_x
print ASA_Corrdinates_y
print ASA_Score
print Helix_Corrdinates_x
print Helix_Corrdinates_y
print Strand_Corrdinates_x
print Strand_Corrdinates_y
print Coil_Corrdinates_x
print Coil_Corrdinates_y
print PsiPred_Score
print SignalP_Corrdinates_x
print SignalP_Corrdinates_y
print SignalP_Score
print Vsl_Corrdinates_x
print Vsl_Corrdinates_y
print Vsl_Score
print DFL_Corrdinates_x
print DFL_Corrdinates_y
print DFL_Score
print DisoProt_Corrdinates_x
print DisoProt_Corrdinates_y
print DisoProt_Score
print Scribber_Corrdinates_x
print Scribber_Corrdinates_y
print Scribber_Score
print Morf_Corrdinates_x
print Morf_Corrdinates_y
print Morf_Score
print DisoDNA_Corrdinates_x
print DisoDNA_Corrdinates_y
print DisoDNA_Score
print DrnaDNA_Corrdinates_x
print DrnaDNA_Corrdinates_y
print DrnaDNA_Score
print DisoRNA_Corrdinates_x
print DisoRNA_Corrdinates_y
print DisoRNA_Score
print DrnaRNA_Corrdinates_x
print DrnaRNA_Corrdinates_y
print DrnaRNA_Score
print SS_Code





