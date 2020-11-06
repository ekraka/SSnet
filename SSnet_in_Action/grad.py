from math import log10, floor
from keras.models import Sequential,Model
import tensorflow as tf
from keras.layers import Dense, Activation, Input,RepeatVector,Embedding, Flatten, Concatenate,Dropout
from keras.models import Model
from keras.utils.vis_utils import model_to_dot
#from sklearn import metrics as mt
from keras import metrics
#from sklearn.metrics import confusion_matrix
#from sklearn.metrics import accuracy_score
#import matplotlib.pyplot as plt
#from sklearn.model_selection import train_test_split
from keras.models import Sequential
from keras.layers import average, concatenate,RepeatVector,Lambda,add,subtract
from keras import backend as K
from keras.layers import Conv2D, MaxPooling2D,Conv1D,GlobalMaxPooling1D,MaxPooling1D,Reshape,Add
from keras.layers import Conv1D, GlobalAveragePooling1D, MaxPooling1D,AveragePooling1D
from keras.regularizers import l2 
import numpy as np

def get_heatm(conv, x, model, ref = 0):

    predicted_class = int(round(model.predict(x)[0][0]))

    predicted_class_output = model.output#[:, predicted_class]

    last_conv_layer = model.get_layer(conv)

    # This is the gradient of the predicted class with regard to
    # the output feature map of `block5_conv3`
    grads = K.gradients(predicted_class_output, last_conv_layer.output)[0]
    
    #print (grads.shape)
    a, s = grads.shape[-2:]

    # This is a vector of shape (512,), where each entry
    # is the mean intensity of the gradient over a specific feature map channel
    ##pooled_grads = grads
    pooled_grads = K.mean(grads, axis=(0,1))
    
    #print (pooled_grads.shape)

    # This function allows us to access the values of the quantities we just defined:
    # `pooled_grads` and the output feature map of `block5_conv3`,
    # given a sample image
    #print (last_conv_layer.output.shape)
    iterate = K.function([model.input[0], model.input[1]], [pooled_grads, last_conv_layer.output[0]])

    # These are the values of these two quantities, as Numpy arrays,
    # given our sample image

    #print (iterate(x))
    pooled_grads_value, conv_layer_output_value = iterate(x)
    


    # We multiply each channel in the feature map array
    # by "how important this channel is" with regard to the predicted class
    for i in range(s):
        conv_layer_output_value[:,i] *= pooled_grads_value[i]

    # The channel-wise mean of the resulting feature map
    # is our heatmap of class activation
    #print (conv_layer_output_value.shape)
    heatmap = np.mean(conv_layer_output_value, axis=(-1,))
    heatmap = np.maximum(heatmap, 0)
    if np.max(heatmap) > 0:
        heatmap /= np.max(heatmap)
    return heatmap

def shape_heatmap(heatmap):
    f = int(9000 / heatmap.shape[0])
    print ('Factor:', f)
    l = []
    for i in heatmap:
        l += [i]*f
    if len(l) < 9000:
        l += [0]*(9000 - len(l))
    return np.array(l).reshape((6, 1500))


def round_sig(x, sig=3):
    if x < 1e-2:
        return 0.0
    return round(x, sig-int(floor(log10(abs(x))))-1)

def create_pdb(pdb, heatmap):
    gami, gama = np.amin(heatmap), np.amax(heatmap)
    for i in range (len(heatmap)):
        a = heatmap[i]
        #gami, gama = a.min(), a.max()
        heatmap[i] = np.interp(a, (gami, gama), (0.0, 100.0))
        
        #print (np.max(heatmap))
        
    f = open(pdb, 'r')
    lines = f.readlines()
    f.close()
    
    st, ref = [], []
    mod = 0 # model
    res = [] # residue
    for line in lines:
        k = line.strip().split()
        if len(k) < 1:
            continue
        if line.strip().split()[0] in ['MODEL', 'ENDMDL']:
            st.append(line.strip())
            ref.append(None)
        if 'TER' in line[:3]:
            st.append(line.strip())
            ref.append(None)
            mod += 1
            res = []
        if 'ATOM' in line[:4]:
            nres = int(line[22:26])
            if nres not in res:
                res.append(nres)
            st.append(line[:61])
            ref.append(round_sig(heatmap[mod][len(res) - 1]))
    print ('Number of chains:', mod)
    k_ref = []
    for i in ref:
        if i is not None:
            k_ref.append(i)
        else:
            #print (ref)
            k_ref.append(0)#ref[0])
            
    k_ref = np.array(k_ref)
    k_ref = np.interp(k_ref, (k_ref.min(), k_ref.max()), (0.0, 100.0))
    
    fst = ''
    
    for i in range (len(st)):
        if ref[i] is not None:
            ans = str(round_sig(k_ref[i]))
            #print (ans)
            
            fst += st[i][:-1]+' '*(6-len(ans))+ans + '\n' 
            #print (fst)
            #stop()
        else:
            fst += st[i]+ '\n'
        #print (ref[i])
    '''
    clusters = m_box.job(fst, 1)
    count = i
    for i in clusters:
        print (i)
        l = ['HETATM', str(count), 'Zn', 'Zn', 'A', str(count), round(i[0]), round(i[1]), round(i[2]), '1.00', '0.00', 'Zn']
        fst += "{:>6}{:>5}{:>4} {:>4}{:>2}{:>4}{:>12}{:>8}{:>8}{:>6}{:>6}{:>12}".format(*l)+'\n'
        count += 1
    '''
    g = open(pdb[:-4] + '_GRAD.pdb', 'w')
    #print (fst)
    g.write(fst)
    g.close()

    print ('File written as', pdb[:-4] + '_GRAD.pdb')

def job(proteins, X, pdb, model):
    l = ['conv1d_3', 'conv1d_6', 'conv1d_9', 'conv1d_12']
    heatmap = None
    x = [proteins.reshape((-1, 9000, 2)), X.reshape((-1, 512))]
    for i in l:
            if heatmap is None:
                    heatmap = get_heatm(i, x, model, 1)/len(l)
            else:
                    heatmap += get_heatm(i, x, model)/len(l)

    heatmap = shape_heatmap(heatmap)
    create_pdb(pdb, heatmap)
