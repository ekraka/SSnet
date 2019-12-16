import random
import pickle
import sys
import timeit

import numpy as np
from sklearn.metrics import roc_curve

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim

from sklearn.metrics import roc_auc_score, precision_score, recall_score


class CompoundProteinInteractionPrediction(nn.Module):
    def __init__(self):
        super(CompoundProteinInteractionPrediction, self).__init__()
        self.embed_fingerprint = nn.Embedding(n_fingerprint, dim)
        #proteins = []
        #f = np.load('../../../final_data_dude_rd2_opp.npy').item()
        #self.proteins = torch.from_numpy(np.array(f['proteins']).reshape((-1,)))
        self.w1 = self.branch_conv(5, 2, 28)
        self.w11 = self.branch_conv(5, 28, 64)
        self.w111 = self.branch_conv(5, 64, 128)
        self.w2 = self.branch_conv(10, 2, 28)
        self.w22 = self.branch_conv(10, 28, 64)
        self.w222 = self.branch_conv(10, 64, 128)
        self.w3 = self.branch_conv(15, 2, 28)
        self.w33 = self.branch_conv(15, 28, 64)
        self.w333 = self.branch_conv(15, 64, 128)
        self.w4 = self.branch_conv(20, 2, 28)
        self.w44 = self.branch_conv(20, 28, 64)
        self.w444 = self.branch_conv(20, 64, 128)
        self.embed_word = nn.Linear(512*2, dim)
        #self.bc3 = self.branch_conv(15)
        #print (self.bc1.shape, self.bc2.shape, self.bc3.shape)
        self.W_gnn = nn.ModuleList([nn.Linear(dim, dim)
                                    for _ in range(layer_gnn)])
        self.W_cnn = nn.ModuleList([nn.Conv2d(
                     in_channels=1, out_channels=1, kernel_size=2*window+1,
                     stride=1, padding=window) for _ in range(layer_cnn)])
        self.W_attention = nn.Linear(dim, dim)
        self.W_out = nn.ModuleList([nn.Linear(2*dim, 2*dim)
                                    for _ in range(layer_output)])
        self.W_interaction = nn.Linear(2*dim, 2)

    def branch_conv(self, kernel_size = 5, in_channel = 1, out_channel = 1):
        return nn.Conv1d(
                     in_channels=in_channel, out_channels=out_channel, kernel_size=kernel_size,
                     stride=1).cuda() # remove cuda when using CPU !!

    def gnn(self, xs, A, layer):
        for i in range(layer):
            hs = torch.relu(self.W_gnn[i](xs))
            xs = xs + torch.matmul(A, hs)
        # return torch.unsqueeze(torch.sum(xs, 0), 0)
        return torch.unsqueeze(torch.mean(xs, 0), 0)

    def attention_cnn(self, x, xs, layer):
        """The attention mechanism is applied to the last layer of CNN."""

        xs = torch.unsqueeze(torch.unsqueeze(xs, 0), 0)
        for i in range(layer):
            xs = torch.relu(self.W_cnn[i](xs))
        xs = torch.squeeze(torch.squeeze(xs, 0), 0)

        h = torch.relu(self.W_attention(x))
        hs = torch.relu(self.W_attention(xs))
        weights = torch.tanh(F.linear(h, hs))
        ys = torch.t(weights) * hs

        # return torch.unsqueeze(torch.sum(ys, 0), 0)
        return torch.unsqueeze(torch.mean(ys, 0), 0)

    def forward(self, inputs):

        fingerprints, adjacency, words = inputs

        """Compound vector with GNN."""
        fingerprint_vectors = self.embed_fingerprint(fingerprints)
        compound_vector = self.gnn(fingerprint_vectors, adjacency, layer_gnn)
        #print (compound_vector.shape)

        """Protein vector with attention-CNN."""
        #word_vectors = self.embed_word(words)
        #print (words.shape)
        words = words.reshape((-1, 2, 9000))
        #words = torch.tensor(words).to(device)
        w1 = self.w1(words)
        w1 = self.w11(w1)
        w1 = self.w111(w1)
        #w1 = self.bc1(w1)
        w_1 = F.avg_pool1d(w1, kernel_size=w1.size()[-1]).view(-1,)
        #print (w1.shape)
        u = w_1.repeat(w1.shape[-1]).view(-1, w1.shape[-1])
        wu_vec = (u-w1)**2
        wu = F.avg_pool1d(wu_vec, kernel_size=wu_vec.size()[-1]).view(-1,) 
        w1 = torch.cat((w_1, wu), 0).view(-1, 1)
        #print (w1.shape)
        #stop()

        #words = words.reshape((-1, 1, 18000))
        w2 = self.w2(words)
        w2 = self.w22(w2)
        w2 = self.w222(w2)
        #w1 = self.bc1(w1)
        w_2 = F.avg_pool1d(w2, kernel_size=w2.size()[-1]).view(-1,)
        #print (w1.shape)
        u = w_2.repeat(w2.shape[-1]).view(-1, w2.shape[-1])
        wu_vec = (u-w2)**2
        wu = F.avg_pool1d(wu_vec, kernel_size=wu_vec.size()[-1]).view(-1,)
        w2 = torch.cat((w_2, wu), 0).view(-1, 1)

        #w2 = F.avg_pool1d(w2, kernel_size=w2.size()[-1]).view(-1, 1)

        w3 = self.w3(words)
        w3 = self.w33(w3)
        w3 = self.w333(w3)
        #w1 = self.bc1(w1)
        w_3 = F.avg_pool1d(w3, kernel_size=w3.size()[-1]).view(-1,)
        #print (w1.shape)
        u = w_3.repeat(w3.shape[-1]).view(-1, w3.shape[-1])
        wu_vec = (u-w3)**2
        wu = F.avg_pool1d(wu_vec, kernel_size=wu_vec.size()[-1]).view(-1,)
        w3 = torch.cat((w_3, wu), 0).view(-1, 1)

        #w3 = F.avg_pool1d(w3, kernel_size=w3.size()[-1]).view(-1, 1)

        w4 = self.w4(words)
        w4 = self.w44(w4)
        w4 = self.w444(w4)
        #w1 = self.bc1(w1)
        w_4 = F.avg_pool1d(w4, kernel_size=w4.size()[-1]).view(-1,)
        #print (w1.shape)
        u = w_4.repeat(w4.shape[-1]).view(-1, w4.shape[-1])
        wu_vec = (u-w4)**2
        wu = F.avg_pool1d(wu_vec, kernel_size=wu_vec.size()[-1]).view(-1,)
        w4 = torch.cat((w_4, wu), 0).view(-1, 1)

        #w4 = F.avg_pool1d(w4, kernel_size=w4.size()[-1]).view(-1, 1)

        #print (w1.shape, w2.shape)
        word_vectors = torch.cat((w1, w2, w3, w4), 0).view(-1,)
        #print (word_vectors.shape)
        #words_vectors = self.embed_word(word_vectors)
        protein_vector = self.embed_word(word_vectors)#word_vectors
        protein_vector = protein_vector.view(1, dim)
        #print (protein_vector.shape)
        #protein_vector = self.attention_cnn(compound_vector,
        #                                    word_vectors, layer_cnn)

        """Concatenate the above two vectors and output the interaction."""
        cat_vector = torch.cat((compound_vector, protein_vector), 1)
        for j in range(layer_output):
            cat_vector = torch.relu(self.W_out[j](cat_vector))
        interaction = self.W_interaction(cat_vector)

        return interaction

    def __call__(self, data, train=True):

        inputs, correct_interaction = data[:-1], data[-1]
        predicted_interaction = self.forward(inputs)

        if train:
            loss = F.cross_entropy(predicted_interaction, correct_interaction)
            return loss
        else:
            correct_labels = correct_interaction.to('cpu').data.numpy()
            ys = F.softmax(predicted_interaction, 1).to('cpu').data.numpy()
            predicted_labels = list(map(lambda x: np.argmax(x), ys))
            predicted_scores = list(map(lambda x: x[1], ys))
            return correct_labels, predicted_labels, predicted_scores

def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()


class Trainer(object):
    def __init__(self, model):
        self.model = model
        self.optimizer = optim.Adam(self.model.parameters(),
                                    lr=lr, weight_decay=weight_decay)

    def train(self, dataset):
        prog = printProgressBar(0, len(dataset), prefix = 'Progress:', suffix = 'Complete', length = 50)
        np.random.shuffle(dataset)
        N = len(dataset)
        loss_total = 0
        count = 1
        for data in dataset:
            printProgressBar(count, len(dataset), prefix = 'Progress:', suffix = 'Complete', length = 50)
            count += 1
            loss = self.model(data)
            self.optimizer.zero_grad()
            loss.backward()
            self.optimizer.step()
            loss_total += loss.to('cpu').data.numpy()
        return loss_total


class Tester(object):
    def __init__(self, model):
        self.model = model

    def test(self, dataset):
        global min_AUC
        N = len(dataset)
        T, Y, S = [], [], []
        for data in dataset:
            (correct_labels, predicted_labels,
             predicted_scores) = self.model(data, train=False)
            T.append(correct_labels)
            Y.append(predicted_labels)
            S.append(predicted_scores)
        AUC = roc_auc_score(T, S)

        y_true, y_score = T, S
        fpr, tpr, thresholds = roc_curve(y_true, y_score)
        if AUC > min_AUC:
                print ('Saving best model !')
                min_AUC = AUC
                np.save('roc_data_final.npy', {'fpr': fpr, 'tpr' : tpr, 'thres': thresholds, 'y_true' : T, 'y_score': S, 'pred_label' : Y})


        precision = precision_score(T, Y)
        recall = recall_score(T, Y)
        return AUC, precision, recall

    def save_AUCs(self, AUCs, filename):
        with open(filename, 'a') as f:
            f.write('\t'.join(map(str, AUCs)) + '\n')

    def save_model(self, model, filename):
        torch.save(model.state_dict(), filename)


def load_tensor(file_name, dtype):
    return [dtype(d).to(device) for d in np.load(file_name + '.npy')]


def load_pickle(file_name):
    with open(file_name, 'rb') as f:
        return pickle.load(f)


def shuffle_dataset(dataset, seed):
    np.random.seed(seed)
    np.random.shuffle(dataset)
    return dataset


def split_dataset(dataset, ratio):
    n = int(ratio * len(dataset))
    dataset_1, dataset_2 = dataset[:n], dataset[n:]
    return dataset_1, dataset_2

def expand_train(dataset_train):
    a,b,c,d = list(map(list, list(zip(*dataset_train))))
    e,f,g,h = [], [], [], []
    for i in range (len(a)):

        e.append(a[i])
        f.append(b[i])
        h.append(d[i])
        k = c[i].cpu().reshape((2, 6, 1500))
        ik = random.shuffle(list(range(6)))
        #k = k[ik].reshape((9000))
        i, j = k[0], k[1]
        i = i[ik]
        j = j[ik]

        cp = np.concatenate((i,j)).reshape((2, 9000))

        #print (c.shape)
        g.append(torch.from_numpy(cp).float().to(device))
        #if i%10000 ==0 and i !=0:
        #    break
            
    a += e#np.concatenate((a, np.array(e)))
    b += f#np.concatenate((b, np.array(f)))
    c += g#np.concatenate((c, np.array(g)))
    d += h#np.concatenate((d, np.array(h)))
    
    return list(zip(a,b,c,d))

if __name__ == "__main__":
    min_AUC = -99999
    """Hyperparameters."""
    (DATASET, radius, ngram, dim, layer_gnn, window, layer_cnn, layer_output,
     lr, lr_decay, decay_interval, weight_decay, iteration,
     setting) = sys.argv[1:]
    (dim, layer_gnn, window, layer_cnn, layer_output, decay_interval,
     iteration) = map(int, [dim, layer_gnn, window, layer_cnn, layer_output,
                            decay_interval, iteration])
    lr, lr_decay, weight_decay = map(float, [lr, lr_decay, weight_decay])

    """CPU or GPU."""
    if torch.cuda.is_available():
        device = torch.device('cuda')
        print('The code uses GPU...')
    else:
        device = torch.device('cpu')
        print('The code uses CPU!!!')

    """Load preprocessed data."""
    dir_input = ('../dataset/' + DATASET + '/input/'
                 'radius' + radius + '_ngram' + ngram + '/')
    compounds = load_tensor(dir_input + 'compounds', torch.LongTensor)
    adjacencies = load_tensor(dir_input + 'adjacencies', torch.FloatTensor)
    proteins = load_tensor(dir_input + 'proteins', torch.FloatTensor)

    interactions = load_tensor(dir_input + 'interactions', torch.LongTensor)
    fingerprint_dict = load_pickle(dir_input + 'fingerprint_dict.pickle')
    
    n_fingerprint = len(fingerprint_dict)


    """Create a dataset and split it into train/dev/test."""
    dataset = list(zip(compounds, adjacencies, proteins, interactions))
    dataset = shuffle_dataset(dataset, 1234)
    dataset_train, dataset_ = split_dataset(dataset, 0.8)
    dataset_dev, dataset_test = split_dataset(dataset_, 0.5)

    dataset_train = expand_train(dataset_train)

    print (len(dataset_train), len(dataset_test))

    """Set a model."""
    torch.manual_seed(1234)
    model = CompoundProteinInteractionPrediction().to(device)
    trainer = Trainer(model)
    tester = Tester(model)

    """Output files."""
    file_AUCs = '../output/result/AUCs--' + setting + '.txt'
    file_model = '../output/model/' + setting
    print (setting)
    AUCs = ('Epoch\tTime(sec)\tLoss_train\tAUC_dev\t'
            'AUC_test\tPrecision_test\tRecall_test')
    with open(file_AUCs, 'w') as f:
        f.write(AUCs + '\n')

    """Start training."""
    print('Training...')
    print(AUCs)
    start = timeit.default_timer()
    best = -9999

    for epoch in range(1, iteration):

        if epoch % decay_interval == 0:
            trainer.optimizer.param_groups[0]['lr'] *= lr_decay

        loss_train = trainer.train(dataset_train)
        AUC_dev = tester.test(dataset_dev)[0]
        AUC_test, precision_test, recall_test = tester.test(dataset_test)

        end = timeit.default_timer()
        time = end - start

        AUCs = [epoch, time, loss_train, AUC_dev,
                AUC_test, precision_test, recall_test]
        tester.save_AUCs(AUCs, file_AUCs)
        if AUC_test > best:
                best = AUC_test
                tester.save_model(model, file_model)
        AUCs = [round(i, 3) for i in AUCs]
        print('\t'.join(map(str, AUCs)))
