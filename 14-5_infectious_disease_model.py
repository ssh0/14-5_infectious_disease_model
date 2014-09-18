#! /usr/bin/env python
# -*- coding:utf-8 -*-
#
# written by Shotaro Fujimoto, June 2014.

from Tkinter import *
import numpy as np
import sys
import matplotlib.pyplot as plt


class Percolation:

    def __init__(self, L=61, p=0.5927):
        if L % 2 == 0:
            raise ValueError("lattice size L must be odd number.")
        self.sub = None
        self.L = L # lattice size
        self.p = p

    def perc_cluster(self, p):
        if p > 1 or p < 0:
            raise ValueError("site occupation probability must be 0 <= p <= 1")
        self.p = p
        self.lattice = np.zeros([self.L+2, self.L+2], dtype=int)
        self.lattice[:1,:] = self.lattice[:, :1] = -1
        self.lattice[self.L+1:,:] = self.lattice[:, self.L+1:] = -1
        center = (L/2) + 1
        self.lattice[center, center] = 1
        nextseed = [(center, center)]
        if self.sub is None or not self.sub.winfo_exists():
            lattice = self.lattice
            rn = np.random.random
            ne = [(0, -1), (0, 1), (-1, 0), (1, 0)]
            nnsite = set([(center+nx, center+ny) for nx, ny in ne])
            t = [0] # time
            S = [4] # a number of sites can be infected
            N = [1] # a number of infected sites
            percolate = False
            l = set([])
            while len(nnsite) != 0 and percolate == False:
                nextseed = []
                for nn in nnsite:
                    if rn() < p:
                        lattice[nn] = 1
                        nextseed.append(nn)
                    else:
                        lattice[nn] = -1
                nnsite = set([])
                for i, j in nextseed:
                    nnsite = nnsite | set([(i+nx, j+ny) for nx, ny in ne
                                            if lattice[i+nx, j+ny] == 0])
                    if i == 1:
                        l = l | set(['top'])
                    if i == self.L:
                        l = l | set(['bottom'])
                    if j == 1:
                        l = l | set(['left'])
                    if j == self.L:
                        l = l | set(['right'])
                
                if ('top' in l and 'bottom' in l) or \
                   ('right' in l and 'left' in l):
                    percolate = True
                
                t.append(t[-1]+1)
                S.append(len(nnsite))
                N.append(np.sum(lattice == 1))
            self.lattice = lattice[1:-1, 1:-1]
        
        return t, S, N

    def draw_canvas(self, rect, L):
        default_size = 640 # default size of canvas
        r = int(default_size/(2*L))
        fig_size = 2*r*L
        margin = 10
        sub = Toplevel()

        sub.title('figure  '+'(p=%s)' % str(self.p))
        self.canvas = Canvas(sub, width=fig_size+2*margin,
                             height=fig_size+2*margin)
        self.canvas.create_rectangle(margin, margin,
                                    fig_size+margin, fig_size+margin,
                                    outline='black', fill='white')
        self.canvas.pack()

        c = self.canvas.create_rectangle

        site = np.where(rect == 1)
        for m, n in zip(site[0], site[1]):
            c(2*m*r+margin, 2*n*r+margin,
              2*(m+1)*r+margin, 2*(n+1)*r+margin,
              outline='', fill='black')


class TopWindow:

    def quit(self):
        self.root.destroy()
        sys.exit()

    def show_window(self, pr, pushed, b4_pushed, auto):
        self.root = Tk()
        self.root.title('Percolation')
        f = Frame(self.root)
        self.label = Label(f, text='p =')
        self.label.pack(side='left')
        self.entry = Entry(f, width=20)
        self.entry.pack(side='left')
        self.entry.delete(0, END)
        self.entry.insert(0, 0.5927)
        self.entry.focus_set()
        
        b5 = Button(f, text='auto', command=auto)
        b5.pack(side='left', expand=YES, fill='x')

        b1 = Button(f, text='run', command=pushed)
        b1.pack(side='left', expand=YES, fill='x')

        b4 = Button(f, text='plot graph', command=b4_pushed)
        b4.pack(side='left', expand=YES, fill='x')

        b2 = Button(f, text='write canvas to sample.eps', command=pr)
        b2.pack(side='left', expand=YES, fill='x')

        b3 = Button(f, text='quit', command=self.quit)
        b3.pack(side='right', expand=YES, fill='x')

        f.pack(fill='x')

        self.root.mainloop()

def plot_graph(x_data, y_data, x_labels, y_labels,
               xscale, yscale, aspect):
    """ Plot the graph about y_data for each x_data.
    """
    d = len(y_data)
    if not len(x_data) == len(y_data) == len(x_labels) == len(y_labels)\
           == len(xscale) == len(yscale) == len(aspect):
        raise ValueError("Arguments must have the same dimension.")
    if d == 0:
        raise ValueError("At least one data for plot.")
    if d > 9:
        raise ValueError("""So much data for plot in one figure.
                            Please divide two or more data sets.""")

    fig = plt.figure(figsize=(9, 8))
    subplot_positioning = ['11', '21', '22', '22', '32', '32', '33', '33', '33']
    axes = []
    for n in range(d):
        lmn = int(subplot_positioning[d-1] + str(n+1))
        axes.append(fig.add_subplot(lmn))

    for i, ax in enumerate(axes):
        ymin, ymax = min(y_data[i]), max(y_data[i])
        ax.set_aspect(aspect[i])
        ax.set_xscale(xscale[i])
        ax.set_yscale(yscale[i])
        ax.set_xlabel(x_labels[i], fontsize=16)
        ax.set_ylabel(y_labels[i], fontsize=16)
        ax.set_ymargin(0.05)
        ax.plot(x_data[i], y_data[i], 'o-')

    fig.subplots_adjust(wspace=0.2, hspace=0.5)
    fig.tight_layout()
    plt.show()

if __name__ == '__main__':
    L = 61
    top = TopWindow()
    per = Percolation(L=L)
    count = 1

    def pr():
        global count
        p = float(top.entry.get())
        d = per.canvas.postscript(file="figure_%d(p=%s).eps" % (count, str(p)))
        print "saved the figure to a eps file"
        count += 1

    def pushed():
        global t, S, N
        p = float(top.entry.get())
        t, S, N = per.perc_cluster(p)
        per.draw_canvas(per.lattice, L)

    def b4_pushed():
        x_data = [t[1:-1]]*2
        y_data = [N[1:-1], S[1:-1]]
        x_labels = [r'$t$', r'$t$']
        y_labels = [r'$N$', r'$S$']
        plot_graph(x_data, y_data, x_labels, y_labels,
                    ['log', 'linear'], ['log', 'linear'], ['auto']*2)

    def auto():
        trial = 200
        p = np.linspace(0.3, 1.0, 50)
        N_p = []
        S_p = []
        perc_rate = []
        for _p in p:
            N_p_ = []
            S_p_ = []
            perc = 0
            for i in range(trial):
                t, S, N = per.perc_cluster(_p)
                if S[-1] != 0:
                    perc += 1
                N_p_.append(N[-1])
                S_p_.append(S[-1])
            N_p.append(np.average(np.array(N_p_)))
            S_p.append(np.average(np.array(S_p_)))
            perc_rate.append(float(perc)/trial)
        plot_graph([p], [perc_rate], [r'$p$'], [r'$P(p)$'],
                    ['linear'], ['linear'], ['auto'])
        
    top.show_window(pr, pushed, b4_pushed, auto)

