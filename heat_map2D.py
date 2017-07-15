def plot(x,y,dx,dy,z, labels, title):
    xlabel, ylabel, zlabel = labels
    plt.figure()
    plt.errorbar(x,y,xerr=dx, yerr=dy,fmt='None',ecolor='k',elinewidth=1.,
                 barsabove=False, zorder=1)
    plt.scatter(x,y,c=z,cmap='viridis',zorder=3,marker='s')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.loglog()
    plt.colorbar().set_label(zlabel)



import matplotlib.pyplot as plt


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Plot 2D heat map from columnar data in file')
    parser.add_argument('-f','--file',help='file: x, y, z, dx, dy, dz')
    args = vars(parser.parse_args())

    # get data
    x = []; y = []; z = []; dx = []; dy = []; dz = []
    xlabel, ylabel, zlabel = '','',''
    with open(args['file']) as f:
        for line in f:
            if line.startswith('#'):
                xlabel, ylabel, zlabel = line.split()[1:4]
            else:
                a,b,c,d,e,f = map(float,line.split()[:6])
                x.append(a)
                y.append(b)
                z.append(c)
                dx.append(d)
                dy.append(e)
                dz.append(f)
    plot(x,y,dx,dy,z,[xlabel,ylabel,zlabel],'z-mean')
    plot(x,y,dx,dy,dz,[xlabel,ylabel,zlabel],'z-error')
    plt.show()
