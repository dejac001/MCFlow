from MCFlow.runAnalyzer import getRhoByMCC
from MCFlow.getData import outputDB
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from MCFlow.runAnalyzer import calc95conf, checkRun
import shelve

def gaussian(x,alpha,peakdens,variation):

    #Input: An array of x-values, and three parameters for a Gaussian fit
    #Output: An array of y-values consistent with the specified Gaussian distribution

    return alpha*np.exp( (-(x-peakdens)**2) / (2.0*(variation**2)))


class DPD:
    def __init__(self, **kwargs):

        self.indep = kwargs['indep']
        self.mol = kwargs['molNum']
        self.feeds = kwargs['feeds']
        self.fracPeakMax = 0.75 # only fit data 75 % of peak maximum
        self.conversionFactor = 1/1000
        self.numBins = 2000
        self.verbosity = kwargs['verbosity']
        self.guessStart = kwargs['guessStart']
        self.path = kwargs['path']
        self.type = kwargs['type']
        self.interval = kwargs['interval']


    def initializeData(self):
        if self.verbosity > 0: print('-' * 12 + 'Dir is %s' % self.feed + '-' * 12)
        self.rho_by_MCC_by_seed = []
        for seed in self.indep:
            self.rho_by_MCC_by_seed.append(
                    getRhoByMCC(self.feed, seed, self.path, self.type, self.guessStart,
                                self.interval, self.verbosity)
                )



    def plot(self, density_probability_distributions):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        densities_overall = {'low dens':[], 'high dens':[]}
        for i, dpd in enumerate(density_probability_distributions):
            for key in 'low dens','high dens':
                bins = dpd[key]['bins']
                vals = dpd[key]['vals']
                gauss_x = np.linspace(dpd[key]['min rho gauss'], dpd[key]['max rho gauss'])
                gauss_y = gaussian(gauss_x, *dpd[key]['parameters'])
                kwargs = {'color':'C%i'%i}
                densities_overall[key].append(dpd[key]['parameters'][1])
                ax.plot(bins, vals,ls='dotted', **kwargs)
                if key == 'high dens': kwargs['label'] = 'run %i'%self.indep[i]
                ax.plot(gauss_x, gauss_y, ls='solid', **kwargs)
        ax.set_ylabel('$P(\\rho)$')
        ax.set_xlabel('$\\rho\mathrm{[\;molec/\AA^3]}$')
        ax.legend()
        for key, val in densities_overall.items():
            with shelve.open('rho-data.db',writeback=True) as db:
                run = [i for i in list(db[self.feed].keys()) if self.type in i][0]
                db[self.feed][run][self.mol][key] = {'mean':np.mean(val),
                                            '95conf': calc95conf(np.std(val), len(self.indep))}
        fig.savefig('dpd-%s.pdf'%(self.feed.replace('/','-')))

    def main(self):
        for self.feed in self.feeds:
            all_data = []
            self.initializeData()
            for i in range(len(self.indep)):
                my_hist, my_bins = self.dpd(i)
                all_data.append(self.fit_gaussian_2peaks(my_hist, my_bins))
            self.plot(all_data)

    def dpd(self, indep_seed):
        # Input: A Python list of densities
        # Output: Two coordinated Python lists: a list of densities
        # and a corresponding list of normalized probabilities

        # make density list
        all_box_densities = []
        for box, val in self.rho_by_MCC_by_seed[indep_seed][self.mol].items():
            all_box_densities += [i* self.conversionFactor for i in val]

        densities = np.array(all_box_densities)
        hist, edges = np.histogram(densities, bins=self.numBins, density=True)
        bin_means = [(edges[i] + edges[i+1])/2 for i in range(len(edges)-1)]
        return hist, bin_means

    def fit_gaussian_2peaks(self, probs, dens_bins):

        # Input: A Python list of densities
        # Output: The parameters for two Gaussian fits of the density

        assert len(dens_bins) == len(probs), 'Inconsistend input densities'

        middle_density = (dens_bins[-1] - dens_bins[0]) / 2.0

        data = {'low dens':{'bins':[],'vals':[]}, 'high dens':{'bins':[],'vals':[]}}
        for i, rho in enumerate(dens_bins):
            if rho < middle_density:
                key = 'low dens'
            else:
                key = 'high dens'
            data[key]['bins'].append(rho)
            data[key]['vals'].append(probs[i])

        for key, values in data.items():
            indices,  = np.where(values['vals'] >= self.fracPeakMax*np.max(values['vals']))
            x = [data[key]['bins'][i] for i in indices]
            y = [data[key]['vals'][i] for i in indices]
            data[key]['min rho gauss'] = np.min(x)
            data[key]['max rho gauss'] = np.max(x)
            initial_guess = [np.max(y), np.sum(i*j for i,j in zip(x,y))/np.sum(y), np.std(y)]
            params, errors = curve_fit(gaussian, x, y, p0=initial_guess, maxfev=16000)
            data[key]['parameters'] = params
        return data

import numpy as np

if __name__ == '__main__':
    from parser import Main
    parser = Main()
    parser.other()
    parser.parser.add_argument('-m','--molNum',help='mol number',type=str)
    args = vars(parser.parse_args())
    I = DPD(**args)
    I.main()