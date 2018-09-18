def initializeData(originalData, storageData):
    for key, value in originalData.items():
        if isinstance(value, dict):
            storageData[key] = {}
            initializeData(originalData[key], storageData[key])
        else:
            storageData[key] = []
    return storageData

class AnyProperty():
    def __init__(self, values):
        '''
        values must be nested dictionaries
        '''
        self.averaging = False
        self.data = initializeData(values, {})
        self.nblocks = 5
        self.block_mean = initializeData(values, {})
        self.block_std = initializeData(values, {})

    def addVals(self, newData):
        '''
        self.data is a nested dictionary ending in a list containing a value for each independent simulation.
        The mean of this simulation is the data point for the simulation.
        '''
        def appendValues(newRunData, storageData, blockMean, blockStdev):
            for key, value in newRunData.items():
                if isinstance(value, dict):
                    if key not in storageData.keys():
                        storageData[key] = {}
                        blockMean[key] = {}
                        blockStdev[key] = {}
                    appendValues(newRunData[key], storageData[key], blockMean[key], blockStdev[key])
                elif isinstance(value, list):
                    if key not in storageData.keys():
                        storageData[key] = []
                        blockMean[key] = []
                        blockStdev[key] = []
                    if len(value) > 0: # if non-empty
                        storageData[key].append(np.mean(value))
                        means, stdevs = [], []
                        block_length = math.floor(len(value)/self.nblocks)
                        for iblock in range(1,self.nblocks + 1):
                            start = int((iblock-1)*block_length)
                            stop = int(iblock*block_length)
                            block_data = value[start:stop]
                            means.append(np.mean(block_data))
                            stdevs.append(np.std(block_data))
                        blockMean[key].append(means)
                        blockStdev[key].append(stdevs)
                    else:
                        print('tried to take mean of empty slice; ignoring')
                elif isinstance(value, int) or isinstance(value, float):
                    if key not in storageData.keys(): storageData[key] = []
                    storageData[key].append(value)
                else:
                    print('Data type not accounted for in AnyProperty.addVals')
                    print(value)
        appendValues(newData, self.data, self.block_mean, self.block_std)

#    def getRunAvg(self, newData, newKey):
#        '''
#        self.data is a nested dictionary ending in a list containing a value for each independent simulation.
#        The mean of this simulation is the data point for the simulation.
#        '''
#        def averageValues(dataToAverage, allAverages):
#            for key, value in dataToAverage.items():
#                if isinstance(value, dict):
#                    allAverages[key] = {}
#                    averageValues(dataToAverage[key], allAverages[key])
#                elif isinstance(value, list):
#                    mean = np.mean(value)
#                    stdev = np.std(value)
#                    allAverages[key] = {'mean':mean, 'stdev':stdev}
#                    dataToAverage[key].clear()
#        if not self.averaging: # if havent started averaging need to make new dict
#            self.averages = {}
#            self.averaging = True
#        self.averages[newKey] = {}
#        averageValues(newData, self.averages[newKey])

    def avgVals(self, newKey):
        '''
        avgVals averages all indep sims for one simulation.
        self.averages contains the average values for each data point with stdevs.
        all of these averages then could make a plot
        '''
        def averageValues(dataToAverage, allAverages, blockMean, blockStd):
            for key, value in dataToAverage.items():
                if isinstance(value, dict):
                    allAverages[key] = {}
                    averageValues(dataToAverage[key], allAverages[key], blockMean[key], blockStd[key])
                elif isinstance(value, list):
                    mean = np.mean(value)
                    stdev = np.std(value)
                    allAverages[key] = {'mean':mean, 'stdev':stdev, 'raw':value[:]}
                    if key in blockMean.keys() and key in blockStd.keys():
                        allAverages[key]['block means'] = blockMean[key]
                        allAverages[key]['block std'] = blockStd[key]
                    dataToAverage[key].clear()
        if not self.averaging: # if havent started averaging need to make new dict
            self.averages = {}
            self.averaging = True
        self.averages[newKey] = {}
        averageValues(self.data, self.averages[newKey], self.block_mean, self.block_std)

    def __str__(self):
        def printValues(data, keysList = [], representation=''):
            nonlocal lineLength
            for key, value in sorted(data.items()):
                if key in self.averages.keys(): keysList = []
                if 'mean' in value.keys():
                    line = ''
                    if len(keysList) > lineLength:
                        keysList[-2] = keysList[-1]
                        keysList.pop(-1)
                    lineLength = len(keysList)
                    if value['stdev'] > 0.:
                        for k in keysList:
                            line += '%-10s  '%k
                        representation += '   %s    %-10s%-12e +/- %-12e\n'%(line, key, value['mean'], value['stdev'])
                else:
                    keysList.append(key)
                    representation = printValues(data[key], keysList=keysList, representation=representation)
            return representation
        lineLength = 8
        rep = printValues(self.averages)

        return rep


class Property(AnyProperty):
    '''
    Initialize and do operations on Properties that are stored as one total value 
    per box. For example: boxlength, x_BuOH, rho_BuOH_RE
    '''

    def __init__(self, nbox):
        '''
        Initialize each value in each box = 0
        '''
        self.data = {'box%i' % boxNum: [] for boxNum in range(1, nbox + 1)}

    def avgOverRuns(self, weights):  # should probably make unbound method so function
        self.avgOverRuns = {}
        for box in list(self.data.keys()):
            for (val, weight) in zip(self.data[box], weights):
                if box not in self.avgOverRuns.keys():
                    self.avgOverRuns[box] = np.multiply(val, weight)
                else:
                    self.avgOverRuns[box] += np.multiply(val, weight)
        return self.avgOverRuns

class MolProperty(Property, AnyProperty):
    '''
    Initialize and do operations on Properties that are stored as one value
    per molecule type per box. For example, N mlcls, number density
    '''

    def __init__(self, nbox, nmolty):
        '''
        Initialize each value for each mol in each box = 0
        '''
        self.data = {str(mol): {'box%i' % boxNum: [] for boxNum in range(1, nbox + 1)} for mol in range(1, nmolty + 1)}
        self.averaging = False

    def avgOverRuns(self, weights):  # should probably make unbound method so function
        assert abs(sum(weights)-1) <= 0.00001, 'Weights not calculated correctly, sum = {}'.format(sum(weights))
        self.avgOverRuns = {}
        for mlcl in list(self.data.keys()):
            self.avgOverRuns[mlcl] = {}
            for box in list(self.data[mlcl].keys()):
                self.avgOverRuns[mlcl][box] = np.sum([val * weight
                                                      for (val, weight) in zip(self.data[mlcl][box], weights)])
        return self.avgOverRuns
    def addVals(self, values):
        for outerKey in values.keys():
            for innerKey in values[outerKey].keys():
                if type(self.data[outerKey][innerKey]) == list:
                    self.data[outerKey][innerKey].append(np.mean(values[outerKey][innerKey]) )
                else:
                    self.data[outerKey][innerKey].append(values[outerKey][innerKey])

class AllProperty(AnyProperty):
    def __init__(self, outerKeys, innerKeys):
        self.data = {}
        self.averaging = False
        for outerKey in outerKeys:
            self.data[outerKey] = {}
            for innerKey in innerKeys:
                self.data[outerKey][innerKey] = []


class MolTransProperty(MolProperty):
    '''
    Initialize and do operations on Properties that are stored as one value
    per molecule type but only between boxes. For example, dG_trans
    '''

    def __init__(self):
        pass


import numpy as np
import math
