from MCFlow.rdf import RDF, G, writeXY

class LocRDF(RDF):

    def __init__(self, file_name, *args):
        RDF.__init__(self, file_name, *args)

    def calculateRDF(self, m1, b1, c2_xyz, sbox):
        self.box = 'box%s'%sbox
        mol1 = 'mol%s'%m1
        self.ignored_frames = 0
        for iframe, FRAME_DATA in enumerate(self.frame_data):
            try:
                if (iframe+1)%(self.nframes//4) == 0:
                    print('%5.1f %% of RDF frames completed'%(100*(iframe+1)/self.nframes))
            except ZeroDivisionError:
                print('%5.1f %% of RDF frames completed'%(100*(iframe+1)/self.nframes))
            c1_xyz = []
            for i, each_mol in enumerate(FRAME_DATA[self.box][mol1]):
                for bead, coords in each_mol.items():
                    if bead == b1:
                        for each_bead in coords:
                            c1_xyz.append( list(map(float,each_bead) ) )
            c1_xyz = np.array(c1_xyz)
            if len(c1_xyz) == 0:
                self.ignored_frames += 1
                continue
            frame_boxlengths = self.boxlengths[iframe][self.box]
            self.boxLengths.append(frame_boxlengths)
            self.g_frame( c2_xyz, c1_xyz, frame_boxlengths)
        if self.ignored_frames > 0:
            print('Out of %i total frames, %i were ignored b/c of lack of mlcl'%(
                    self.nframes,self.ignored_frames))
        self.averageG()

class LocationG(G):
    def __init__(self):
        G.__init__(self)

    def checks(self):
        assert (self.args['bead1'] and self.args['mol1']
                 and self.args['bins'] and self.args['rmax']), 'Necessary input missing'
        assert self.args['boxes'], 'box needed for rdf'
        assert self.args['box'] == None, 'Use multiple boxes argument'
        assert self.args['mol2'] == None, 'Only 1 mol used for this file'
        assert self.args['bead2'] == None, 'Only 1 mol used in this files analysis'
        self.analysis_class = LocRDF

    def myCalcs(self, D, m1, b1):
        for box in self.args['boxes']:
            D.calculateRDF(m1, b1, self.getPositions(), box)
            file_path = '%s/%s/'%(self.args['path'],self.feed)
            file_info = 'box%s_mol%s-%s.dat'%(box, m1, b1)
            if len(self.args['dimensions']) == 3:
                writeXY(D.radii, D.g_average, file_path + 'rdf-' + file_info)
                writeXY(D.radii, D.n_average, file_path + 'nint-' + file_info)
            else:
                writeXY(D.radii, D.g_average, file_path +
                        'rdf%s-'%(''.join(self.args['dimensions'])) + file_info)
                writeXY(D.radii, D.n_average, file_path +
                        'nint%s-'%(''.join(self.args['dimensions'])) + file_info)

    def getPositions(self):
        raise NotImplemented

    def main(self):
        for feed in self.args['feeds']:
            self.feed = feed
            if self.args['verbosity'] > 0: print('-'*12 + 'Dir is %s'%self.feed + '-'*12)
            analysis = self.read_movies(self.args['rmax'], self.args['bins'], self.args['dimensions'])
            for m1, b1, in zip(self.args['mol1'], self.args['bead1']):
                self.myCalcs(analysis, m1, b1)

class MFI_Intersections(LocationG):
    def __init__(self):
        LocationG.__init__(self)

    def getPositions(self):
        a, b, c = [20.022, 19.899, 13.383]
        center = [10.011,4.9748,0.]
        x, y, z = center
        uc_sphere = [
                    [x,y,z],
                    [a/2-x,b/2+y,c/2+z],
                    [x,b/2-y,z],
                    [a/2+x,y,c/2-z],
                    [a-x,b-y,c-z],
                    [a/2+x,b/2-y,c/2-z],
                    [a-x,b/2+y,c-z],
                    [a/2-x,b-y,c/2+z]
                ]
        # replicate positions:
        sphere = []
        for i in range(2):
            for j in range(2):
                for k in range(3):
                    for center in uc_sphere:
                        x,y,z = center
                        sphere.append(
                            [
                                x+i*a, y+j*b, z+k*c
                            ]
                        )
        A, B, C = [a*2, b*2, c*3]
        indices_to_remove = []
        for i in range(len(sphere)):
            for j in range(i+1,len(sphere)):
                ii = sphere[i]
                jj = sphere[j]
                if ((abs(pbc(ii[0],A)-pbc(jj[0],A)) < 1e-5) and
                    (abs(pbc(ii[1],B)-pbc(jj[1],B)) < 1e-5) and
                    (abs(pbc(ii[2],C)-pbc(jj[2],C)) < 1e-5)):
                    # then they are the same:
                    indices_to_remove.append(j)
        return np.array([sphere[i] for i in range(len(sphere)) if i not in indices_to_remove])

import numpy as np
from MCFlow.parser import Results
from MCFlow.calc_tools import pbc

if __name__ == '__main__':
    M = MFI_Intersections()
    M.main()
