from calc_tools import calculate_torsion, get_vector, determine_gauche

class NoTorsions(BaseException):
    def __init__(self):
        BaseException.__init__(self)

def get_tors(data, abc):
    '''
    data: molecule data {'H':[[,,]...],'O:[[]' etc.
    '''
    tol = 1e-5
    torsion_molec_data = {}
    lengths = [1.43*1.43, 1.54*1.54]
    for b1, v1 in data.items():
        for xyz1 in v1:
            for b2, v2 in data.items():
                for xyz2 in v2:
                    my_vector1 = np.array(get_vector(xyz1, xyz2, abc))
                    v_2 = sum(r*r for r in my_vector1)
                    tests = [abs(v_2-i)<tol for i in lengths]
                    if True in tests:
                        for b3, v3 in data.items():
                            for xyz3 in v3:
                                my_vector2 = np.array(get_vector(xyz2, xyz3, abc))
                                if np.linalg.norm(my_vector1+my_vector2,2)<tol: continue
                                v_3 = sum(r*r for r in my_vector2)
                                tests = [abs(v_3-i)<tol for i in lengths]
                                if True in tests:
                                    for b4, v4 in data.items():
                                        for xyz4 in v4:
                                            my_vector3 = np.array(get_vector(xyz3, xyz4, abc))
                                            if np.linalg.norm(my_vector2+my_vector3,2)<tol: continue
                                            v_4 = sum(r*r for r in my_vector3)
                                            tests = [abs(v_4-i)<tol for i in lengths]
                                            if True in tests:
                                                my_tor = b1 + b2 + b3 + b4
                                                vectors = [my_vector1, my_vector2,
                                                                    my_vector3]
                                                torsion = calculate_torsion(vectors)
                                                tors_name = ''.join(sorted(my_tor))
                                                if tors_name not in torsion_molec_data.keys():
                                                    torsion_molec_data[tors_name] = []
                                                torsion_molec_data[tors_name].append(torsion)
    return torsion_molec_data

from MCFlow.structure_analysis import Struc, Movie
from MCFlow.channelLocation import FindRegion

class Torsion(Movie):
    def __init__(self, file_name, *args):
        Movie.__init__(self, file_name, *args)

    def getTorsions(self, uc_vectors):
        '''
        :param mlcl: String molecule number.
        :param box: String box number.
        of all molecules in this box
        '''
        regions = FindRegion(uc_vectors)
        assert self.nframes == len(self.frame_data), 'Error in adding frames'
        print('Total amount of frames analyzed was %i'%self.nframes)
        torsion_histogram = {}
        all_trans_histogram = {}
        no_torsion_mols = []
        for iframe, FRAME_DATA in enumerate(self.frame_data):
            try:
                if (iframe+1)%(self.nframes//4) == 0:
                    print('%5.1f %% of frames analyzed for HB'%(100*(iframe+1)/self.nframes))
            except ZeroDivisionError:
                print('%5.1f %% of frames analyzed for HB'%(100*(iframe+1)/self.nframes))
            seed = self.frame_seed[iframe]
            for box in FRAME_DATA.keys():
                for mlcl, data in FRAME_DATA[box].items():
                    if mlcl in no_torsion_mols: continue
                    try:
                        for each_molecule in data:
                            COM = each_molecule.pop('COM')
                            if box == 'box1':
                                abc = uc_vectors
                                COM = COM[0]
                                name = regions.getChannel(COM)
                            else:
                                name = box
                                abc = self.boxlengths[iframe][box]
                            if name not in torsion_histogram.keys():
                                torsion_histogram[name] = {}
                                all_trans_histogram[name] = {}
                            torsion_data = get_tors(each_molecule,abc)
                            if len(torsion_data.keys()) == 0:
                                raise NoTorsions
                            if mlcl not in torsion_histogram[name].keys():
                                torsion_histogram[name][mlcl] = [[] for
                                            i in range(len(set(self.frame_seed)))]
                                all_trans_histogram[name][mlcl] = [[] for
                                            i in range(len(set(self.frame_seed)))]
                            # combine torsion types
                            all_tors_data = []
                            for torsion_type, tdata in torsion_data.items():
                                all_tors_data += tdata
                            l_gauche = [determine_gauche(i) for i in all_tors_data]
                            torsion_histogram[name][mlcl][seed-1] += all_tors_data
                            if l_gauche.count(False) == len(l_gauche):
                                all_trans_histogram[name][mlcl][seed-1].append(True)
                            else:
                                all_trans_histogram[name][mlcl][seed-1].append(False)
                    except NoTorsions:
                        print('No torsions found for %s'%mlcl)
                        no_torsion_mols.append(mlcl)
                        continue
        return torsion_histogram, all_trans_histogram

    def getFracGauche(self, uc_vectors):
        torsion_histogram, all_trans_h = self.getTorsions(uc_vectors)
        fraction_gauche = {}
        for box, d1 in torsion_histogram.items():
            if box not in fraction_gauche.keys(): fraction_gauche[box]  = {}
            for mlcl, d2 in d1.items():
                if mlcl not in fraction_gauche[box].keys():
                    fraction_gauche[box][mlcl] = {'raw data':[]}
                hist = []
                for seed_data in d2:
                    hist += seed_data
                    l_gauche = [determine_gauche(i) for i in seed_data]
                    if len(l_gauche) > 0:
                        f_gauche = l_gauche.count(True)/len(l_gauche)
                        fraction_gauche[box][mlcl]['raw data'].append(f_gauche)
                f_all_trans = []
                for seed_data in all_trans_h[box][mlcl]:
                    if len(seed_data) > 0:
                        f_all_trans.append(seed_data.count(True)/len(seed_data))
                fraction_gauche[box][mlcl]['mean'] = np.mean(fraction_gauche[box][mlcl]['raw data'])
                fraction_gauche[box][mlcl]['all trans mean'] = np.mean(f_all_trans)
                fraction_gauche[box][mlcl]['all trans stdev'] = np.std(f_all_trans)
                fraction_gauche[box][mlcl]['stdev'] = np.std(fraction_gauche[box][mlcl]['raw data'])
#               fraction_gauche[box][mlcl]['hist'] = hist
        return fraction_gauche

class T(Struc):
    def __init__(self):
        my_parser = Results()
        my_parser.parser.add_argument('-ID','--name',help='Name of output db',
                               type=str)
        my_parser.parser.add_argument('-abc','--vectors', help='unit cell vectors.'+ 
                                        'For folding coordinates into a unit cell',
                                      type = float, nargs = '+',
                                        default = [20.022,19.899,13.383])
        my_args = vars(my_parser.parse_args())
        self.args = my_args
        self.checks()

    def checks(self):
        assert self.args['vectors'], 'unit cell vectors needed for folding'
        self.analysis_class = Torsion


    def myCalcs(self, D ):
        D.foldMovieToUC(self.args['vectors'])
        D.averages = {self.feed : D.getFracGauche(self.args['vectors'])}
        if self.args['name']:
            outputDB(self.args['path'],[self.feed],self.args['type'],{self.args['name']: D } )

    def main(self):
        for feed in self.args['feeds']:
            self.feed = feed
            if self.args['verbosity'] > 0: print('-'*12 + 'Dir is %s'%self.feed + '-'*12)
            analysis = self.read_movies()
            self.myCalcs(analysis)


# from MCFlow.runAnalyzer import what2Analyze
from MCFlow.file_formatting.reader import Movie
# from MCFlow.file_formatting.writer import xyz
from MCFlow.parser import Results
import numpy as np
from MCFlow.getData import outputDB
# from MCFlow.calc_tools import calculate_distance2

if __name__ == '__main__':
    M = T()
    M.main()
