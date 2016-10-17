def sortMolKeys(my_molecules):
    if ' and ' not in list(my_molecules.keys())[0]:
        my_sort = ['' for i in range(len(my_molecules))]
        nums = [int(i.split()[0]) for i in my_molecules]
        nums_sorted = ['%i'%i for i in sorted(nums)]
        for key in my_molecules:
            index = nums_sorted.index(key.split()[0])
            my_sort[index] = key
    else:
        my_sort = []
        nums = {}
        for i in my_molecules.keys():
            mol1 = int(i.split()[0])
            mol2 = int(i.split()[-2])
            if mol1 not in nums.keys(): nums[mol1] = []
            nums[mol1].append(mol2)
        for molStart in sorted(nums.keys()):
            for molFinish in sorted(nums[molStart]):
                for key in my_molecules:
                    if (key.split()[0] == '%i'%molStart) and (key.split()[-2] == '%i'%molFinish):
                        my_sort.append( key )
    return my_sort