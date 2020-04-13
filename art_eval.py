import sys
from  src import solveHIsing
import pickle
from src.solveHIsing import Helper
import copy
from itertools import combinations
import argparse


def add_term_ho(ho, term, coeff):
    if len(term) == 0:
        #print("this should not happen")
        ho[2] += coeff
    elif len(term) == 1:
        if(term[0] in ho[0]):
            ho[0][term[0]] += coeff
        else:
            ho[0][term[0]] = coeff
    else:
        if term in ho[1]:
            ho[1][term] += coeff
        else:
            ho[1][term] = coeff
    return

def removeZeros(ho):
    for x in list(ho[0].keys()):
        if ho[0][x] == 0.0:
            ho[0].pop(x)

    for x in list(ho[1].keys()):
        if ho[1][x] == 0.0:
            ho[1].pop(x)
    return

def hobo_to_ho_ising(hobo):
    ho_ising = [{},{},0]
    ho_ising[2] += hobo[2]
    #print("hobo2ising Initial Const: {}".format(ho_ising[2]))

    for kv in hobo[0].items():
        (linterm, coeff) = kv
        ho_ising[2] += coeff/2.0
        #print("hobo2ising Current Const (linterms): {}".format(ho_ising[2]))
        add_term_ho(ho_ising, (linterm,), coeff/2.0)

    for kv in hobo[1].items():
        (term, coeff) = kv
        coeff = float(coeff)/(2**len(term))
        for r in range(1,len(term)+1):
            for tt in combinations(term, r):
                add_term_ho(ho_ising, tt, coeff)
        ho_ising[2] += coeff
        #print("hobo2ising Current Const (higherterms): {}".format(ho_ising[2]))

    removeZeros(ho_ising)
    return ho_ising



def main(args):
    ho_ising = pickle.load(open(args.fpath, 'rb'))
    print("Term counts original")
    print(Helper.get_term_counts(ho_ising))
    print("Number of Variables original: {}".format(len(Helper.get_variable_counts(ho_ising)[0])))


    hs1 = solveHIsing.HIsing()
    hs1.should_scale = False
    hs1.errorPercent = 0
    hs1.initialize(ho_ising)
    hs1.preprocess()
    print("Term counts after preprocessing")
    print(Helper.get_term_counts(hs1._converted_problem))
    print("Number of Variables: {}".format(len(Helper.get_variable_counts(hs1._converted_problem)[0])))


    hs2 = copy.deepcopy(hs1)
    hs3 = copy.deepcopy(hs1)
    hs4 = copy.deepcopy(hs1)

    hs2.isHeurSumDegree = True
    hs3.isIsingConversionFirst = False
    hs4.isIsingConversionFirst = False
    hs4.isHeurSumDegree = True

    hs1.hoising_to_qubo()
    print("Term counts after Conversion. Ising conversion, Algorithm 1")
    print(Helper.get_term_counts(hs1._converted_problem))
    print("Number of Variables: {}".format(len(Helper.get_variable_counts(hs1._converted_problem)[0])))

    hs2.hoising_to_qubo()
    print("Term counts after Conversion. Ising conversion, Algorithm 2")
    print(Helper.get_term_counts(hs2._converted_problem))
    print("Number of Variables: {}".format(len(Helper.get_variable_counts(hs2._converted_problem)[0])))

    hs3.hoising_to_qubo()
    print("Term counts after Conversion. Boolean conversion, Algorithm 1")
    print(Helper.get_term_counts(hs3._converted_problem))
    print("Number of Variables: {}".format(len(Helper.get_variable_counts(hs3._converted_problem)[0])))

    hs4.hoising_to_qubo()
    print("Term counts after Conversion. Boolean conversion, Algorithm 2")
    print(Helper.get_term_counts(hs4._converted_problem))
    print("Number of Variables: {}".format(len(Helper.get_variable_counts(hs4._converted_problem)[0])))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "fpath", type=str, nargs='?', help="path to HOBO-Ising file")
    args = parser.parse_args()
    main(args)
