from src import htoq_binary
from itertools import combinations
from enum import Enum
import random
"""
double underscore in front of variable/function: private variable/function, not accessible outside the class
single underscore in front of variable/function: semi-private variable/function, should not access from outside the class, but accessible if you really want to
"""
class Helper():
    @staticmethod
    def get_coeff_list_by_degree(ho_ising):
        to_ret = {}
        if(len(ho_ising) > 2):
            to_ret[0] = [ho_ising[2]]
        to_ret[1] = ho_ising[0].values()
        for kv in ho_ising[1].items():
            (k,v) = kv
            if len(k) in to_ret:
                to_ret[len(k)].append(v)
            else:
                to_ret[len(k)] = [v]
        return to_ret
    
    @staticmethod
    def get_term_counts(ho_ising):
        termcount_dict = {}
        termcount_dict[1] = len(ho_ising[0])
        for key in ho_ising[1].keys():
            if len(key) in termcount_dict:
                termcount_dict[len(key)] +=1
            else:
                termcount_dict[len(key)] = 1
        return termcount_dict

    @staticmethod
    def get_variable_counts(ho_ising):
        variable_counts = {}
        neg_variables = {}
        for key in ho_ising[0].keys():
            variable_counts[key] = 1
            if(ho_ising[0][key] < 0):
                neg_variables[key] = 1
        
        for key in ho_ising[1].keys():
            for variable in key:
                if variable in variable_counts.keys():
                    variable_counts[variable] += 1
                else:
                    variable_counts[variable] = 1
            if(ho_ising[1][key] < 0):
                for variable in key:
                    if variable in neg_variables.keys():
                        neg_variables[variable] += 1
                    else:
                        neg_variables[variable] = 1

        return (variable_counts, neg_variables)


    @staticmethod
    ##for HBoolean we need another function RandomBooleanSol
    def createRandomIsingSol(variable_list):
        to_ret = {}
        for variable in variable_list:
            to_ret[variable] = 2*random.randint(0,1) -1
        return to_ret

    @staticmethod
    def createRandomBooleanSol(variable_list):
        to_ret = {}
        for variable in variable_list:
            to_ret[variable] = bool(random.getrandbits(1))
        return to_ret

    @staticmethod
    ##for HBoolean we need another function computeBoolean
    def computeIsing(ho_ising, solution):
        answer = ho_ising[2]
        for key in ho_ising[0].keys():
            if(solution[key] == 1):
                answer += ho_ising[0][key]
            else:
                answer -= ho_ising[0][key]
        
        for key in ho_ising[1].keys():
            literalvalue = 1
            for k in key:
                if(solution[k] == -1):
                    literalvalue = -literalvalue
            answer += ho_ising[1][key] * literalvalue     
        return answer

    @staticmethod
    def computeBoolean(hobo, solution):
        answer = hobo[2]
        for key in hobo[0].keys():
            if(solution[key]):
                answer += hobo[0][key]
        
        for key in hobo[1].keys():
            literalvalue = True
            for k in key:
                if(not solution[k]):
                    literalvalue = False
                    break
            if literalvalue:
                answer += hobo[1][key]
        return answer

class SolverState(Enum):
    NOT_INITIALIZED = 1
    INITIALIZED = 2
    PREPROCESSED = 3
    QUBOCONVERTED = 4
    QUBOSOLVED = 5
    POSTPROCESSED = 6

class HIsing(object):
    def __init__(self):
        self.should_scale = True
        self.errorPercent = 0.01
        self.isHeurSumDegree = False
        ##Do not need isIsingConversionFirst for HBoolean
        self.isIsingConversionFirst = True

        self._initial_sol_dict_converted_qubo = {}
        self._sol_qubo_dict = {}
        self._state = SolverState.NOT_INITIALIZED
        self._converted_problem = None
        self._bestKnownMinima = 0
        
        ##Change variable name of __ho_ising_orig
        self.__ho_ising_orig = None
        self._sol_dict = None
        self.__constraint_dict = None


    def initialize(self, ho_ising):
        if(len(ho_ising)!=3):
            raise ValueError("Ising file must be a three element array, First element is a dictionary " 
            +"containing linear terms, second element is a dictionary containing non linear terms, "
            +"third element is the constant")
        
        if self._state != SolverState.NOT_INITIALIZED:
            raise ValueError("HIsing has already been initialized.")
        self.__ho_ising_orig = ho_ising
        self._converted_problem = [ho_ising[0].copy(), ho_ising[1].copy(), ho_ising[2], 1.0]
        #self._sol_dict = self.__createAlloneSol()
        #The line below also updates self._sol_dict accordingly
        upperboundEstimate = self.__estimateUpperboundMinima()
        if upperboundEstimate < self._bestKnownMinima :
            self._bestKnownMinima = upperboundEstimate
        self._state = SolverState.INITIALIZED

    ##This function should probably be createAllTrueSol for HBoolean    
    def __createAlloneSol(self):
        ho_ising = self._converted_problem
        sol_dict ={}
        for key in ho_ising[0].keys():
            if(ho_ising[0][key]>0):
                sol_dict[key] = -1
            else:
                sol_dict[key] = 1
        for key in ho_ising[1].keys():
            for k in key:
                if k not in sol_dict:
                    sol_dict[k] = -1
        return sol_dict

    #The function below also updates self._sol_dict
    ##This function will change for HBoolean
    #This function generates the intial solution
    def __estimateUpperboundMinima(self):
        self._sol_dict = self.__createAlloneSol()
        to_ret = Helper.computeIsing(self.__ho_ising_orig, self._sol_dict)
        for _ in range(100):
            r_sol_dict = Helper.createRandomIsingSol(self._sol_dict.keys())
            new_estimate = Helper.computeIsing(self.__ho_ising_orig, r_sol_dict)
            if(new_estimate < to_ret):
                to_ret = new_estimate
                self._sol_dict = r_sol_dict
        return to_ret

    ##This function might change for HBoolean
    def __estimateLowerboundMinima(self):
        return - sum(map(abs, self._converted_problem[0].values())) - sum(map(abs, self._converted_problem[1].values())) + self._converted_problem[2]

    def __trimhobo_errorlevel(self, errorlevel):
        merged = list(self._converted_problem[0].items()) + list(self._converted_problem[1].items())
        error = 0
        for key, value in sorted(merged, key = lambda kv: abs(kv[1])):
            error += abs(value)
            if(error < errorlevel):
                self._converted_problem[0].pop(key, None)
                self._converted_problem[1].pop(key, None)
            else:
                break
        return


    def __trimhobo(self):
        upperbound = min(self._bestKnownMinima, 0)
        errorlevel = abs(upperbound) * self.errorPercent
        self.__trimhobo_errorlevel(errorlevel)

    #Assume ho_ising is a 3 element tuple, dictionary of linear terms, dictionary of nonlinear terms, constant
    ##This function should change for HBoolean
    def __update_hoising(self, variable_val, literal_term_dict):
        variable = variable_val[0]
        val = variable_val[1]
        if variable in literal_term_dict:
            for term in literal_term_dict[variable]:
                term_coeff = self._converted_problem[1].pop(term)
                term_as_list = list(term)
                term_as_list.remove(variable)
                new_term = tuple(term_as_list)
                if(len(new_term) > 1):
                    if(new_term in self._converted_problem[1]):
                        self._converted_problem[1][new_term] += term_coeff * val
                    else:
                        self._converted_problem[1][new_term] = term_coeff * val
                else:
                    if(new_term[0] in self._converted_problem[0]):
                        self._converted_problem[0][new_term[0]] += term_coeff * val
                    else:
                        self._converted_problem[0][new_term[0]] = term_coeff * val
        linear_coef = self._converted_problem[0].pop(variable) #As we are fixing a diagnoal term, it must be present as a linear term
        self._converted_problem[2] += linear_coef * val
        return 


    def __build_literal_term_dict(self):
        literal_term_dict = {}
        for key in self._converted_problem[1].keys():
            for k in key:
                if k in literal_term_dict:
                    literal_term_dict[k].append(key)
                else:
                    literal_term_dict[k] = [key]
        return literal_term_dict

    ##This function should change for HBoolean
    def __fix_Diagonal(self):
        literal_term_dict = self.__build_literal_term_dict()
        variable_list = tuple(self._converted_problem[0].keys())
        able_to_fix = False
        for k in variable_list:
            off_diag_sum = 0
            if(k in literal_term_dict):
                for term in literal_term_dict[k]:
                    off_diag_sum += abs(self._converted_problem[1][term])
            if(abs(self._converted_problem[0][k]) < off_diag_sum):
                continue
            able_to_fix = True
            if(self._converted_problem[0][k] > 0):
                self._sol_dict[k] = -1
            else:
                self._sol_dict[k] = 1
            self.__update_hoising((k, self._sol_dict[k]), literal_term_dict)
            literal_term_dict = self.__build_literal_term_dict()
        if(able_to_fix):
            self.__fix_Diagonal()
        return

    ##Change Constant namespace in HBoolean
    def preprocess(self):
        if(self._state != SolverState.INITIALIZED):
            raise ValueError("Calling preprocess, but state is not INITIALIZED")
        self.__trimhobo()
        self.__fix_Diagonal()
        self._state = SolverState.PREPROCESSED
        return
    
    ##don't need this for HBoolean
    def __2ising_to_qubo(self):
        for key in self._converted_problem[0].keys():
            coeff = self._converted_problem[0][key]
            self._converted_problem[0][key] = 2 * coeff
            self._converted_problem[2] += -coeff

        for key in self._converted_problem[1].keys():
            coeff = self._converted_problem[1][key]
            self._converted_problem[1][key] = 4*coeff
            htoq_binary.add_linterm_hoising(self._converted_problem, key[0], -2*coeff)
            htoq_binary.add_linterm_hoising(self._converted_problem, key[1], -2*coeff)
            self._converted_problem[2] += coeff
        
        for key in self._sol_dict.keys():
            if(self._sol_dict[key] == 1):
                self._sol_dict[key] = True
            else:
                self._sol_dict[key] = False
        return
    
    ##don't need this for HBoolean
    def __ho_ising_to_hobo(self):
        for key in self._converted_problem[0]:
            coeff = self._converted_problem[0][key]
            self._converted_problem[0][key] = 2 * coeff
            self._converted_problem[2] += -coeff
        temp_ho_dict = {}
        for key in self._converted_problem[1]:
            coeff = self._converted_problem[1][key]
            sign = 1
            for r in reversed(range(2,len(key)+1)):
                for term in combinations(key, r):
                    if(term in temp_ho_dict):
                        temp_ho_dict[term] += coeff * (2**r) * sign
                    else:
                        temp_ho_dict[term] = coeff * (2**r) * sign
                sign = -sign
            for k in key:
                if(k in self._converted_problem[0]):
                    self._converted_problem[0][k] += coeff * 2 * sign
                else:
                    self._converted_problem[0][k] = coeff * 2 * sign
            sign = -sign
            self._converted_problem[2] += sign * coeff
        self._converted_problem[1] = temp_ho_dict

        for key in self._sol_dict:
            if(self._sol_dict[key] == 1):
                self._sol_dict[key] = True
            else:
                self._sol_dict[key] = False
        return

    ##Check: probably remains same for HBoolean
    def __scale_qubo(self):
        count_zero = 0
        max_16_bit = 32678.0
        max_25_bit = 16777216.0
        scale_diag = 0
        scale_offdiag = 0
        
        if(len(self._converted_problem[0]) > 0):
            max_diag = float(max(map(abs, self._converted_problem[0].values())))
            scale_diag = max_25_bit/ max_diag

        if(len(self._converted_problem[1]) > 0):
            max_off_diag = float(max(map(abs, self._converted_problem[1].values())))
            scale_offdiag = max_16_bit/max_off_diag

        scale = min(scale_diag, scale_offdiag)
        if scale_diag == 0:
            scale = scale_offdiag
        if scale_offdiag == 0:
            scale = scale_diag
        
        diag_terms = list(self._converted_problem[0].keys())
        for key in diag_terms:
            self._converted_problem[0][key] = int(self._converted_problem[0][key]*scale)
            if(self._converted_problem[0][key] == 0):
                count_zero+=1
                self._converted_problem[0].pop(key)

        off_diag_terms = list(self._converted_problem[1].keys())
        for key in off_diag_terms:
            self._converted_problem[1][key] = int(self._converted_problem[1][key]*scale)
            if(self._converted_problem[1][key] == 0):
                count_zero+=1
                self._converted_problem[1].pop(key)
        self._converted_problem[2] = int(self._converted_problem[2]*scale)
        self._converted_problem[3] = scale
        #print "scale is: " + str(scale)
        return count_zero

    def __generate_initial_sol_vector(self):
        variable_count_dict = Helper.get_variable_counts(self._converted_problem)[0]
        for variable in variable_count_dict.keys():
            self._initial_sol_dict_converted_qubo[variable] = self._sol_dict[variable]
        
    ##don't need this for HBoolean
    def hoising_to_qubo(self):
        if(self._state != SolverState.PREPROCESSED):
            raise ValueError("Calling hoising_to_qubo, but state is not PREPROCESSED")
        multiplier = abs(self.__estimateLowerboundMinima())
        
        if self.isIsingConversionFirst:
            self.__constraint_dict = htoq_binary.convert_to_2ising(self._converted_problem, 
                multiplier, self._sol_dict, self.isHeurSumDegree)
            self.__2ising_to_qubo()
        else:
            self.__ho_ising_to_hobo()
            self.__constraint_dict = htoq_binary.convert_to2boolean(self._converted_problem,
                multiplier, self._sol_dict, self.isHeurSumDegree)

        if(self.should_scale):
            self.__scale_qubo()

        self.__generate_initial_sol_vector()
        self._state = SolverState.QUBOCONVERTED
        return

    def get_variable_count(self):
        return len(Helper.get_variable_counts(self._converted_problem)[0])

    ##don't need this for HBoolean
    def __check_constraints(self):
        for key in self.__constraint_dict:
            x1 = self._sol_dict[key[0]]
            x2 = self._sol_dict[key[1]]
            y = self._sol_dict[self.__constraint_dict[key][0]]
            z = self._sol_dict[self.__constraint_dict[key][1]]
            if(x1 * x2 != y ):
                raise ValueError("Constraint mismatch")
            if(x1 + x2 == -2):
                if(z != -1):
                    raise ValueError("constrain encoding is positive because of wrong free variable")
            else:
                if(z != 1):
                    raise ValueError("constrain encoding is positive because of wrong free variable")
        return

    def __check_constraints_boolean(self):
        for key in self.__constraint_dict:
            x1 = (self._sol_dict[key[0]] == 1)
            x2 = (self._sol_dict[key[1]] == 1)
            y = (self._sol_dict[self.__constraint_dict[key]] == 1)
            if((x1 and x2) != y ):
                raise ValueError("Constraint mismatch")      
        return

    ##This will change for HBoolean, don't need to convert from True to +1, False to -1
    def __merge_convert_solution2ising(self):
        for k, v in self._sol_qubo_dict.items():
            self._sol_dict[k] = v
        
        for k, v in self._sol_dict.items():
            if v:
                self._sol_dict[k] = 1
            else:
                self._sol_dict[k] = -1
            
    ##This will change for HBoolean        
    def postprocess(self):
        if(self._state != SolverState.QUBOSOLVED):
            raise ValueError("Calling postprocess, but state is not QUBOSOLVED")
        self.__merge_convert_solution2ising()
        if(self.isIsingConversionFirst):
            self.__check_constraints()
        else:
            self.__check_constraints_boolean()
        self._state = SolverState.POSTPROCESSED
    


    ##Name of this will change for HBoolean
    def get_ising_problem(self):
        return self.__ho_ising_orig
    
    def get_solution(self):
        if(self._state != SolverState.POSTPROCESSED):
            raise ValueError("Calling get_solution, but state is not POSTPROCESSED")
        return self._sol_dict
