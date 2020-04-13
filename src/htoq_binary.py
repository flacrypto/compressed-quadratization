def _get_max_variable(ho_ising):
    max_key = 0
    if len(ho_ising[0]) > 0:
        max_key = max(ho_ising[0])
    if(len(ho_ising[1]) == 0):
        return max_key
    max_key2 = max(map(lambda x : max(x), ho_ising[1]))
    return max(max_key, max_key2)

def _xy_to_pairxy(x,y):
    if (x < y):
        return (x,y)
    else:
        print("We Should never be here, some thing wrong")
        return (y,x)

def _term_to_pair_list(term):
    ell = len(term)
    if(ell <= 2):
        return []
    list_to_ret = []
    for i in range(ell):
        for j in range(i+1, ell):
            list_to_ret.append(_xy_to_pairxy(term[i], term[j]))
    return list_to_ret

def _add_term_bipartite(term, bp_graph):
    pair_list = _term_to_pair_list(term)
    for xij in pair_list:
        if xij in bp_graph:
            bp_graph[xij].append(term)
        else:
            bp_graph[xij] = [term]
    return

def _delete_term_bipartite(term, bp_graph, pair_ij):
    pair_list = _term_to_pair_list(term)
    pair_list.remove(pair_ij)
    for xij in pair_list:
        bp_graph[xij].remove(term)
        if(len(bp_graph[xij]) == 0):
            bp_graph.pop(xij)
    return
    
def _build_ho_bipartite(ho_ising):
    to_ret = {}
    for term in ho_ising[1].keys():
        _add_term_bipartite(term, to_ret)
    return to_ret

def _get_sum_degrees(list_of_terms):
    ret = 0
    for term in list_of_terms:
        ret += len(term) - 1
    return ret

def _find_best_pair(bp_graph, isHeurSumDegree):
    if isHeurSumDegree:
        return max(bp_graph.items(), key=lambda kv: _get_sum_degrees(kv[1]))[0]
    else:
        return max(bp_graph.items(), key=lambda kv: len(kv[1]))[0]
    

def _add_linterm_hoising(ho_ising, term , coef):
    if term in ho_ising[0]:
        ho_ising[0][term] += coef
    else:
        ho_ising[0][term] = coef
    return

def _add_nonlinterm_hoising(ho_ising, term , coef):
    if term in ho_ising[1]:
        ho_ising[1][term] += coef
    else:
        ho_ising[1][term] = coef
    return

def _add_constraint(pair_ij, ho_ising, new_val, multiplier):
    ho_ising[2] += 4*multiplier
    _add_linterm_hoising(ho_ising, pair_ij[0], multiplier)
    _add_linterm_hoising(ho_ising, pair_ij[1], multiplier)
    ho_ising[0][new_val] = -multiplier
    ho_ising[0][new_val+1] = -2*multiplier
    _add_nonlinterm_hoising(ho_ising, pair_ij, multiplier)
    ho_ising[1][(new_val, new_val+1)] = 2*multiplier
    ho_ising[1][(pair_ij[1], new_val)] = -multiplier
    ho_ising[1][(pair_ij[0], new_val)] = -multiplier
    ho_ising[1][(pair_ij[0], new_val+1)] = -2*multiplier
    ho_ising[1][(pair_ij[1], new_val+1)] = -2*multiplier  
    return

def _add_constraint_boolean(pair_ij, ho_boolean, new_val, multiplier):
    _add_linterm_hoising(ho_boolean, new_val, 3*multiplier)
    _add_nonlinterm_hoising(ho_boolean, pair_ij, multiplier)
    ho_boolean[1][(pair_ij[0], new_val)] = -2*multiplier
    ho_boolean[1][(pair_ij[1], new_val)] = -2*multiplier
    return

def _replace_pair_ho_ising(pair_ij, bp_graph, ho_ising, new_va, multiplier):
    for term in bp_graph[pair_ij]:
        coeff = ho_ising[1].pop(term)
        termL = list(term)
        termL.remove(pair_ij[0])
        termL.remove(pair_ij[1])
        termL.append(new_va)
        ho_ising[1][tuple(termL)] = coeff
    _add_constraint(pair_ij, ho_ising, new_va, multiplier)
    return

def _replace_pair_ho_boolean(pair_ij, bp_graph, ho_boolean, new_va, multiplier):
    for term in bp_graph[pair_ij]:
        coeff = ho_boolean[1].pop(term)
        termL = list(term)
        termL.remove(pair_ij[0])
        termL.remove(pair_ij[1])
        termL.append(new_va)
        ho_boolean[1][tuple(termL)] = coeff
    _add_constraint_boolean(pair_ij, ho_boolean, new_va, multiplier)
    return

def _replace_pair_bpgraph(pair_ij, bp_graph, ho_ising, new_va):
    for term in bp_graph[pair_ij]:
        _delete_term_bipartite(term, bp_graph, pair_ij)
        termL = list(term)
        termL.remove(pair_ij[0])
        termL.remove(pair_ij[1])
        termL.append(new_va)
        _add_term_bipartite(tuple(termL), bp_graph)
    bp_graph.pop(pair_ij)
    return

def _replace_pair(pair_ij, bp_graph, ho_ising, new_var, multiplier):
    _replace_pair_ho_ising(pair_ij, bp_graph, ho_ising, new_var, multiplier)
    _replace_pair_bpgraph(pair_ij, bp_graph, ho_ising, new_var)
    return

def _replace_pair_boolean(pair_ij, bp_graph, ho_ising, new_var, multiplier):
    _replace_pair_ho_boolean(pair_ij, bp_graph, ho_ising, new_var, multiplier)
    _replace_pair_bpgraph(pair_ij, bp_graph, ho_ising, new_var)
    return

def _is_sorted(tt):
    for i in range(len(tt)-1):
        if(tt[i] >= tt[i+1]):
            return False
    return True

def _convert_to_sorted(ho_ising):
    for key in list(ho_ising[1].keys()):
        if(not _is_sorted(key)):
            coef = ho_ising[1].pop(key)
            sorted_tuple = tuple(sorted(key))
            if( sorted_tuple in ho_ising[1] ):
                raise ValueError("duplicate entry in ho_ising: " + str(key))
            ho_ising[1][sorted_tuple] = coef
    return

def _are_tuples_sorted(ho_ising):
    #print "Type of ho_ising: " + str(type(ho_ising))
    #print "Type of ho_ising[1]: " + str(type(ho_ising[1]))
    for key in ho_ising[1].keys():
        if(_is_sorted(key)==False):
            return False
    return True

def _update_initial_values(initial_val_dict, pair_ij, new_var):
    initial_val_dict[new_var] = initial_val_dict[pair_ij[0]] * initial_val_dict[pair_ij[1]]
    if(initial_val_dict[pair_ij[0]] + initial_val_dict[pair_ij[1]] == -2):
        initial_val_dict[new_var + 1] = -1
    else:
        initial_val_dict[new_var + 1] = 1
    return

def _update_initial_values_boolean(initial_val_dict, pair_ij, new_var):
    initial_val_dict[new_var] = initial_val_dict[pair_ij[0]] and initial_val_dict[pair_ij[1]]
    return

def convert_to_2ising(ho_ising, abs_lower_bound, initial_val_dict, isHeurSumDegree):
    if(not _are_tuples_sorted(ho_ising)):
        #print("Tuples not sorted. Now sorting")
        _convert_to_sorted(ho_ising)
        if(not _are_tuples_sorted(ho_ising)):
            raise ValueError("something wrong in sorting")
        
    bp_graph = _build_ho_bipartite(ho_ising)
    new_var = _get_max_variable(ho_ising) + 1
    constraint_dict = {}
    while(len(bp_graph) > 0):
        best_pair = _find_best_pair(bp_graph, isHeurSumDegree)
        _replace_pair(best_pair, bp_graph, ho_ising, new_var, abs_lower_bound)
        _update_initial_values(initial_val_dict, best_pair, new_var)
        constraint_dict[best_pair] = (new_var, new_var + 1)
        new_var += 2
    return constraint_dict

def convert_to2boolean(ho_boolean, abs_lower_bound, initial_val_dict, isHeurSumDegree):
    if(not _are_tuples_sorted(ho_boolean)):
        #print("Tuples not sorted. Now sorting")
        _convert_to_sorted(ho_boolean)
        if(not _are_tuples_sorted(ho_boolean)):
            raise ValueError("something wrong in sorting")
        
    bp_graph = _build_ho_bipartite(ho_boolean)
    new_var = _get_max_variable(ho_boolean) + 1
    constraint_dict = {}
    while(len(bp_graph) > 0):
        best_pair = _find_best_pair(bp_graph, isHeurSumDegree)
        _replace_pair_boolean(best_pair, bp_graph, ho_boolean, new_var, abs_lower_bound)
        _update_initial_values_boolean(initial_val_dict, best_pair, new_var)
        constraint_dict[best_pair] = new_var
        new_var += 1
    return constraint_dict    

def add_linterm_hoising(ho_ising, term , coef):
    _add_linterm_hoising(ho_ising, term , coef)