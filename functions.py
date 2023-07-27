import numpy as np
from sympy import simplify, Symbol
from random import randint
from random import random
from math import isclose
from sys import stdout
from time import time
from itertools import combinations
from quivers import generate_quiver
'''
gr_36_quiver=np.array(
  [[0,1,1,-1,-1,0,0,0,0,0],
   [-1,0,0,1,0,1,-1,0,0,0],
   [-1,0,0,1,0,0,0,1,-1,0],
   [1,-1,-1,0,0,0,1,0,1,-1],
   [1,0,0,0,0,0,0,0,0,0],
   [0,-1,0,0,0,0,0,0,0,0],
   [0,1,0,-1,0,0,0,0,0,0],
   [0,0,-1,0,0,0,0,0,0,0],
   [0,0,1,-1,0,0,0,0,0,0],
   [0,0,0,1,0,0,0,0,0,0] ]
)
gr_36_vertices=['a236','a235','a136','a356','a123','a234','a345','a126','a156','a456',]

all_plucker_coords_gr36=gr_36_vertices+list(set([f'a{c[0]}{c[1]}{c[2]}' for c in combinations(list(range(1,7)),3)])-set(gr_36_vertices))

all_plucker_coords_gr36_s=[Symbol(x) for x in all_plucker_coords_gr36]

all_plucker_coords_gr36_s+=[Symbol('ax1'),Symbol('ax2')]

gr_36_vertices_s=[Symbol(x) for x in gr_36_vertices]


dummy_matrix=np.array([[n+6*i+random() for n in range(6)] for i in range(3)])

#These are the associated numbers to each vertex as listed above
gr_36_vertices_num=[np.linalg.det(np.array([ [row[int(i)-1] for i in list(plucker[1:])] for row in dummy_matrix ] )) for plucker in gr_36_vertices]

gr_36_vertices_num_all=[np.linalg.det(np.array([ [row[int(i)-1] for i in list(plucker[1:])] for row in dummy_matrix ] )) for plucker in all_plucker_coords_gr36]

gr_36_vertices_num_all+=[(Symbol('a236')*Symbol('a145')-Symbol('a123')*Symbol('a456')).subs([z for z in zip(all_plucker_coords_gr36_s,gr_36_vertices_num_all)]), ((Symbol('a123')*Symbol('a156')*Symbol('a346')+Symbol('a126')*Symbol('a134')*Symbol('a356'))*(Symbol('a136')**(-1))).subs([z for z in zip(all_plucker_coords_gr36_s,gr_36_vertices_num_all)]) ]'''

new_variables_dictionary = {}

######################################################################

def overprint(s):
  stdout.write(f'\r{s}')
  stdout.flush()


def same_cluster_test(cluster1, cluster2):
  if len(cluster1) != len(cluster2):
    return False

  threshhold = len(cluster1)
  equalities = 0

  for i in range(len(cluster1)):
    for j in range(len(cluster2)):
      if isclose(cluster1[i], cluster2[j]):
        equalities += 1

  if threshhold == equalities:
    return True
  if threshhold < equalities:
    print(cluster1, cluster2)
    print('Threshhold < equalities')
    return True
  else:
    return False


def mutation(vertex, quiver):
  #Function that mutates a given quiver at the given vertex
  dummy_quiver = quiver.copy()
  for i in range(len(quiver)):
    for j in range(i):

      #Piecewise function as described in equation (1.1)
      if vertex == i or vertex == j:
        quiver[i][j] = -dummy_quiver[i][j]

      elif dummy_quiver[i][vertex] * dummy_quiver[vertex][j] <= 0:
        quiver[i][j] = dummy_quiver[i][j]

      else:
        quiver[i][j] = dummy_quiver[i][j] + abs(
          dummy_quiver[i][vertex]) * dummy_quiver[vertex][j]

      quiver[j][i] = -quiver[i][j]

  return quiver


def quiver_mutations(mutations, quiver, testmode=False):
  all_quivers = set()
  mutation_so_far = []
  mutation_count = 0
  for k in mutations:
    mutation_so_far += [k]
    #Update console
    mutation_count += 1
    overprint('Mutation: {0}/{1}'.format(mutation_count, (len(mutations))))

    #Mutates the quiver
    quiver = mutation(k, quiver)

    if testmode == True:
      hashable_quiver = tuple(
        tuple(quiver[i][j] for i in range(len(quiver)))
        for j in range(len(quiver)))
    else:
      hashable_quiver = (tuple(
        tuple(quiver[i][j] for i in range(len(quiver)))
        for j in range(len(quiver))), tuple(mutation_so_far))
    all_quivers.add(hashable_quiver)

  return all_quivers


def coordinates(mutations,
                quiver,
                vars=0,
                num_vars=0,
                mutables=0,
                xcoords=False,
                acoords=False,
                find_clusters=False,
                find_x1x=False,
                new_coordinate_number=0,
                testmode=False,
                muted=False):
  #This function returns a tuple of all A_coords and X_coords. Mutations are a list of vertices at which this function will perform mutations on the given quiver. The quiver's vertices are labeled by the variables in vars, which will be x0,...x(n-1) by default.

  #Currently there is nothing that differentiates between mutable and non-mutable vertices, and instead this differentiation is expected to be accounted for by the list of mutations inputed into the function.

  coordinatescomputetimetest = time()

  if testmode == True:
    mutations_so_far = []
    coords_with_mutations = []

  #Sets variables if none given
  if vars == 0:
    vars = ['x{0}'.format(i + 1) for i in range(len(quiver))]

  #turns variables in symbols to be used by Sympy
  vertices = [Symbol(vars[v]) for v in range(len(quiver))]

  #Set of all A-coordinates, starts with all the vertices.
  A_coords = {vertex for vertex in vertices}
  if testmode == True:
    coords_with_mutations += [(a, '') for a in A_coords]

  #Set of all X-coordinates. This block of code is designed to give all the X_coords generated from the first quiver.
  #I think this should be edited so that it only affects mutable vertices, but I'm not sure.
  if xcoords == True:
    X_coords = set()
    for k in range(mutables + 1):
      product_in = 1
      product_out = 1
      for i in range(len(quiver)):
        if quiver[i][k] > 0:
          product_in *= quiver[i][k] * vertices[i]
        elif quiver[i][k] < 0:
          product_out *= quiver[i][k] * vertices[i] * (-1)
      X_coords.add(product_out / product_in)

  if find_x1x == True:
    x1x_coords = set()

  if find_clusters == True:
    clusters = {frozenset([vertex for vertex in vertices])}

  mutation_count = 0
  #Performs the mutations.
  for mut in mutations:
    if testmode == True:
      mutations_so_far += [mut]
    #Update console
    mutation_count += 1
    if not muted:
      overprint('Mutation: {0}/{1}'.format(mutation_count, (len(mutations))))

    #Mutates the quiver
    quiver = mutation(mut, quiver)

    #Mutates the associated vertex
    product_in = 1
    product_out = 1
    for i in range(len(quiver)):
      if quiver[i][mut] > 0:
        product_in *= quiver[i][mut] * vertices[i]
      elif quiver[i][mut] < 0:
        product_out *= quiver[i][mut] * vertices[i] * (-1)

    #New A-coords from plucker list.
    new_a_coord = (vertices[mut]**(-1)) * ((product_in) + (product_out))

    newvars = [Symbol(x) for x in vars]

    evaluation_a_coord = new_a_coord.subs(list(zip(newvars, num_vars)))

    found_old = False
    for i in range(len(num_vars)):
      if isclose(evaluation_a_coord, num_vars[i]):
        new_a_coord = newvars[i]
        found_old = True

    if found_old == False:
      num_vars += [evaluation_a_coord]
      vars += [f'aI{new_coordinate_number}']
      new_variables_dictionary[f'aI{new_coordinate_number}'] = new_a_coord
      new_a_coord = Symbol(f'aI{new_coordinate_number}')

      new_coordinate_number += 1

    if find_x1x == True:
      x1x_coords.add(
        (product_out / product_in, (vertices[mut] * new_a_coord) / product_in))

    vertices[mut] = new_a_coord

    if testmode == True:
      coords_with_mutations.append(
        (simplify(vertices[mut]), mutations_so_far.copy()))

    #Adds the X coordinates and A coordinates of the mutated vertex.
    if xcoords == True:

      X_coords.add(product_out / product_in)

      for k in range(mutables + 1):
        product_in = 1
        product_out = 1
        for i in range(len(quiver)):
          if quiver[i][k] > 0:
            product_in *= quiver[i][k] * vertices[i]
          elif quiver[i][k] < 0:
            product_out *= quiver[i][k] * vertices[i] * (-1)
        X_coords.add(product_out / product_in)

    if acoords == True:
      A_coords.add(vertices[mut])

    if find_clusters == True:
      clusters.add(frozenset([simplify(vertex) for vertex in vertices]))

  if testmode == True:
    print('\nTestmode output, coords with mutations: ', coords_with_mutations)

  if xcoords == True and acoords == True:
    result = A_coords, X_coords
  elif acoords == True:
    result = A_coords
  elif xcoords == True:
    result = X_coords
  elif find_clusters == True:
    result = clusters
  elif find_x1x == True:
    result = x1x_coords

  if not muted:                
    overprint(
    f'\nTime to compute coordinates: {time()-coordinatescomputetimetest}\n')

  return result, new_coordinate_number, vars, num_vars, new_variables_dictionary


def main(user_input,
         num_of_mutations=1,
         quiver_data='',
         find_xcoords=False,
         find_acoords=False,
         find_clusters=False,
         find_x1x=False,
         new_coordinate_number=0,
         testmode=False,
         muted=False):

  quiver = quiver_data[0].copy()
  verts = quiver_data[1].copy()
  mutables = quiver_data[2]
  num_verts = quiver_data[3].copy()

  mutations_list = []
  lastint = -1
  for x in range(int(num_of_mutations)):
    randomint = randint(0, mutables)
    while randomint == lastint and mutables > 0:
      randomint = randint(0, mutables)
    mutations_list += [randomint]
    lastint = randomint

  if user_input == 'qu':

    return quiver_mutations(mutations_list, quiver, testmode)

  if user_input == 'gr':

    return coordinates(mutations_list, quiver, verts, num_verts, mutables,find_xcoords, find_acoords, find_clusters, find_x1x, new_coordinate_number, testmode, muted)


def numerical_test(iter,
                   quiver_data,
                   final_result,
                   do_find_clusters,
                   total_duplicates,
                   old_length,
                   skip_numeric_test_error,
                   find_coords,
                   keep_looping,
                   test_mode=False):

  vertices_numerical = quiver_data[3].copy()
  vertices = [Symbol(x) for x in quiver_data[1]]
  try:
    time_to_complete_numeric_test = time()
    final_result_list = list(x for x in final_result)
    to_remove = set()
    pairs_of_dupes_found = 0

    if not do_find_clusters:

      evals_plucker = [
        expr.subs([z for z in zip(vertices, vertices_numerical)])
        for expr in final_result_list
      ]

      overprint(f'Time to complete numeric test: ')
      for j in range(len(evals_plucker)):
        for i in range(j + 1, len(evals_plucker)):

          if j != i and isclose(evals_plucker[i], evals_plucker[j]) == True:
            pairs_of_dupes_found += 1
            to_remove.add(final_result_list[j])

    elif do_find_clusters:

      evals_plucker_clusters = [
        list(
          expr.subs([z for z in zip(vertices, vertices_numerical)])
          for expr in cluster) for cluster in final_result_list
      ]

      overprint('Time to complete numeric test: ')
      for j in range(len(evals_plucker_clusters)):
        for i in range(j + 1, len(evals_plucker_clusters)):
          if j != i and same_cluster_test(evals_plucker_clusters[j],
                                          evals_plucker_clusters[i]) == True:
            pairs_of_dupes_found += 1
            to_remove.add(final_result_list[j])

    else:
      print('\nError: Something weird happened with do_find_clusters. \n')
      raise KeyboardInterrupt

    #Testmode stuff
    if pairs_of_dupes_found == len(to_remove) and test_mode == True:
      print('Test dupes = amt removed passed')
    elif test_mode == True:
      print('Test dupes = amt rmoved failed', pairs_of_dupes_found,
            len(to_remove))
    ###

    total_duplicates += len(to_remove)

    final_result = final_result - to_remove

    overprint(
      f'Time to complete numeric test: {time()-time_to_complete_numeric_test}\n'
    )

  except Exception as e:
    if skip_numeric_test_error == False:
      print(
        '\n Error in numeric test. Maybe vertices are not plucker coordinates.  \n'
      )
      print('Error caused by: ', e)
      print()
      keep_looping = False
    pass

  overprint(
    f'\n{iter}: Found {-old_length+len(final_result)} new and {len(to_remove)} duplicate {find_coords.capitalize()}-coordinates. Total: {len(final_result)}\n'
  )

  last_coordinates_found = len(final_result) - old_length

  return final_result, last_coordinates_found, keep_looping

def time_conv(seconds):
  if seconds<60:
    return round(seconds,2), 'seconds'
  elif seconds<3600:
    return round(seconds/60,2), 'minutes'
  elif seconds<86400:
    return round(seconds/3600,2), 'hours'
  else:
    return round(seconds/86400,2), 'days'
