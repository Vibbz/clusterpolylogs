import numpy as np
from sympy import *
from random import randint
from random import random
from math import isclose
from sys import stdout
from time import time
from quivers import generate_quiver


gr_36_quiver=np.array([[0,1,1,-1,-1,0,0,0,0,0],[-1,0,0,1,0,1,-1,0,0,0],[-1,0,0,1,0,0,0,1,-1,0],[1,-1,-1,0,0,0,1,0,1,-1],[1,0,0,0,0,0,0,0,0,0],[0,-1,0,0,0,0,0,0,0,0],[0,1,0,-1,0,0,0,0,0,0],[0,0,-1,0,0,0,0,0,0,0],[0,0,1,-1,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0,0] ])
gr_36_vertices=['a236','a235','a136','a356','a123','a234','a345','a126','a156','a456']
gr_36_vertices_s=[Symbol(x) for x in gr_36_vertices]
dummy_matrix=np.array([[n+6*i+random() for n in range(6)] for i in range(3)])


#These are the associated nubmers to each vertex as listed above
gr_36_vertices_num=[np.linalg.det(np.array([ [row[int(i)-1] for i in list(plucker[1:])] for row in dummy_matrix ] )) for plucker in gr_36_vertices]

def mutation(vertex,quiver):
#Function that mutates a given quiver at the given vertex
  dummy_quiver=quiver.copy()
  for i in range(len(quiver)):
    for j in range(i):
      
      #Piecewise function as described in equation (1.1)
      if vertex == i or vertex == j:
        quiver[i][j]=-dummy_quiver[i][j]

      elif dummy_quiver[i][vertex]*dummy_quiver[vertex][j]<=0:
        quiver[i][j]=dummy_quiver[i][j]

      else:
        quiver[i][j]=dummy_quiver[i][j]+abs(dummy_quiver[i][vertex])*dummy_quiver[vertex][j]

      quiver[j][i]=-quiver[i][j]

  return quiver

def quiver_mutations(mutations,quiver,testmode=False):
  all_quivers=set()
  mutation_so_far=[]
  mutation_count=0
  for k in mutations:
    mutation_so_far+=[k]
    #Update console
    mutation_count+=1
    sys.stdout.write('\r Mutation: {0}/{1}'.format(mutation_count,(len(mutations))))
    sys.stdout.flush()
    
    #Mutates the quiver
    quiver=mutation(k,quiver)

    if testmode==True:
      hashable_quiver=tuple(tuple(quiver[i][j] for i in range(len(quiver))) for j in range(len(quiver)))
    else:
      hashable_quiver=(tuple(tuple(quiver[i][j] for i in range(len(quiver))) for j in range(len(quiver))), tuple(mutation_so_far))
    all_quivers.add(hashable_quiver)

  return all_quivers

def coordinates(mutations,quiver,vars=0,mutables=0,xcoords=0,acoords=0,clusters=0,testmode=False):
  #This function returns a tuple of all A_coords and X_coords. Mutations are a list of vertices at which this function will perform mutations on the given quiver. The quiver's vertices are labeled by the variables in vars, which will be x0,...x(n-1) by default.
  
  #Currently there is nothing that differentiates between mutable and non-mutable vertices, and instead this differentiation is expected to be accounted for by the list of mutations inputed into the function.

  if testmode==True:
    mutations_so_far=[]
    coords_with_mutations=[]
  
  #Sets variables if none given
  if vars==0:
    vars=['x{0}'.format(i+1) for i in range(len(quiver))]
  
  #turns variables in symbols to be used by Sympy
  vertices=[Symbol(vars[v]) for v in range(len(quiver))]

  #Set of all A-coordinates, starts with all the vertices.
  A_coords={vertex for vertex in vertices}
  if testmode==True:
    coords_with_mutations+=[(a,'') for a in A_coords]
  
  #Set of all X-coordinates. This block of code is designed to give all the X_coords generated from the first quiver. 
  #I think this should be edited so that it only affects mutable vertices, but I'm not sure.
  if xcoords==True:
    X_coords=set()
    for k in range(mutables+1):
      product_in=1
      product_out=1
      for i in range(len(quiver)):
        if quiver[i][k]>0:
          product_in*=quiver[i][k]*vertices[i]
        elif quiver[i][k]<0:
          product_out*=quiver[i][k]*vertices[i]*(-1)
      X_coords.add(product_out/product_in)
  
  mutation_count=0
  #Performs the mutations.
  for mut in mutations:
    if testmode==True:
      mutations_so_far+=[mut]
    #Update console
    mutation_count+=1
    sys.stdout.write('\r Mutation: {0}/{1}'.format(mutation_count,(len(mutations))))
    sys.stdout.flush()
    
    #Mutates the quiver
    quiver=mutation(mut,quiver)

    #Mutates the associated vertex
    product_in=1
    product_out=1
    for i in range(len(quiver)):
      if quiver[i][mut]>0:
        product_in*=quiver[i][mut]*vertices[i]
      elif quiver[i][mut]<0:
        product_out*=quiver[i][mut]*vertices[i]*(-1)

    vertices[mut]=(vertices[mut]**(-1))*((product_in)+(product_out)) 
    
    if testmode==True:
        coords_with_mutations.append((simplify(vertices[mut]), mutations_so_far.copy()))
    
    #Adds the X coordinates and A coordinates of the mutated vertex.
    if xcoords==True:
      X_coords.add(simplify(product_out/product_in))
    if acoords==True:
      A_coords.add(simplify(vertices[mut]))

  if testmode==True:
    print('\n Testmode output, coords with mutations: ', coords_with_mutations)
  
  if xcoords==True and acoords==True:
    result=A_coords, X_coords
  elif acoords==True:
    result=A_coords
  elif xcoords==True:
    result=X_coords
  
  return result

def main(user_input,num_of_mutations=1,quiver_data='',find_xcoords=False,find_acoords=False,find_clusters=False,testmode=False):


  quiver=quiver_data[0].copy()
  verts=quiver_data[1].copy()
  mutables=quiver_data[2]

  #print('\n\n\n Testing \n\n\n ')
  #print(quiver,verts,mutables)
  #print(gr_36_quiver,gr_36_vertices,3)
  #Create quiver data
  if quiver_data=='gr36':
    quiver=np.array([[0,1,1,-1,-1,0,0,0,0,0],[-1,0,0,1,0,1,-1,0,0,0],[-1,0,0,1,0,0,0,1,-1,0],[1,-1,-1,0,0,0,1,0,1,-1],[1,0,0,0,0,0,0,0,0,0],[0,-1,0,0,0,0,0,0,0,0],[0,1,0,-1,0,0,0,0,0,0],[0,0,-1,0,0,0,0,0,0,0],[0,0,1,-1,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0,0] ])
    verts=['a236','a235','a136','a356','a123','a234','a345','a126','a156','a456']
    mutables=3

  elif quiver_data=='ex2':
    quiver=np.array([[0,1],[-1,0]]) #This is the basic two-vertex quiver with one edge going from one vertex to the other: 1->2
    mutables=1
    verts=0

  elif quiver_data=='ex3':
   quiver=np.array([[0,1,0],[-1,0,1],[0,-1,0]] ) #This is the basic three-vertex quiver with two edges, 1 -> 2 -> 3
   mutables=2
   verts=0

  if testmode==True:
    test_quiver=generate_quiver('gr36')[0]
    print('Testing generate_quiver: ', test_quiver==quiver)
  
  mutations_list=[]
  lastint=-1
  for x in range(int(num_of_mutations)):
    randomint=randint(0,mutables)
    while randomint==lastint:
      randomint=randint(0,mutables)
    mutations_list+=[randomint]
    lastint=randomint
  
  if user_input=='qu':
    
    return quiver_mutations(mutations_list,quiver,testmode)

  if user_input=='gr':
    #print(mutations_list,quiver,verts,find_xcoords,find_acoords,find_clusters,testmode)
    return coordinates(mutations_list,quiver,verts,mutables,find_xcoords,find_acoords,find_clusters,testmode)
