import numpy as np
from sympy import *
from random import randint
from random import random
from math import isclose
import sys
from time import time
from time import sleep
from functions import mutation, quiver_mutations, coordinates, main
from quivers import generate_quiver

#Here we consider quivers as nxn skew-symmetric matrices where the entries e_ij correspond to the number of edges going from vertex i into vertex j, with the number considered negative if the direction is reversed. Mutation is described on quivers in this manner in equation (1.1) at the bottom of the first page on Zickert's notes on Cluster Algebras.


#Defining gr(3,6) and codifying plucker relations

gr_36_quiver=np.array([[0,1,1,-1,-1,0,0,0,0,0],[-1,0,0,1,0,1,-1,0,0,0],[-1,0,0,1,0,0,0,1,-1,0],[1,-1,-1,0,0,0,1,0,1,-1],[1,0,0,0,0,0,0,0,0,0],[0,-1,0,0,0,0,0,0,0,0],[0,1,0,-1,0,0,0,0,0,0],[0,0,-1,0,0,0,0,0,0,0],[0,0,1,-1,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0,0] ])
gr_36_vertices=['a236','a235','a136','a356','a123','a234','a345','a126','a156','a456']
gr_36_vertices_s=[Symbol(x) for x in gr_36_vertices]

#This is the matrix from which the plucker relations will be evaluated
#dummy_matrix=np.array([[1,2,3,4,5,6],[7,8,9,10,11,12],[13,14,16,17,18,19]])

dummy_matrix=np.array([[n+6*i+random() for n in range(6)] for i in range(3)])


#These are the associated nubmers to each vertex as listed above
gr_36_vertices_num=[np.linalg.det(np.array([ [row[int(i)-1] for i in list(plucker[1:])] for row in dummy_matrix ] )) for plucker in gr_36_vertices]


#Later, I'll make code that actually just generates this stuff in general.


# Below is the main program



###########################
# Main program code block #
###########################
print('''This program works by doing some number of random mutations repeatedly for some number of iterations. The program can be interrupted at any time with keyboard interrupt to print the result so far.  \n
Enter "co" to compute coordinates. \n
Enter "qu" to compute quivers resulting from mutations. \n
Enter "test" to enter test mode (not useful). \n''')

while True:
  option1=''
  test_mode=False
  while option1 not in ['gr','qu']:
    option1=input('Option: ')
    if option1=='test':
      test_mode=True
      print('Test mode enabled.')
    if option1=='co':
      option1='gr'
  
  quiver_data_s=input('Enter the quiver to use.\n gr36 uses the quiver for gr(3,6) with plucker coordinate vertices. \n gr2n uses the quiver for gr(2,n) with plucker coordinate vertices. \n ex2 uses the quiver 1->2 with default vertex labeling. \n ex3 uses the quiver 1->2->3 with default vertex labeling. \n Enter: ')
  quiver_data=generate_quiver(quiver_data_s)
  
  if option1=='gr':
    find_coords=input('Enter "x" to find X-coords, "a" to find A-coords, or "cl" to find all the clusters: ')
    if find_coords=="x":
      do_find_x_coords=1
      do_find_a_coords=0
    else:
      do_find_a_coords=1
      do_find_x_coords=0
    if find_coords=="cl":
      do_find_clusters=1
    else:
      do_find_clusters=0
  else:
    do_find_a_coords=0
    do_find_x_coords=0
  
  num_of_mutations=int(input('Enter the number of random mutations (not recommended to be larger than 20): '))
  iteration_limit=input('Enter the iteration limit or leave blank or enter 0 to iterate forever: ')
  
  if iteration_limit=='':
    iteration_limit=0
  iteration_limit=int(iteration_limit)
  
  condition=True
  
  universal_start_time=time()
  
  final_result=set()
  iter=0
  
  try:
    while condition==True:
      iter+=1

      if test_mode==True:
        print('\n General test mode output...')
        print('Option 1: {0}'.format(option1))
        print('Number of mutations: {0}'.format(num_of_mutations))
        print('Quiver: ',quiver_data[0])
        print('Vertices: ',quiver_data[1])
        print('Mutables: ',quiver_data[2])
        print('Numerical vertices: ', quiver_data[3])
        print()

      
      ####################### 
      #Actually calling main#
      new_results=main(option1,num_of_mutations,[quiver_data[0].copy(),quiver_data[1].copy(),quiver_data[2],quiver_data[3].copy()],do_find_x_coords,do_find_a_coords,do_find_clusters,test_mode)
      
      old_length=len(final_result)
      
      for co in new_results:
        final_result.add(co)
      
      if option1=='gr':
        #Numeric test for equality modulo plucker relations
        #if statement is standin
        
        vertices_numerical=quiver_data[3].copy()
        vertices=[Symbol(x) for x in quiver_data[1]]

        
        if True==True:
          try:
            to_remove=[]
            pairs_of_dupes_found=0
            
            evals_plucker=[expr.subs([z for z in zip(vertices,vertices_numerical)]) for expr in final_result]
            
            for j in range(len(evals_plucker)):
              for i in range(j,len(evals_plucker)):
                  if j!=i:
                    if isclose(evals_plucker[i],evals_plucker[j]):
                      pairs_of_dupes_found+=1
  
                      #Testmode information
                      #if test_mode==True:
                      #  print('\n Test mode output, finding duplicates: \n', i, j, list(final_result)[i],list(final_result)[j], '\n', list(final_result)[i].subs([z for z in zip(vertices,vertices_numerical)]), list(final_result)[j].subs([z for z in zip(vertices,vertices_numerical)]),evals_plucker[i],evals_plucker[j])
                      #####################
                      to_remove+=[list(final_result)[j]]
  
            if pairs_of_dupes_found==len(to_remove) and test_mode==True:
              print('Test dupes = amt removed passed')
            elif test_mode==True:
              print('Test dupes = amt rmoved failed', pairs_of_dupes_found, len(to_remove))

            for x in to_remove:
              final_result.remove(x)
            final_result=set(final_result)
          except:
            print('\n Error in numeric test, or vertices are not plucker coordinates.  \n')
            pass
            
          sys.stdout.write('\r \n {1}: Removed {4} duplicates. Found {3} new {2}-coordinates. Found {0} total {2}-coordinates. \n'.format(len(final_result),iter,find_coords.capitalize(),-old_length+len(final_result),len(to_remove)))
          sys.stdout.flush()
        else:
          sys.stdout.write('\r \n {1}: Found {3} new {2}-coordinates, {0} total {2}-coordinates. \n'.format(len(final_result),iter,find_coords.capitalize(),-old_length+len(final_result)))
          sys.stdout.flush()
          
        if option1=='qu':
          sys.stdout.write('\r \n Found {0} quivers on iteration {1}. \n'.format(len(final_result),iter))
          sys.stdout.flush()
  
      
      
      if iteration_limit!=0:
        condition= iter<iteration_limit
  except KeyboardInterrupt:
    pass
  
  
  #Test sameness of coordinates
  if option1=='gr':  
    option3=input('Enter "y" to test equivalency among {0}-coordinates symbolically: '.format(find_coords))
    if option3=='y':
          final_result=list(final_result)
          #print(final_result[len(final_result)-1])
          #Symbolic test for sameness
          to_remove=set()
          for i in range(len(final_result)-1):
            for j in range(i+1,len(final_result)-1):
              sys.stdout.write('\r Checking equivalence of {0},{1}    '.format(i,j))
              sys.stdout.flush()   
              same_test=simplify(final_result[i]-final_result[j])
              if same_test==0:
                print('Found equivalence at {0},{1}'.format(i,j))
                if i!=j:
                  print('Found equivalence between expr at position {0} and expr at position {1}'.format(i,j))
                  to_remove.add(i)
          tempvar=0
          for index in to_remove:
            final_result.pop(index-tempvar)
            tempvar+=1
          final_result=set(final_result) 
          print(final_result)
          #print('Time for iteration:', time()-start_time)
          sys.stdout.write('\r'+' '*20)

  
  print('\n \n Final result: \n', final_result, '\n Number found: ',len(final_result))
  print('\n Total computation time: ', time()-universal_start_time)
  
  
  
  if option1=='qu':
    if test_mode==True:
      test='Passed'
      print('Test skew-symmetry: ')
      e=set()
      for x in final_result:
        x=list(list(t) for t in x)
        y=x.pop(-1)
        x=np.asarray(x[0])
        skewtranspose=-np.transpose(x)
        comparison=x==skewtranspose
        if comparison.all()==False:
          print(False)
          test='Failed'
        print(x,'\n','Mutations: ',y)
      print('\n',test)
  
  print('\n Run again:')

