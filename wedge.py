import itertools
from sympy import Symbol, Add, Mul
from sympy.abc import x,y,z,a,b,c,d,e
import scipy as sp
import numpy as np
from sympy import Matrix
from functions import overprint, time_conv
from time import time
import re
from time import time

from wolframclient.evaluation import WolframLanguageSession
from wolframclient.language import wl, wlexpr

a12=Symbol('a12')
a13=Symbol('a13')
a14=Symbol('a14')
a15=Symbol('a15')
a16=Symbol('a16')
a26=Symbol('a26')
a36=Symbol('a36')
a46=Symbol('a46')
a56=Symbol('a56')
a23=Symbol('a23')
a24=Symbol('a24')
a25=Symbol('a25')
a34=Symbol('a34')
a35=Symbol('a34')
a45=Symbol('a45')

def mathematica_test():
  try:
    session=WolframLanguageSession()
    test=session.evaluate(wlexpr('Range[5]'))
    session.terminate()
    if (test==[1,2,3,4,5]).all():
      print('Mathematica successfully loaded.')
  except Exception as exception:
    print(exception)
    print('Mathematica failed to load.')
    
  

def matrix_to_mathematica(matrix):
  generating_string='{'
  for v in matrix:
    vstring='{'
    for x in v:
      vstring+=str(x)
      vstring+=','
    vstring=vstring[:-1]
    vstring+='}'
    generating_string+=vstring + ','

  if generating_string[-1]==',':
    generating_string=generating_string[:-1]
  generating_string+='}'

  return generating_string


def mta(expression):
  result=[]
  if len(expression.args)==0:
    return expression
  else:
    for i in range(len(expression.args)):
      try:
        result+=[expression.args[i].args[0]*expression.args[i].args[1]]
      except:
        result+=[expression.args[i]]
    
    return Add(*tuple(result))
  

def sym_wedge(v1,v2,variables,wedge_vars=[]):

  v1={var:mta(v1).coeff(var) for var in variables}
  
  v2={var:mta(v2).coeff(var) for var in variables}



  result=0
  combos=list(itertools.combinations(variables,2))
  for c in combos:
    result+=(v1[c[0]]*v2[c[1]]-v1[c[1]]*v2[c[0]])*(Symbol(str(c[0])+'^'+str(c[1])))

  
  if len(wedge_vars)==0:
    wedge_vars=[Symbol(str(c[0])+'^'+str(c[1])) for c in itertools.combinations(variables,2)]

  return result,tuple(result.coeff(va) for va in wedge_vars)


def basis_of_CL2(list_x1x,vars):
  #Input should be a list of tuples containg X and 1+X, X being X-coords. This is the output of the program when entering x1x for the type of coordinate.
  computation_time=time()
  overprint('Creating the generating set...')

  basis_variables=[Symbol(str(c[0])+'^'+str(c[1])) for c in itertools.combinations(vars,2)]
 

  generating_set=[sym_wedge(it[0],it[1],vars,basis_variables) for it in list_x1x]

  
  generating_vectors=np.transpose(np.array([ list(g[1][i] for g in generating_set) for i in range(len(generating_set[0][1])) ]).astype(int))


  generating_string=matrix_to_mathematica(generating_vectors)

  generating_string=f'LatticeReduce[{generating_string}]'

  overprint(f'Completed generating set at {time_conv(time()-computation_time)[0]} {time_conv(time()-computation_time)[1]}, computing the basis...' + '\n')
 
  session=WolframLanguageSession()
 
  basis=np.transpose(np.array(session.evaluate(wlexpr(generating_string))))
  
  session.terminate()

  dimension=basis.shape[1]


  overprint(f'Completed basis at {time_conv(time()-computation_time)[0]} {time_conv(time()-computation_time)[1]}, computing the symbolic basis...' + '\n')
  
  symbolic_basis=[]
  for vec in np.transpose(basis):
    basis_elem=0
    vindex=0
    for j in vec:
      basis_elem+=j*basis_variables[vindex]
      vindex+=1

    symbolic_basis+=[basis_elem]

  if dimension!=len(symbolic_basis):
    print('\n Dimension error when computing CL2 \n')
  
  print(generating_vectors.shape,basis.shape)
  overprint(' '*100)
  overprint(f'Computed basis in {time_conv(time()-computation_time)[0]} {time_conv(time()-computation_time)[1]}.')
  return basis,symbolic_basis,dimension


def cl_intb_map(tensor,variables):
  if len(tensor)<=2:
    print('Tuple length is {0}'.format(len(tensor)) + ', which is too short.')
  else:
    result=[]
    for i in range(len(tensor)-1):
      list_tensor=list(tensor)
      list_tensor[i]=sym_wedge(tensor[i],tensor[i+1],variables)[0]
      list_tensor.pop(i+1)
      result+=[tuple(list_tensor)]
    return tuple(result) 


def expected_dimension_of_cln_gr2m(n,m):
  sum=0
  for k in range(3,n+2):
    sum+=sp.special.comb(m-1,k)
  return int(sum)


def basis_of_CLn(clusters,variables,cl2,n=3):
  basis_variables=[Symbol(str(c[0])+'^'+str(c[1])) for c in itertools.combinations(variables,2)]

  compute_time=time()
  
  if n<3:
    print('Error: n<3 for CLn')
    return None
  print('n=3')
  basis_An=set()

  iter=0
  for c in clusters:
    iter+=1
    c=list(c)
    for p in itertools.product(c,repeat=n):
      basis_An.add(p)

  
  basis_An=list(basis_An)

  

  print('\n',f'Legnth of A_n basis: {len(basis_An)}','\n')

  overprint('Creating image_of_An')
  image_of_An=[cl_intb_map(x,variables) for x in basis_An]

  overprint(f'Created image_of_An at {time_conv(-compute_time+time())[0]} {time_conv(-compute_time+time())[1]}' + '      \n')

  print('\n',f'Length of image_of_An: {len(image_of_An)}','\n')


  iter=0
  mat=[]
  for x in image_of_An:
    iter+=1
    overprint(f'Adding to matrix from image_of_An: {iter}/{len(image_of_An)}')
    v=[]
    for i in range(n-1):
      if i==0:
        j=1
      if i==1:
        j=0
      for bvar in basis_variables:
        for var in variables:
          v+=[x[i][i].coeff(bvar)*x[i][j].coeff(var)]
    mat+=[v]

  overprint(f'Created A_n image matrix at {time_conv(-compute_time+time())[0]} {time_conv(-compute_time+time())[1]}' '    \n')


  basis_of_W=[]
  doubles=[]
  for w in cl2:
    for a in variables:
      doubles+=[(w,Symbol(a))]
  for d in doubles:
    basis_of_W+=[(d,(0,0))]+[((0,0),(d[1],d[0]))]

  print('\n',f'Length of basis_of_W: {len(basis_of_W)}','\n')

  iter=0
  for x in basis_of_W:
    iter+=1
    overprint(f'Adding to matrix from W: {iter}/{len(basis_of_W)}')
    v=[]
    for i in range(n-1):
      if i==0:
        j=1
      if i==1:
        j=0
      for bvar in basis_variables:
        for var in variables:
          try:
            v+=[x[i][i].coeff(bvar)*x[i][j].coeff(var)]
          except AttributeError:
            v+=[0]
    mat+=[v]

  mat=np.transpose(np.array(mat))

  overprint(f'Finished creating matrix at {time_conv(-compute_time+time())[0]} {time_conv(-compute_time+time())[1]}'+'        \n')

  try:
    print('Size: ', mat.shape)
  except Exception as e:
    print(e)
    pass

  try:
    with open('matrix.txt','w') as file:
      #print(str(sp.sparse.csr_array(np.array(mat).astype(int))))
      file.write(str(np.array(mat)))
      file.close()
  except Exception as e:
    print(e)
    pass

  print(f'Computing null space...')

  generating_string=matrix_to_mathematica(mat)
  generating_string=f'NullSpace[{generating_string}]'

  
  session=WolframLanguageSession()

  null_space=np.array(session.evaluate(wlexpr(generating_string)))

  session.terminate()


  try:

    with open('matrix2.txt','w') as file:
      file.write(str(sp.sparse.csr_array(null_space)))
      file.close()
  except Exception as e:
    print(e)
    pass

  '''
  loops=0
  symb=[]
  for vect in null_space:
    loops+=1
    index=0
    print(vect)
    print()
    part1=0
    part2=0
    part3=0
    part4=0
    for i in range(n-1):
      for j in range(len(basis_variables)):
        for k in range(len(variables)):
          if (i*len(basis_variables)*len(variables)+j*len(variables)+k-index) not in [-1,0,1]:
            print('ijk: ',i,j,k, 'index: ',i*len(basis_variables)*len(variables)+j*len(variables)+k)
          index=i*len(basis_variables)*len(variables)+j*len(variables)+k
          coeff=vect[index]
          if i==0:
            if coeff!=0:
              part1+=basis_variables[j]
              part2+=coeff*Symbol(variables[k])
              print('1,2',part1,part2)
              print('i,j,k,index,coeff: ',i,j,k,index,coeff)
          else:
            if coeff!=0:
              part3+=coeff*Symbol(variables[i])
              part4+=basis_variables[j]
              print('3,4:',part3,part4)
    symb+=[ ((part1,part2),  (part3,part4)) ]
    print(symb)
          
  print('\nsymb:\n')
  print(symb)
  '''

  

  print('Null space shape:', null_space.shape)

  null_space=list(null_space)

  num=0
  to_remove=[]
  for vec in null_space:
      #print(vec,num)
      #print('range:', list(range(len(image_of_An),len(image_of_An)+len(basis_of_W))))
      #print([vec[j] for j in range(len(image_of_An),len(image_of_An)+len(basis_of_W))])
      #print([0 for j in range(len(image_of_An),len(image_of_An)+len(basis_of_W))])

      if ( [vec[j] for j in range(len(image_of_An),len(image_of_An)+len(basis_of_W))] == [0 for j in range(len(image_of_An),len(image_of_An)+len(basis_of_W))] ):
        #print('NS:', [vec[j] for j in range(len(image_of_An),len(image_of_An)+len(basis_of_W))])
        to_remove+=[num]
      num+=1

  removed=0
  for num in to_remove:
    null_space=np.delete(null_space,num-removed,0)
    removed+=1

  preimage=[]
  for vec in null_space:
    vec=vec[:len(image_of_An)]
    pre_vec=0
    terms=0
    for j in range(len(vec)):
      pre_vec+=vec[j]*Symbol(str(basis_An[j]))
      if vec[j]!=0:
        terms+=1
    
    print(pre_vec,terms,'\n')
    print()
    preimage+=[pre_vec]




  print(f'Done at {time_conv(-compute_time+time())[0]} {time_conv(-compute_time+time())[1]}')

  return preimage,null_space.shape[0]

