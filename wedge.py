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
import pickle

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
  variables=[str(v) for v in variables]
  variables.sort()
  variables=[Symbol(v) for v in variables]

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
  vars.sort()
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
  variables=[str(v) for v in variables]
  variables.sort()
  variables=[Symbol(v) for v in variables]

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

def null_space(matrix):
  
  generating_string=matrix_to_mathematica(matrix)
  generating_string=f'NullSpace[{generating_string}]'

  
  session=WolframLanguageSession()

  nspace_sss=np.array(session.evaluate(wlexpr(generating_string)))

  session.terminate()

  return nspace_sss

def lattice_reduce(matrix):

  generating_string=matrix_to_mathematica(matrix)
  generating_string=f'LatticeReduce[{generating_string}]'

  
  session=WolframLanguageSession()

  lreduce=np.array(session.evaluate(wlexpr(generating_string)))

  session.terminate()

  return lreduce


def basis_of_CLn(clusters,variables,cl2,n=3):
  variables.sort()
  basis_variables=[Symbol(str(c[0])+'^'+str(c[1])) for c in itertools.combinations(variables,2)]

  compute_time=time()

  print(variables)
  print(basis_variables)
  
  if n<3:
    print('Error: n<3 for CLn')
    return None
  print('n=3')


  basis_An=[]
  for c in clusters:
    c=list(str(x) for x in c)
    for p in itertools.product(c,repeat=n):
      basis_An+=[p]

  basis_An.sort()

  basis_An=[(Symbol(x[0]),Symbol(x[1]),Symbol(x[2]))for x in basis_An]

  print(type(basis_An),type(basis_An[0]),type(basis_An[0][0]))

  to_remove=[]  
  for i in range(len(basis_An)):
    if basis_An[len(basis_An)-i-1] in basis_An[0:len(basis_An)-i-1]:
        to_remove+=[len(basis_An)-i-1]

  for ind in to_remove:
    basis_An.pop(ind)

 

  print('\n',f'Legnth of A_n basis: {len(basis_An)}','\n')

  overprint('Creating image_of_An')
  image_of_An=[cl_intb_map(x,variables) for x in basis_An]


  #print(image_of_An[0],image_of_An[1],image_of_An[2])

  overprint(f'Created image_of_An at {time_conv(-compute_time+time())[0]} {time_conv(-compute_time+time())[1]}' + '      \n')

  print('\n',f'Length of image_of_An: {len(image_of_An)}','\n')

  try:
    #with open('temp_An_matrix.pk','rb') as fi:
    #  mat=pickle.load(fi)
      print(1+'s')
  except:

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
      
    try:
      #with open('temp_An_matrix.pk','wb') as fi:
      #  pickle.dump(mat,fi)
        print(1+'s')
    except:
      pass

  print(np.count_nonzero(mat))

  basis_of_W=[]
  doubles=[]
  for w in cl2:
    for a in variables:
      doubles+=[(w,Symbol(a))]
  for d in doubles:
    basis_of_W+=[(d,(0,0))]+[((0,0),(d[1],d[0]))]

  #print(basis_of_W[11])

  print('\n',f'Length of basis_of_W: {len(basis_of_W)}','\n')

  print(basis_of_W)

  iter=0
  for x in basis_of_W:
    iter+=1
    #print('\n',i,iter,'\n')
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
            #print(i,j,x,bvar,var,x[i][i].coeff(bvar)*x[i][j].coeff(var),x[i][i].coeff(bvar),x[i][j].coeff(var))
          except AttributeError:
            if x[i][i]!=0:
              print('x[i][i]',x[i][i])
            if x[i][j]!=0:
              print('x[i][j]',x[i][j])
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

  nspace=np.array(session.evaluate(wlexpr(generating_string)))

  session.terminate()

  
  #new_m=np.array([vec[186:] for vec in null_space])

  #print(new_m,new_m.shape)

  #print(np.linalg.matrix_rank(new_m)
  #print(np.linalg.matrix_rank(np.transpose(new_m)))


  try:
    with open('matrix2.txt','w') as file:
      file.write(str(sp.sparse.csr_array(nspace)))
      file.close()
  except Exception as e:
    print(e)
    pass

  

  print('Null space shape:', nspace.shape)

  nspace=list(nspace)

  num=0
  to_remove=[]
  for vec in nspace:
      if ( [vec[j] for j in range(len(image_of_An),len(image_of_An)+len(basis_of_W))] == [0 for j in range(len(image_of_An),len(image_of_An)+len(basis_of_W))] ):
        to_remove+=[num]
      num+=1

  removed=0
  for num in to_remove:
    nspace=np.delete(nspace,num-removed,0)
    removed+=1

  #mat2=np.array([v[:185] for v in mat])
  #print(mat2.shape)
  #mat3=np.transpose(np.array([v[:185] for v in nspace]))
  #print(mat3.shape)
  '''

  for vector in nspace:
    expr=0
    terms=0
    v_dic={}
    for k in range(n-1):
      for i in range(len(basis_variables)):
        for j in range(len(variables)):
          if k==1:
            expr+=vector[k*len(basis_variables)*len(variables) + i*len(variables) + j]*Symbol(str((basis_variables[i],Symbol(str(variables[j])))) )
            if (vector[k*len(basis_variables)*len(variables) + i*len(variables) + j]*Symbol(str((basis_variables[i],Symbol(str(variables[j])))) ))!=0:
                terms+=1
                try:
                  v_dic[str(variables[j])+',r']+=vector[k*len(basis_variables)*len(variables) + i*len(variables) + j]*Symbol(str(basis_variables[i]))
                except:
                  v_dic[str(variables[j])+',r']=vector[k*len(basis_variables)*len(variables) + i*len(variables) + j]*Symbol(str(basis_variables[i]))
          if k==0:
            expr+=vector[k*len(basis_variables)*len(variables) + i*len(variables)+j]*Symbol(str((Symbol(str(variables[j])),basis_variables[i])) )
            if (vector[k*len(basis_variables)*len(variables) + i*len(variables) + j]*Symbol(str((basis_variables[i],Symbol(str(variables[j])))) ))!=0:
                terms+=1
            try:
              v_dic[str(variables[j])+',l']+=vector[k*len(basis_variables)*len(variables) + i*len(variables) + j]*Symbol(str(basis_variables[i]))
            except:
              v_dic[str(variables[j])+',l']=vector[k*len(basis_variables)*len(variables) + i*len(variables) + j]*Symbol(str(basis_variables[i]))
    print('\n',expr,'\n\n',v_dic,terms,'\n')

 
  test_preimage=[]
  for vec in null_space:
    vec=vec[:len(image_of_An)]
    pre_vec=0
    terms=0
    vector_dictionary={}
    for j in range(len(vec)):
      clintbtensor=cl_intb_map(basis_An[j],variables)
      try:
        vector_dictionary[str(clintbtensor[0][1])+',r']=(vector_dictionary[str(clintbtensor[0][1])+',r']+clintbtensor[0][0])*vec[j]
      except Exception as exception:
        vector_dictionary[str(clintbtensor[0][1])+',r']=clintbtensor[0][0]*vec[j]
      try:
        vector_dictionary[str(clintbtensor[1][0])+',l']=(vector_dictionary[str(clintbtensor[1][0])+',l']+clintbtensor[1][1])*vec[j]
      except Exception as exception:
        vector_dictionary[str(clintbtensor[1][0])+',l']=clintbtensor[1][1]*vec[j]

      pre_vec+=vec[j]*Symbol(str(cl_intb_map(basis_An[j],variables)))
      if vec[j]!=0:
        terms+=1
    
    transformed_vector=[]
    for var in variables:
      if vector_dictionary[var+',r']!=0:
        transformed_vector+=[(Symbol(var),vector_dictionary[var+',r'])]
      if vector_dictionary[var+',l']!=0:
        transformed_vector+=[((vector_dictionary[var+',l'],Symbol(var)),Symbol(var))]
    print(transformed_vector,terms)
  '''
  print()
  

  print('\n')

  preimage=[]
  term_list=[]
  for vec in nspace:
    vec=vec[:len(image_of_An)]
    pre_vec=0
    terms=0
    for j in range(len(vec)):
      pre_vec+=vec[j]*Symbol(str(basis_An[j]))
      if vec[j]!=0:
        terms+=1
    
    print(pre_vec,terms,'\n')
    print()
    term_list+=[terms]
    preimage+=[pre_vec]

  print()
  print(term_list)

  print(f'Done at {time_conv(-compute_time+time())[0]} {time_conv(-compute_time+time())[1]}')

  return preimage,nspace.shape[0]

#e1=cl_intb_map((a,b,b),[a,b])
#e2=cl_intb_map((b,a,b),[a,b])
#e3=cl_intb_map((b,b,a),[a,b])

#print(e1,e2,e3)

def new_basis_CLn(clusters,variables,cl2,n=3):
  variables.sort()
  basis_variables=[Symbol(str(c[0])+'^'+str(c[1])) for c in itertools.combinations(variables,2)]

  compute_time=time()

  print(variables)
  print(basis_variables)
  
  if n<3:
    print('Error: n<3 for CLn')
    return None
  print('n=3')

  basis_An=[]
  for c in clusters:
    c=list(str(x) for x in c)
    for p in itertools.product(c,repeat=n):
      basis_An+=[p]

  basis_An.sort()

  basis_An=[(Symbol(x[0]),Symbol(x[1]),Symbol(x[2])) for x in basis_An]


  to_remove=[]  
  for i in range(len(basis_An)):
    if basis_An[len(basis_An)-i-1] in basis_An[0:len(basis_An)-i-1]:
        to_remove+=[len(basis_An)-i-1]

  for ind in to_remove:
    basis_An.pop(ind)


  print(len(basis_An))

  try:
    with open('temp_An_matrix.pk','rb') as fi:
      mat=pickle.load(fi)
      print('loaded matrix, shape, nonzero:', mat.shape, np.count_nonzero(mat))
  except:

    index=0
    mat=np.zeros((len(variables)*len(basis_variables)*2,len(basis_An)+len(variables)*len(cl2)*2)).astype(int)

    print('mat shape:', mat.shape)


    for k in range(len(basis_An)):
      x=cl_intb_map(basis_An[k],variables)
      overprint(f'Adding to matrix from image_of_An: {k+1}/{len(basis_An)}')
      for i in range(len(basis_variables)):
        for j in range(len(variables)):
          coeff=(x[0][0].coeff(basis_variables[i])) * (x[0][1].coeff(variables[j]))
          index=i*len(variables) + j
          mat[index][k]=coeff
          
          coeff=x[1][1].coeff(basis_variables[i]) * x[1][0].coeff(variables[j])
          index=i*len(variables) + j + len(basis_variables)*len(variables)
          mat[index][k]=coeff
      #print(np.transpose(mat)[k],basis_An[k])
        


    overprint(f'Created A_n image matrix at {time_conv(-compute_time+time())[0]} {time_conv(-compute_time+time())[1]}' '    \n')

    print('Nonzero before adding W:', np.count_nonzero(mat))

    try:
      with open('temp_An_matrix.pk','wb') as fi:
        pickle.dump(mat,fi)
        print(f'Saved matrix')
    except:
      pass
  
  #phi_map=mat[:,:186]
  #phi_map_op=mat[:,186:]


  k=0
  i=0
  j=0

  basis_of_W=np.zeros((len(basis_variables)*len(variables)*2,len(variables)*len(cl2)*2))
  cl2_vectors=[[x.coeff(bvar) for bvar in basis_variables] for x in cl2]

  print(cl2_vectors[0],basis_variables)
  print()
  

  for j in range(len(cl2_vectors)):
    for k in range(len(variables)):
      col=len(basis_An) + k + j*len(variables)
      for i in range(len(cl2_vectors[j])):
        coeff=cl2_vectors[j][i]
        row=i*len(variables)+k
        mat[row][col]=coeff
        basis_of_W[row][col-len(basis_An)]=coeff
        #print(index,k,basis_of_W[index][k],coeff)


        row=i*len(variables)+k
        mat[row + len(variables)*len(basis_variables)][col+len(variables)]=coeff
        basis_of_W[row + len(variables)*len(basis_variables)][col-len(basis_An)+len(variables)]=coeff
        #print(index + len(variables)*len(basis_variables),k+len(variables),basis_of_W[index + len(variables)*len(basis_variables)][k+len(variables)],coeff)

  #print('basis_of_W info: \n')
  #print(basis_of_W)
  print(basis_of_W.shape,np.count_nonzero(basis_of_W),np.linalg.matrix_rank(basis_of_W))

  print('\n')

  #print('Nonzero after adding W:', np.count_nonzero(mat))

  print('Computing Null Space...')
  
  nspace1=null_space(mat)


  nspace1=list(nspace1)

  num=0
  to_remove=[]
  for vec in nspace1:
      if ( [vec[j] for j in range(len(basis_An),len(basis_An)+len(np.transpose(basis_of_W)))] == [0 for j in  range(len(basis_An),len(basis_An)+len(np.transpose(basis_of_W)))] ):
        to_remove+=[num]
      num+=1

  removed=0
  for num in to_remove:
    nspace1=np.delete(nspace1,num-removed,0)
    removed+=1

  #print('Shape of NS1: ', nspace1, nspace1.shape)

  basis_of_image_space=[Symbol( str( (wedge,Symbol(acoord))    )   ) for wedge in basis_variables for acoord in variables ]+[Symbol( str( (Symbol(acoord),wedge))    )    for wedge in basis_variables for acoord in variables ]
  basis_of_image_space_s=[(wedge,Symbol(acoord))  for wedge in basis_variables for acoord in variables ]+[(Symbol(acoord),wedge)   for wedge in basis_variables for acoord in variables ]
  #print(basis_of_image_space)    

  truncated_nspace=np.transpose(nspace1)[:186,:]
  #print(truncated_nspace.shape)

  truncated_nspace_op=np.transpose(nspace1)[186:,:]

  #image_tnspace=np.transpose(np.matmul(phi_map,truncated_nspace))

  #image_tnspace_op=np.transpose(np.matmul(phi_map_op,truncated_nspace_op))

  #print(image_tnspace)
  #print(image_tnspace.shape)

  #print(image_tnspace_op)
  #print(image_tnspace_op.shape)

  #image_list=[]
  #image_list_op=[]

  #image_list_dict=[]

  '''

  for v in image_tnspace:
    v_expr=0
    v_dict={}
    for i in range(len(v)):
      v_expr+=v[i]*basis_of_image_space[i]
      if i<90:
        try:
          v_dict[str(basis_of_image_space_s[i][1])+',r']+=v[i]*basis_of_image_space_s[i][0]
        except:
          v_dict[str(basis_of_image_space_s[i][1])+',r']=v[i]*basis_of_image_space_s[i][0]
      if i>=90:
        try:
          v_dict[str(basis_of_image_space_s[i][0])+',l']+=v[i]*basis_of_image_space_s[i][1]
        except:
          v_dict[str(basis_of_image_space_s[i][0])+',l']=v[i]*basis_of_image_space_s[i][1]
        
    print('\n',v_dict,'\n')
    image_list+=[v_expr]
    image_list_dict+=[v_dict]

  
    
  print('\n\n\n')
  
  for v in image_tnspace_op:
    v_expr=0
    for i in range(len(v)):
      v_expr+=v[i]*basis_of_image_space[i]
      if i<90:
        try:
          v_dict[str(basis_of_image_space_s[i][1])+',r']+=v[i]*basis_of_image_space_s[i][0]
        except:
          v_dict[str(basis_of_image_space_s[i][1])+',r']=v[i]*basis_of_image_space_s[i][0]
      if i>=90:
        try:
          v_dict[str(basis_of_image_space_s[i][0])+',l']+=v[i]*basis_of_image_space_s[i][1]
        except:
          v_dict[str(basis_of_image_space_s[i][0])+',l']=v[i]*basis_of_image_space_s[i][1]
    image_list_op+=[v_expr]

  for z in zip(image_list_op,image_list):
    print(z[0]+z[1])


  is_cl2=0
  for z in zip([1,-1,-1,1,0,1,1,0,-1,0,-1,1,-1,1,-1],basis_variables):
    is_cl2+=z[0]*z[1]

  not_failed=True
  for v in image_list_dict:
    for value in list(v.values()):
      test1=value-is_cl2
      test2=value+is_cl2
      test3=2*value-is_cl2
      test4=2*value+is_cl2
      test5=value+2*is_cl2
      test6=value-2*is_cl2
      if value!=0:
        if test1==0 or test2==0 or test3==0 or test4==0 or test5==0 or test6==0:
          pass
        else:
          not_failed=False
      print('\n')

  print(not_failed)

  '''

  
  for nsv in np.transpose(truncated_nspace):
    tensorprod=0
    mathematica_string='{'
    for i in range(len(nsv)):
      tensorprod+=nsv[i]*Symbol(str(basis_An[i]))
      mathematica_string+=f'{{{basis_An[i][0]},{basis_An[i][1]},{basis_An[i][2]},{nsv[i]}}},'
    mathematica_string=mathematica_string[:-1] + '}'
    print('\n',mathematica_string,'\n')


  '''
  for v in image_list_dict:
    print(v)
    lis=[]
    for x in v.keys():
      if v[x]!=0:
        lis+=[x]
    print(lis)
    print('\n\n')
  print(np.linalg.matrix_rank(np.array([v[186:] for v in nspace1])))
'''

  #print(np.array([v[186:] for v in nspace1]))

  return nspace1.shape

       
          
 