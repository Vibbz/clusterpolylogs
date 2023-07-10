import numpy as np
from random import random

def grassmanian_vertex_indices(i,j,p,q,n):
  index=[]
  if i<=j:
    index+=list(range(i+1,p+1))
    index+=list(range(n+1-j,n+i-j+1))
  if i>j:
    index+=list(range(1,i-j+1))
    index+=list(range(i+1,p+1))
    index+=list(range(n+1-j,n+1))
  return index

def generate_quiver(user_input):
  while user_input not in ['ex2','ex3','gr36'] and user_input.startswith('gr2')==False:
    user_input=input('Please enter an available quiver: ')
  quiver=[[0]]
  verts=0
  mutables=0
  verts_num=[[0]]
  
  if user_input=='ex2':
    quiver=np.array([[0,1],[-1,0]]) #This is the basic two-vertex quiver with one edge going from one vertex to the other: 1->2
    mutables=1
    verts=0
    
  if user_input=='ex3':
    quiver=np.array([[0,1,0],[-1,0,1],[0,-1,0]] ) #This is the basic three-vertex quiver with two edges, 1 -> 2 -> 3
    mutables=2
    verts=0
    verts_num=[random(),random()*100,random()*10000]

  if user_input=='gr36':
    quiver=np.array([[0,1,1,-1,-1,0,0,0,0,0],[-1,0,0,1,0,1,-1,0,0,0],[-1,0,0,1,0,0,0,1,-1,0],[1,-1,-1,0,0,0,1,0,1,-1],[1,0,0,0,0,0,0,0,0,0],[0,-1,0,0,0,0,0,0,0,0],[0,1,0,-1,0,0,0,0,0,0],[0,0,-1,0,0,0,0,0,0,0],[0,0,1,-1,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0,0] ])
    verts=['a236','a235','a136','a356','a123','a234','a345','a126','a156','a456']
    mutables=3
    
    dummy_matrix=np.array([[n+6*i+random() for n in range(6)] for i in range(3)])

    verts_num=[np.linalg.det(np.array([ [row[int(i)-1] for i in list(plucker[1:])] for row in dummy_matrix ] )) for plucker in verts]

  if user_input.startswith('gr2'):
    print('\nGenerating quiver for Gr(2,n)')
    try:
      n=int(user_input[3])
    except:
      n=int(input('Enter dimension of space (n): '))
    p=2
    q=n-p
  
    print('Computing the quiver for Gr({0},{1})...'.format(p,n))
    
    quiver_size=p*q+1
    mutables=q-2
    
    quiver=np.zeros((quiver_size,quiver_size)).astype(int)
  
    for i in range(mutables+1):
      quiver[i][i+1], quiver[i][i+q], quiver[i][i+q+1] = 1, 1, -1
  
    quiver[len(quiver)-1][0]=1
  
    negtransp=(-1)*np.transpose(quiver)
    for i in range(len(negtransp)):
      for j in range(len(negtransp)):
        if negtransp[i][j]!=0:
          quiver[i][j]=negtransp[i][j]


    #The order of vertices being added is very important!
    verts=[]
    for i in range(1,p+1):
      for j in range(1,q+1):
        index=grassmanian_vertex_indices(i,j,p,q,n)
        verts+=['a{0}{1}'.format(index[0],index[1])]
    verts+=['a12']

    dummy_matrix=np.array([[n+6*i+random() for n in range(6)] for i in range(2)])

    verts_num=[np.linalg.det(np.array([ [row[int(i)-1] for i in list(plucker[1:])] for row in dummy_matrix ] )) for plucker in verts]
    
  if user_input=='gras':
    print('Not yet implemented.')
    pass
    n=int(input('Dimension of the space: '))
    p=int(input('Dimension of subspaces: '))
    q=n-p

    print('Computing the quiver for Gr({0},{1})...'.format(p,q))

    quiver_size=p*q+1
    mutables=(p-1)*(q-1)-1

    quiver=np.zeros((quiver_size,quiver_size))

    for i in range(len(quiver)):
      if i<=mutables:
        if i<=q-3:
          quiver[i][i+1]=1
          quiver[i][i+q-1]=1
        elif i%(q-2)==0:
          quiver[i][mutables+1+i//(q-2)]=1
          quiver[i][i+q-1]=1
        elif i%(q-2)==1 and i<=(p-2)*(q-1)-1:
          quiver[i][i+q-1]=1
          quiver[i][i+1]=1
        elif i%(q-1)==1 and i>(p-2)*(q-1):
          quiver[i][i+1]=1
          quiver[i][mutables+1+p]
        elif i<(mutables-q+3):
          quiver[i][i+1]=1
          quiver[i][i+q-1]=1
          quiver[i][i-q]=1
        else:
          quiver[i][i-q]=1
          quiver[i][i+mutables+2]
      elif i==mutables+1:
        quiver[i][0]=1
      elif mutables+2<i<=mutables+p:
        quiver[i][(i-mutables-2)*(q-1)]
      elif i>mutables+p+1:
        quiver[i][i-p-q]=1

  if user_input in ['ex2','ex3']:
    verts=['x{0}'.format(i+1) for i in range(len(quiver))]
  
  return quiver,verts,mutables,verts_num
