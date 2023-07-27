from quivers import generate_quiver
from wedge import basis_of_CL2
from wedge import cl_intb_map, basis_of_CLn,expected_dimension_of_cln_gr2m
import pickle

def menu(info,results,cl2_dictionary):
    
  do_find_x_coords=False
  do_find_a_coords=False
  do_find_clusters=False
  do_find_x1x=False

  skip_input_1=False
  skip_input_2=False
  skip_input_3=False


  option1=''
  test_mode=False
  skip_numeric_test_error=False
  test_dupes=False

  test_continually=''


  m=4

  muted=False

  option1=''
  while option1 not in ['gr','qu']:
    option1=input('(Main menu) Enter your choice: ')
    if option1=='t':
      print('Resulintg NS dim: ', basis_of_CLn(results['cl,gr24,13'][0],results['cl,gr24,13'][1],cl2_dictionary['x1x,gr24,7'])[1])
      print('Expected: ', expected_dimension_of_cln_gr2m(3,4))
    elif option1=='cl':
      print(cl2_dictionary[input('Input key to retrieve basis for cl2: ')])
    elif option1=='skip':
      print('Ignoring numeric test error.')
      skip_numeric_test_error=True
    elif option1=='skiptest':
      print('Enabling numeric test for duplicates.')
      test_dupes=True
    elif option1=='info':
      print(info)
    elif option1=='debug':
      test_mode=True
    elif option1=='co':
      option1='gr'
    elif option1=='keys':
      print(list(results.keys()))
      print()
    elif option1.startswith('cln'):
      try:
        test_continually=option1[3]
      except:
        pass
      option1='gr'
      skip_input_1=True
      do_find_x1x=True
      find_coords='x1x'
      if test_continually=='c':
        muted=True
        skip_input_2=True
        skip_input_3=True
        num_of_mutations=50
        if 4*m < 100:
          iteration_limit=4*m
        else:
          iteration_limit=0
    elif option1 in results.keys():
      display_result=list(results[option1][0])
      print(f'\nPrinting final result: "{option1}"\n')
      print(display_result,results[option1][1:])
      print()
      new_option=input(f' integer in range (1,{len(display_result)}): return that element. \n list: return each element printed to its own line. \n cln: computes the basis and dimension of CLn (must have chosen to compute X coordinates with 1+X). \n vars: returns the list of variables used.' + '\n'+'Enter your choice: ')
      while new_option!='':
      
        if new_option in [str(x+1) for x in range(len(display_result))]:
          print(f'Result #{new_option}: \n {display_result[int(new_option)-1]}')


        elif new_option=='list':
          for x in display_result:
            print(x)


            
        elif new_option.startswith('cl') and option1.startswith('x1x'):
          if new_option[2]=='2':
            print('Computing CL2...')
            basis_and_dim=basis_of_CL2(display_result,results[option1][1])
            print('\nBasis:' )
            print(basis_and_dim[1])
            print('Dimension: ', basis_and_dim[2])
            cl2_dictionary[option1]=basis_and_dim[1]
            if input('Enter "s" to save CL2 basis: ')=='s':
                with open('cl2_dictionary.pk','wb') as fi1:
                    pickle.dump(cl2_dictionary,fi1)

          elif int(new_option[2:])>2:
            cluster_adjacency=[tuple(x) for x in results[input('Enter the key for the clusters: ')][0]]
            print('Not implimented')

        elif new_option=='vars':
          print('All variables:')
          print(results[option1][1])
          print('New variables in terms of plucker coordinates:')
          print(results[option1][2])
          

        new_option=input('Choose again, or enter nothing to return to the main menu: ')
    else:
      print('Invalid option, try again. \n')


  if not skip_input_3:
    quiver_data_s=input('Enter the quiver to use.\n'+
   'gr3m uses the quiver for gr(3,m) with plucker coordinate vertices. (m=6,7,8)'+
    '\ngr2n uses the quiver for gr(2,n) with plucker coordinate vertices. (n=4,5,6,7,8,9)\n'+
  'ex2 uses the quiver 1->2 with default vertex labeling.\n'+
  'ex3 uses the quiver 1->2->3 with default vertex labeling.'+
  '\nEnter:')  
    quiver_data=generate_quiver(quiver_data_s)
  
  if option1=='gr' and not skip_input_1:
    find_coords=input('Enter "x" to find X-coords, "x1x" to find X-coords with 1+X-coords, "a" to find A-coords, or "cl" to find all the clusters: ')
    if find_coords=="x":
      do_find_x_coords=True
      do_find_a_coords=False
      do_find_x1x=False
      do_find_clusters=False
    elif find_coords=='a':
      do_find_a_coords=True
      do_find_x_coords=False
      do_find_x1x=False
      do_find_clusters=False
    elif find_coords=="cl":
      do_find_clusters=True
      do_find_x_coords=False
      do_find_a_coords=False
      do_find_x1x=False
    elif find_coords=="x1x":
      do_find_x1x=True
      do_find_x_coords=False
      do_find_a_coords=False
      do_find_clusters=False
      
  if not skip_input_2:
    num_of_mutations=int(input('Enter the number of random mutations (not recommended to be larger than 100): '))
    iteration_limit=input('Enter the iteration limit or leave blank or enter 0 to iterate forever: ')
    
  return do_find_x_coords,do_find_a_coords,do_find_x1x,do_find_clusters,test_mode,test_dupes,skip_numeric_test_error,skip_input_1,skip_input_2,skip_input_3,muted,m,option1,num_of_mutations,quiver_data,quiver_data_s,find_coords,iteration_limit
