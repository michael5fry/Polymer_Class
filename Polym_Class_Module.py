import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

class polymer:
    def __init__(self, polymer_name):
        self.polymer_name = polymer_name
        
        
    intro1 = """Part 1: EXTRACTION From Gaussian Output the Bands
       THE ONLY FUNCTION TO CALL IS:
       
       Outpout_to_DataFrame_coverted_shifted(file_name_gaussian, K_points_requested, DELTA_EN, conv)
       
       Needed variables:
       file_name_gaussian ---> file path and name of the Gaussian Output 
       K_points_requested ---> Number of K point needed
       DELTA_EN ---> Shift in energy put to 0 if you don't want it
       conv ---> Conversion of the initial data to the Unit of Measure of your preference
                   It is Multiplied for the Number of interest
    """
    
    #PART 1
    #-----------------
    # GENERAL PART for ALL Gaussian Bands 
    # PART TO TAKE THE BANDS AND PUT THEM IN A SEPATE CONTENT
    # NOW IF YOU WANT YOU CAN WRITE THE FILE
    #------------------
    def finder_mathod_A(self, list_from_fp, string_to_search):
        for line_num in range(len(list_from_fp)):
            data = list_from_fp[line_num].find(string_to_search)
            if data !=-1:
               #print(line_num)
               return line_num
           
    ## Analysis of the Simple one:
    def wrangling_bands_stage0(self):
        with open(self.polymer_name) as fp:   
            start = [line for line in fp]
            
        init_line_index = self.finder_mathod_A(start, ' energies kxyz=     0     0')
        # End of lines of Interest:
        end_line_index = self.finder_mathod_A(start, ' HOCO           at k-point   ')
        # content wanted
        content_wanted = start[init_line_index : end_line_index ]
        
        return content_wanted
    
    
    #----------------------------------------
    # SPECIFIC PART FOR Gaussian Bands WITH MORE THAN 10 BANDS
    #----------------------------------------
    
#Program to divide how MANY lines there are for the individual k point
    
    def division_k1_and_k2(self):
        data = self.wrangling_bands_stage0() ################################################################################
        # take the 1 point that gives you the division
        iterator = -1
        list_of_division_indexes = []
        for line_i in range(len(data)):
            main_point_line = int(data[line_i].split()[5])
           
           # here created an iterator to get the number of lines for specific k point
            if main_point_line > iterator:
                iterator = main_point_line
                list_of_division_indexes.append(line_i)
        
        return list_of_division_indexes
    
    
#--------------------------------
#     VERSION 3
#     ----------------------------

#-------------------------------------
#     WORKING RIGHT NOW
#     ---------------------------------------------
    def gaussian_band_for_1_k_point_v4(self, band_i, k_point_i ):
        # SCRIPT VALID FOR EVERY 48
        energies_n = [k_point_i]
        for line in band_i:
            # here the part of the line in interest
            start_line_interest = line.split('     ')[-1].strip(' eV\n') 
            length_int = len(start_line_interest.split())
            
            # conditions to find the various energies
            if '**' in start_line_interest:
                string_star = start_line_interest.count('*') * '*'
                # first splitting
                string_energies = start_line_interest.split(string_star)[-1]
                # list energies by splitting with the minus
                list_energies = string_energies.split('-')
                #
                list_energies1 = filter(lambda item: item!='',list_energies)
                # append each to energies_n
                for i in list_energies1:
                    #print('Round1:',float(i )*-1)
                    #print(line)
                    energies_n.append(float(i )*-1)
                   
            elif '-' in start_line_interest and start_line_interest.count('-')==5:
                
                list_energies = start_line_interest.split('-')
                
                list_energies1 = filter(lambda item: item!='1' and item!='1 ',list_energies)
                list_energies2 = filter(lambda item: item!='' and item!='1  ',list_energies1)
                list_energies3 = filter(lambda item: item!='0',list_energies2)
                #list_energies4 = filter(lambda item: item!='1 ',list_energies3)
                # append each to energies_n
                for i in list_energies3:
                    #print('Round2:',float( i)*-1)
                    #print(line)
                    energies_n.append(float( i)*-1)# bro
                    #print('Round2:',float(i )*-1)
                    #print(list_energies)
           
            else:
                list_energies = start_line_interest.split()
                list_energies1 = filter(lambda item: item!='1',list_energies)
                list_energies2 = filter(lambda item: item!='0',list_energies1)
               
                for i in list_energies2:
                    energies_n.append(float( i))
                   
        return np.array(energies_n)
    
    
#STRUCTURED WRANGLED DATA:
    def getting_energies_119_gauss(self, K_points_requested):
        # call the function to Wrangle the initial data
        wrangled_data = self.wrangling_bands_stage0() ######################################################################
        number_of_line_for_each_k = self.division_k1_and_k2()[1]
        
        energies_K_point_requested=[]
        for k_point_i in range(1,K_points_requested +1):
            
            begin_set_i = number_of_line_for_each_k * k_point_i
            end_set_i = number_of_line_for_each_k * (1+k_point_i)
            
            structured_wrangled_data_N = wrangled_data[begin_set_i:end_set_i]
           
            # Call the effective function that UnWrangle the Energies Data
            energ_k_i = self.gaussian_band_for_1_k_point_v4( structured_wrangled_data_N, k_point_i ) # HERE CHANFE
            energies_K_point_requested.append(energ_k_i)  
           
        return energies_K_point_requested#[0]
    
    
    
    
    def Outpout_to_DataFrame_coverted_shifted(self, K_points_requested, DELTA_EN, conv): # in eV
        
        all_energies_K_req = self.getting_energies_119_gauss(K_points_requested) ###########################
        # Creation of a DataFrame to read the Gaussian File
       
        df_band_polym_A = pd.DataFrame(all_energies_K_req)#, columns=['K_points'])
        df_band_polym_A  = df_band_polym_A.rename({0:'K_points'}, axis=1) # ---> RENAME JUST 0 as K_points
       
        df_band_polym_A['K_points'] = list(reversed(df_band_polym_A['K_points']))
        for column in range(1, len(df_band_polym_A.columns)-1):
            df_band_polym_A[column]= list(reversed((df_band_polym_A[column] * conv) + DELTA_EN)) # New Origin 
            # I 
    
        return df_band_polym_A
    





"""Part 2: CALCULATION OF THE COUPLING BETWEEN THE 2 SITES
    
    
    Outpout_to_DataFrame_coverted_shifted(file_name_gaussian, K_points_requested, DELTA_EN, conv)
    
    Needed variables:
     --->  
     ---> 
     ---> 
     ---> 
def DF_for_coupling_J( *args): # Input are the Polymers I want with all they need
    
    # All entries
    polym_list = []
    for polym_i in args:
        # Polymer A
        df_band_polym_A =    polym_i[0]
        trigonometric_func = polym_i[1]
        band_i_up =          polym_i[2] 
        band_j_down =        polym_i[3]
        entry_i = extraction_parameters_get_table(df_band_polym_A, trigonometric_func, band_i_up, band_j_down)
        
        polym_list.append(entry_i)
    """ 
def BZ_eigenvalue(df_band_polym_A, band_i): # NO MORE NEED: They need to be taken on the Opposite because the 
                                            # The K point 1 is in BZ and the 50 is at Gamma
    return float(df_band_polym_A[band_i][len(df_band_polym_A) - 1])

def GM_eigenvalue(df_band_polym_A, band_i):
    return float(df_band_polym_A[band_i][0])




def extraction_parameters_get_table(df_band_polym_A, trigonometric_func, band_j_down, band_i_up):
    
    if trigonometric_func == 'cos':
        # eigenvalues found at K = pi/2 ---> First
        val_eps_down = BZ_eigenvalue(df_band_polym_A, band_j_down)
        val_eps_up = BZ_eigenvalue(df_band_polym_A, band_i_up)
        print(val_eps_down)
        delta_eps = val_eps_down - val_eps_up
        # eigenvalues found at K = 0 ---> Second 
        delta_lambda = GM_eigenvalue(df_band_polym_A, band_j_down) - GM_eigenvalue(df_band_polym_A, band_i_up)
        
    elif trigonometric_func == 'sin':
        # eigenvalues found at K = 0 ---> First
        val_eps_down = GM_eigenvalue(df_band_polym_A, band_j_down)
        val_eps_up = GM_eigenvalue(df_band_polym_A, band_i_up)
        
        delta_eps = val_eps_down - val_eps_up
        # eigenvalues found at K = pi/2 ---> Second 
        delta_lambda = BZ_eigenvalue(df_band_polym_A, band_j_down) - BZ_eigenvalue(df_band_polym_A, band_i_up)
        
    # Value of the coupling:
    J = (1/4) * np.sqrt(delta_lambda**2 - delta_eps**2)
    
    return J, val_eps_down, val_eps_up
        

def DataFrame_for_coupling_J( *args): # Input are the Polymers I want with all they need
    band_dict = {'DPPBTz' : [131,132,133,134], 'IDTBT':[229,230,231,232], 'P3HT':[47,48,49,50], 'PBTTT':[80,81,82,83]}

    # All entries
    polym_list = []
    for polym_i in args:
        J_i, val_eps_down_i, val_eps_up_i = extraction_parameters_get_table(polym_i[0], polym_i[1], polym_i[2], polym_i[3])

        for i in ['DPPBTz', 'IDTBT', 'P3HT','PBTTT']:
            if polym_i[2] == band_dict[i][0] or polym_i[2] == band_dict[i][1] :
                name_polym = i + ' VB'
                
            elif polym_i[2] == band_dict[i][2] or polym_i[2] == band_dict[i][3] :
                name_polym = i + ' CB'
                
        entry_i = [name_polym, J_i, val_eps_down_i, val_eps_up_i] # Row_i
        polym_list.append(entry_i) # I have a List of Rows
             
    df_coupling_polym = pd.DataFrame(polym_list)
    df_coupling_polym.columns = ['Polym_name', 'J [eV]', 'eps_A [eV]', 'eps_D [eV]']
    return df_coupling_polym
 
    



#-----------------------------------
# TRUE PLOTTING PART   
#-----------------------------------
    
class plotting_env():

    
    def plotting_faster(self, df_polym_A, band_i, color_linA, color_linB, label_name, more_bands ):
            
        k_point_pi= np.arange(0, np.pi+0.15, 0.0325*2)
        #k_point_rev = list(reversed(list((k_point_pi)) ) ) # Here the List is Reversed
        plt.plot(k_point_pi, df_polym_A[band_i] , color= color_linA, linewidth=2 , label=label_name )
        #if more_bands == 'yes +1':
        plt.plot(k_point_pi, df_polym_A[band_i + 1], color= color_linB, linewidth=2  )
        #elif more_bands == 'yes -1':
        #§plt.plot(k_point_pi, df_polym_A[band_i - 1], color= color_linB, linewidth=2  )
            
    def plotting_faster1(self, df_polym_A, band_i, color_linA, color_linB, label_name, more_bands ):
            
        k_point_pi= np.arange(0, np.pi+0.15, 0.0325*2)
        #k_point_rev = list(reversed(list((k_point_pi)) ) ) # Here the List is Reversed
        plt.plot(k_point_pi, df_polym_A[band_i] , color= color_linA, linewidth=2 , label=label_name )
                
                
                
    def analitical_band(self,epsilon_A, epsilon_D, J_coupling, type_TB_bands):
        k_point_pi= np.arange(0, np.pi+0.1, 0.0325)
        a=1 # Lattice Parameter
    
        sum_A_D= epsilon_A + epsilon_D
        delta_A_D = epsilon_A - epsilon_D
        #Branch 1 --> Acceptor
        if type_TB_bands == 'sin':
            lambda1 = 0.5*(sum_A_D - np.sqrt(delta_A_D**2 + (4 * J_coupling * np.sin(k_point_pi*a*0.5))**2 ) ) 
            lambda2 = 0.5*(sum_A_D + np.sqrt(delta_A_D**2 + (4 * J_coupling * np.sin(k_point_pi*a*0.5))**2 ) ) 
        
        elif type_TB_bands == 'cos':
            lambda1 = 0.5*(sum_A_D - np.sqrt(delta_A_D**2 + (4 * J_coupling * np.cos(k_point_pi*a*0.5))**2 ) ) 
            lambda2 = 0.5*(sum_A_D + np.sqrt(delta_A_D**2 + (4 * J_coupling * np.cos(k_point_pi*a*0.5))**2 ) )
        
        return lambda1, lambda2 #- for HOMO-1 and + for HOMO

            
    def TB_plotting(self, epsilon_A, epsilon_D, J_coupling, type_TB_bands, color_linA, label_name):
        # taking eigenvalues
        lambda_1, lambda_2 = self.analitical_band(epsilon_A, epsilon_D, J_coupling, type_TB_bands)
        k_point_pi= np.arange(0, np.pi+0.1, 0.0325)
    
        # plotting
        plt.plot(k_point_pi, lambda_1, color= color_linA , ls='-.', label=label_name )
        plt.plot(k_point_pi, lambda_2,  color= color_linA , ls='-.')


def plotting_to_MODIFY(title, polym_1, trig_func, VB_band, CB_band):
    plot1 = plotting_env()
    plt.figure(figsize=(7.5,6.5))
    plt.title(title)
    
    value_to_see = 20
    plt.xticks([0,np.pi],['GM', 'BZ'], fontsize= value_to_see-7)
    plt.yticks( fontsize= value_to_see-7)
    plt.xlim(0,np.pi)
    plt.xlabel('K points', fontsize= value_to_see-3)
    plt.ylabel('Energy eV', fontsize=value_to_see-3)

    VB = DataFrame_for_coupling_J([polym_1, trig_func, VB_band-1, VB_band]).iloc[0,3]
    #plt.plot(3, VB, '*')
    VB_1 = DataFrame_for_coupling_J([polym_1, trig_func, VB_band-1, VB_band]).iloc[0,2]
    J_VB = DataFrame_for_coupling_J([polym_1, trig_func, VB_band-1, VB_band]).iloc[0,1]

    plot1.plotting_faster(polym_1, VB_band-1, color_linA='b', color_linB='b', label_name='VB-1 VB', more_bands='yes -1' )
    plot1.TB_plotting(VB_1, VB, J_VB, type_TB_bands=trig_func, color_linA = 'r', label_name = 'VB-1 VB')

    CB = DataFrame_for_coupling_J([polym_1, trig_func, CB_band, CB_band+1]).iloc[0,3]
    #plt.plot(0.05, CB, '*')
    CB_1 = DataFrame_for_coupling_J([polym_1, trig_func, CB_band, CB_band+1]).iloc[0,2]
    J_CB = DataFrame_for_coupling_J([polym_1, trig_func, CB_band, CB_band+1]).iloc[0,1]

    plot1.plotting_faster(polym_1, CB_band, color_linA='b', color_linB='b', label_name='CB CB+1', more_bands='yes +1' )
    plot1.TB_plotting(CB, CB_1, J_CB, type_TB_bands=trig_func, color_linA = 'r', label_name = 'CB CB+1')

    plt.legend(fontsize= 7, ncol=2)
    return plt.show()
    
    
    
def MO_coefficient_calculator(file_name, num_MO, full_basis_set,  basis_set_fragment1):
    # open file and save it as single 1D array 
    with open(file_name) as fp:
        start = [round(float(element_i),12) for line in fp for element_i in line.split()]
    start = np.array(start)
    
    # Transform the file in a 2D matrix with full_basis_set as the Number of Rows and Columns
    # Transpose it, because THE initial full_basis_set rows HAS BECOME A Column
    
    matrix_start = start.reshape(full_basis_set, full_basis_set).transpose() # make the single array in a 526*526 matrix
    # NB as a Matrix first Row and then Column
    
    # Take the orbital of interest num_MO
    # this is the the way to Take the column that correspond to the Molecular Orbital num_MO 
    orbital_of_interest = matrix_start[:,num_MO -1] 
    
    # Find the Coefficient requested
    coeff_1_array = orbital_of_interest[0 : basis_set_fragment1]
    coeff_2_array = orbital_of_interest[basis_set_fragment1 : :]

    coeff_1 = np.sum(coeff_1_array**2)
    coeff_2 = np.sum(coeff_2_array**2)

    print(f' The coefficient for the Fragment1 is: {coeff_1}\n',
          f'The coefficient for the Fragment2 is: {coeff_2}\n',
          f'We can check that the summation bring to 1: {round(coeff_1 + coeff_2, 5)}')
    return round(coeff_1, 2), round(coeff_2, 2)