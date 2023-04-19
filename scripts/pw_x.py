# %% [markdown]
# # DFT SCF Step (PW.X)

##
# @file pw_x.py
# 
# @brief PW_X step input widget creation
# 
# @page PW_X
# @section description_pw Description
# @ref pw_x.py is a jupyter notebook that contains the input widgets for the PW_X step
# 
# @section todo_pw TODO
# 
# @section notes_pw Notes
# @ref Quantum Espresso Input_PW documentation: https://www.quantum-espresso.org/Doc/INPUT_PW.html
#
# @section libraries_pw Libraries/Modules
# - ipywidgets  : https://ipywidgets.readthedocs.io/en/latest/
#  - Tab        : https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20List.html#Tab
#  - Box        : https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20List.html#Box
#  - VBox       : https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20List.html#VBox
#  - HBox       : https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20List.html#HBox
#  - GridBox    : https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20List.html#GridBox
#  - Layout     : https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20List.html#Layout
#  - HTML        : https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20List.html#HTML
#  - Label      : https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20List.html#Label
#  - IntText    : https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20List.html#IntText
#  - FloatText  : https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20List.html#FloatText
#  - Dropdown   : https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20List.html#Dropdown
#  - Text       : https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20List.html#Text
#  - Textarea   : https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20List.html#Textarea
# - os          : https://docs.python.org/3/library/os.html
# - nanohub     : https://nanohub.org/tools/quantumespresso/
# - styles      : styles.py
#  - Runs the styles.ipynb notebook 
# 


# %%
### Import libraries

import numpy as np
from ipywidgets import Tab, Box, VBox,HBox, GridBox, Layout
from ipywidgets import HTML, Label, IntText, FloatText, Dropdown, Text, Textarea, Button

import hublib.use
import os
%run styles.ipynb

# %% [markdown]
# ## Documentation Link

# %%
os

# %% [markdown]
# Check out the documentation on inputs variables needed for this step with this link:
# <a href="https://www.quantum-espresso.org/Doc/INPUT_PW.html" target="_blank">INPUT_PW.html</a>

# %%
documentation_link = HTML(value='<a href="https://www.quantum-espresso.org/Doc/INPUT_PW.html" target="_blank">INPUT_PW.html</a>')

# %% [markdown]
# ## Control Tab

# %%
#### 
# Create Control form
material_prefix = Text(name="Material Prefix", value="si")

restart_mode = Dropdown(name="Restart Mode", value="from_scratch", options=["from_scratch", "restart"])

pseudo_dir = Text(name="Pseudopotential File Directory", value=os.getcwd())

outdir = Text(name="Current Directory", value=os.getcwd())

max_runtime = FloatText(name='Max Run Time (Seconds)', value=10000000)


####
# Creates the info buttons and adds in the description when you hover your cursor over the button

icon_material_prefix = Button( icon='fa-info-circle',
                                tooltip='prepended to input/output filenames: \n prefix.wfc, prefix.rho, etc.', 
                                layout = Layout(width='40px'))

icon_restart_mode = Button( icon='fa-info-circle',
                            tooltip='"from_scratch" : ' '\n' ' From scratch. This is the normal way to perform a PWscf calculation' '\n' '\n'
                               '"restart : "' '\n' 'From previous interrupted run. Use this switch only if you want to' '\n' 
                               'continue, using the same number of processors and parallelization,' '\n'
                               'an interrupted calculation. Do not use to start a new one, or to' '\n'
                               'perform a non-scf calculations.  Works only if the calculation was' '\n'
                               'cleanly stopped using variable max_seconds, or by user request' '\n'
                               'with an "exit file" (i.e.: create a file "prefix".EXIT, in directory' '\n'
                               '"outdir"; see variables prefix, outdir). The default for' '\n'
                               'startingwfc and startingpot is set to "file".'
                               , 
                             layout = Layout(width='40px'))


icon_pseudo_dir = Button( icon='fa-info-circle',
                                tooltip='directory containing pseudopotential files', 
                                layout = Layout(width='40px'))

icon_outdir = Button( icon='fa-info-circle',
                                tooltip='input, temporary, output files are found in this directory' '\n' '\n'
                                         'DEFAULT: current directory (./)', 
                                layout = Layout(width='40px', grid_row_gap='30px', justify_content='flex-start'))


icon_max_runtime = Button( icon='fa-info-circle',
                                tooltip='Jobs stops after max_seconds CPU time. Use this option in conjunction' '\n'
                                          'with option restart_mode if you need to split a job too long to' '\n'
                                          'complete into shorter jobs that fit into your batch queues.' '\n'
                                          'DEFAULT: 1E7, or 150 days, i.e. no time limit', 
                                layout = Layout(width='40px'))

#####


form_items = [
    # I AM ADDING NONBREAKING SPACES TO TRY AND FIX THE ALIGNMENT#
    # IF THERE IS A UNICODE TAB THAT WORKS HERE PLEASE SWAP IT OUT. IT IS CURRENTLY 1:20AM AND I AM EYEBALLING THE AMOUNT OF SPACES
    
    HBox([Label(value='Material Prefix'), 
          Box([material_prefix, icon_material_prefix]),
         ], layout=form_item_layout()),
         
    HBox([Label(value='Restart Mode'), 
         Box([restart_mode, icon_restart_mode]),
        ], layout=form_item_layout()),
    
    HBox([Label(value='Pseudopotential File Directory'), 
          Box([pseudo_dir, icon_pseudo_dir])
         ], layout=form_item_layout()),
    
    HBox([Label(value='Current Directory'), 
          Box([outdir, icon_outdir])
         ], layout=form_item_layout()),
    
    HBox([Label(value='Max Run Time (Seconds)'), 
          Box([max_runtime, icon_max_runtime])
         ], layout=form_item_layout()),
    
    HBox([Label(value='Documentation URL:'), documentation_link], layout=Layout())
    
]

controls_box = Box(form_items, layout=box_layout(60))

                   
controls_box

# %% [markdown]
# ## System Tab

# %%
# Create System form

ibrav = Dropdown(name="Bravals-lattice index", 
                 value=0, 
                 options={"0":0,"1":1,"2":2,"3":3,"-3":-3,"4":4,"5":5,"-5":-5,"6":6,"7":7,"8":8,"9":9,"10":10,"11":11,"12":12,"-12":-12,"13":13,"-13":-13,"14":14},
)
celldm_1 = FloatText(name="Crystallographic Constant",value=10.2094)

nat = FloatText(name="Number of atoms in the cell system" ,value=2)

ntyp = FloatText(name="Number of Types of Atoms",value=1)

ecutwfc = FloatText(name="Wavefunction Kinetic NRG Cutoff (Ry)",min=25.0,max=100.0,value=25.0)


####
#Creates the info buttons and adds in the description when you hover your cursor over the button

icon_ibrav = Button( icon='fa-info-circle',
                                tooltip='crystal structure type, the number maps to a specific structrure that can be looked up in the official documentation''\n''SEE LINK BELOW', 
                                layout = Layout(width='40px'))

icon_celldm_1 = Button( icon='fa-info-circle',
                                tooltip='SEE DOCUMENTATION LINK BELOW', 
                                layout = Layout(width='40px'))

icon_nat = Button( icon='fa-info-circle',
                                tooltip='number of atoms in the unit cell (ALL atoms, except if space_group is set, in which case, INEQUIVALENT atoms)', 
                                layout = Layout(width='40px'))

icon_ntyp = Button( icon='fa-info-circle',
                                tooltip='number of TYPES of atoms in the unit cell', 
                                layout = Layout(width='40px'))

icon_ecutwfc = Button( icon='fa-info-circle',
                                tooltip='kinetic energy cutoff (Ry) for wavefunctions', 
                                layout = Layout(width='40px'))


####



#add icon button to array[Label(value='Bravals-lattice Index'), ibrav]
form_items = [
    HBox([Label(value='Bravals-lattice Index'), 
          Box([ibrav, icon_ibrav])
         ], layout=form_item_layout()),
    HBox([Label(value='Crystallographic Constant'), 
          Box([celldm_1, icon_celldm_1])
         ], layout=form_item_layout()),
    HBox([Label(value='Number of atoms in the cell system'), 
          Box([nat, icon_nat])
         ], layout=form_item_layout()),  
    HBox([Label(value='Number of Types of Atoms'), 
          Box([ntyp, icon_ntyp])
         ], layout=form_item_layout()),
    HBox([Label(value='Wavefunction Kinetic NRG Cutoff (Ry)'), 
          Box([ecutwfc, icon_ecutwfc])
         ], layout=form_item_layout()),
    
    HBox([Label(value='Documentation URL:'), documentation_link], layout=Layout())
]

system_box = Box(form_items, layout=box_layout(65))

#system_box

# %% [markdown]
# ## Electrons

# %%
# Electrons

diagonalization = Dropdown(name="Diagonalization",value="david",
                           options=["david", "cg", "ppcg", "paro", "rmm-davison"]
                          )                        

mixing_beta = FloatText(name="Mixing Beta",value=0.70)


con_threshold = FloatText(name="Convergence Threshold",value=1.0e-13)

atom = Text(name="Atom",value="Si")

atomic_mass = FloatText(name="Atomic Mass",value=28.0855)

pseudo_filename = Text(name="Pseudopotential File Name",
                       value='data/pseudopotentials/sg15_oncv_upf_2jun20/Si_ONCV_PBE-1.2.upf',                 
                      )   
atomic_coord_type = Dropdown(name="Atomic Coordinate Type",
                             options= ['alat','crystal','cartesian'],value= 'alat',
                            )

####
#Creates the info buttons and adds in the description when you hover your cursor over the button

icon_diagonalization = Button( icon='fa-info-circle',
                                tooltip='Available options are:''\n'
                              
                                          '"david" :' '\n'
                                          'Davidson iterative diagonalization with overlap matrix (default). Fast, may in some rare cases fail.''\n' '\n'
                              
                                          '"cg" :' '\n'
                                          'xa0Conjugate-gradient-like band-by-band diagonalization. MUCH slower than "david" but uses less memory and is (a little bit) more robust.' '\n''\n'
                              
                                          '"ppcg" :' '\n'
                                          'xa0PPCG iterative diagonalization' '\n''\n'
                                            
                                          '"paro", "Par0" :' '\n'
                                          'ParO iterative diagonalization''\n''\n'
                              
                                           '"rmm-davidson", "rmm-paro" :' '\n'
                                          'RMM-DIIS iterative diagonalization. To stabilize the SCF loop RMM-DIIS is alternated with calls to Davidson or ParO solvers depending on the string used.''\n'
                                          'Other variables that can be used to tune the behavior of RMM-DIIS are:  diago_rmm_ndim and diago_rmm_conv', 
                                layout = Layout(width='40px'))


icon_mixing_beta = Button( icon='fa-info-circle',
                                tooltip='mixing factor for self-consistency', 
                                layout = Layout(width='40px'))


icon_con_threshold = Button( icon='fa-info-circle',
                                tooltip='Convergence threshold for selfconsistency:''\n'
                                       'estimated energy error < conv_thr' '\n'
                                        '(note that conv_thr is extensive, like the total energy).''\n'

                                        'For non-self-consistent calculations, conv_thr is used to set the default value of the threshold (ethr) for iterative diagonalization: see diago_thr_init', 
                                layout = Layout(width='40px'))


icon_atom = Button( icon='fa-info-circle',
                                tooltip='label of the atom. Acceptable syntax:chemical symbol X (1 or 2 characters, case-insensitive) or chemical symbol plus a number or a letter, as in "Xn" (e.g. Fe1) or "X_*" or "X-*"' '\n'
                                       '(e.g. C1, C_h; max total length cannot exceed 3 characters)', 
                                layout = Layout(width='40px'))



icon_atomic_mass = Button( icon='fa-info-circle',
                                tooltip='mass of the atomic species [amu: mass of C = 12]''\n'
                                        'Used only when performing Molecular Dynamics run or structural optimization runs using Damped MD.''\n'
                                        'Not actually used in all other cases (but stored in data files, so phonon calculations will use these values unless other values are provided)', 
                                layout = Layout(width='40px'))


icon_pseudo_filename = Button( icon='fa-info-circle',
                                tooltip='File containing PP for this species.' '\n'

                                        'The pseudopotential file is assumed to be in the new UPF format. If it doesnt work, the pseudopotential format is determined by the file name:' '\n'

                                        '*.vdb or *.van Vanderbilt US pseudopotential code' '\n'
                                        '*.RRKJ3 Andrea Dal Corso"s code (old format)' '\n'
                                        'none of the above old PWscf norm-conserving format', 
                                layout = Layout(width='40px'))




icon_atomic_coord_type = Button( icon='fa-info-circle',
                                tooltip='Units for ATOMIC_POSITIONS:' '\n'
                                
                                        'alat :' '\n'
                                        'atomic positions are in cartesian coordinates, in units of the lattice parameter (either celldm(1) or A).' '\n'
                                        'If no option is specified, "alat" is assumed; not specifying units is DEPRECATED and will no longer be allowed in the future' '\n''\n'
                                        
                                        'bohr :' '\n'
                                        'atomic positions are in cartesian coordinate, in atomic units (i.e. Bohr radii)' '\n''\n'
                                        
                                        'angstrom :' '\n'
                                        'atomic positions are in cartesian coordinates, in Angstrom' '\n''\n'
                                
                                        'crystal :' '\n'
                                        'atomic positions are in crystal coordinates, i.e.in relative coordinates of the primitive lattice' '\n'
                                        'vectors as defined either in card CELL_PARAMETERS or via the ibrav + celldm / a,b,c... variables' '\n''\n'
                                
                                        'crystal_sg :' '\n'
                                        'atomic positions are in crystal coordinates, i.e. in relative coordinates of the primitive lattice.' '\n'
                                        'This option differs from the previous one because in this case only the symmetry inequivalent atoms are given.' '\n'
                                        'The variable space_group must indicate the space group number used to find the symmetry equivalent atoms.' '\n'
                                        'The other variables that control this option are uniqueb, origin_choice, and rhombohedral.'
                                
                                , 
                                layout = Layout(width='40px'))



####

form_items = [
    HBox([Label(value='Diagonalization'), 
          Box([diagonalization, icon_diagonalization])
         ], layout=form_item_layout()),
    HBox([Label(value='Mixing Beta'), 
          Box([mixing_beta, icon_mixing_beta])
          ], layout=form_item_layout()),
    HBox([Label(value='Convergence Threshold'), 
          Box([con_threshold, icon_con_threshold])
           ], layout=form_item_layout()),
    HBox([Label(value='Atom'), 
          Box([atom, icon_atom])
           ], layout=form_item_layout()),
    HBox([Label(value='Atomic Massa'), 
          Box([atomic_mass, icon_atomic_mass])
           ], layout=form_item_layout()),
    HBox([Label(value='Pseudopotential File name'), 
          Box([pseudo_filename, icon_pseudo_filename])
           ], layout=form_item_layout()),
    HBox([Label(value='Atomic Coordinate Type'),
          Box([atomic_coord_type, icon_atomic_coord_type])
           ], layout=form_item_layout()),
    HBox([Label(value='Documentation URL:'), documentation_link], layout=Layout())
]

## UI
electrons_box = Box(form_items, layout=box_layout(60))

#electrons_box

# %% [markdown]
# ## K Points

# %%
# K Points

# option for method of specifying k point grid
kpoint_type = Dropdown(name="K Point Type", value="tplba", 
                       options=["tplba", "automatic","crystal"]                   
                      )                  

# kpoints in 100 crystal direction, for sampling periodic cell in reciprocal space in the electronic structure calculation
kptx = IntText(name="kptx", value=6, layout=Layout(width='60%'))

# kpoints in 010 crystal direction, for sampling periodic cell in reciprocal space in the electronic structure calculation
kpty = IntText(name="kpty", value= 6, layout=Layout(width='60%'))

# kpoints in 001 crystal direction, for sampling periodic cell in reciprocal space in the electronic structure calculation
kptz= IntText(name="kptz", value= 6, layout=Layout(width='60%'))

# offset for kpoints in 100 crystal direction
offset_x= IntText(name="offsetx", value= 0, layout=Layout(width='60%'))

# offset for kpoints in 010 crystal direction
offset_y= IntText(name="offsety", value= 0, layout=Layout(width='60%'))

# offset for kpoints in 001 crystal direction
offset_z= IntText(name="offsetz", value= 0, layout=Layout(width='60%'))

left_section = VBox(
    [
        
        Box([Label(value='kptx'), kptx], layout=form_item_layout()),
        Box([Label(value='kpty'), kpty], layout=form_item_layout()),
        Box([Label(value='kptz'), kptz], layout=form_item_layout()),
    ])

right_section = VBox(
    [
        Box([Label(value='offsetx'), offset_x], layout=form_item_layout()),
        Box([Label(value='offsety'), offset_y], layout=form_item_layout()),
        Box([Label(value='offsetz'), offset_z], layout=form_item_layout()),
    ])


kpoints_box = VBox(
    [
        HBox([Label(value="K Point Type"), kpoint_type], layout=Layout(width='60%')),
        GridBox([left_section, right_section], 
                           layout=Layout(grid_template_columns="repeat(2,30%)", grid_gap="50px"))
    ])
    
#kpoints_box

# %% [markdown]
# ### Combine UI forms for this section

# %%
pw_x_tabs = Tab()
tab_contents = [controls_box, system_box, electrons_box, kpoints_box]
pw_x_tabs.children = tab_contents

pw_x_tabs.set_title(0, "Control")
pw_x_tabs.set_title(1, "System")
pw_x_tabs.set_title(2, "Electrons")
pw_x_tabs.set_title(3, "K Points")

#pw_x_tabs


# %% [markdown]
# ## NON USER INPUTS

# %%
#input parameters
#&control
wf_collect = '.true'
etot_conv_thr = 1e-4

#&system

nbnd_flag = 0
nbnd = 4
nosym = '.false.'

#&electrons

psName = 'Si_ONCV_PBE-1.2.upf'

cellparam_coordinate_type = 'alat'
cellparam_list = [0.0,0.0,0.0,0.25,0.25,0.25,0.5,0.5,0.5]

atom_list = ['Si',0.0,0.0,0.0,'Si',0.25,0.25,0.25]

atom_array = np.array(atom_list)
atom_splits = np.split(atom_array,nat.value)
join_line = [' '.join(i) for i in atom_splits]
atomic_positions = '\n'.join(join_line)

search_sym = {
    "desc": "option for turning off automatic symmetry search, useful for defect cells",
    "type": "Text",
    "value": 'true' 
}

pw_cores = 20
cores_per_node = 20
pw_walltime = '00:05:00'

# global variable
scf_inputs={}

# %% [markdown]
# ### Bind inputs to outputs

# %%
def bind_PW_X_inputs(self):

    scf_inputs = {
        'material_prefix': material_prefix.value,
        'restart_mode': restart_mode.value,
        'wf_collect': wf_collect,
        'pseudo_dir': pseudo_dir.value,
        'outdir': outdir.value,
        'max_seconds': max_runtime.value,
        'etot_conv_thr': etot_conv_thr,
        'ibrav': ibrav.value,
        'celldm_1': celldm_1.value,
        'nat': nat.value,
        'ntyp': ntyp.value,
        'ecutwfc': ecutwfc.value,
        'nosym': nosym,
        'diagonalization': diagonalization.value,
        'mixing_beta': mixing_beta.value,
        'conv_thr': con_threshold.value,
        'atom': atom.value,
        'atomic_mass': atomic_mass.value,
        'psName': psName,
        'coordinate_type': atomic_coord_type.value,
        'atomic_positions': atomic_positions,
        'kptx': kptx.value,
        'kpty': kpty.value,
        'kptz': kptz.value,
        'offsetx': offset_x.value,
        'offsety': offset_y.value,
        'offsetz': offset_z.value
    }
    
    if ibrav.value == 0:
        cellparam_strings = [str(i) for i in cellparam_list]
        cellparam_array = np.array(cellparam_strings)
        cellparam_splits = np.split(cellparam_array,3)
        join_line = [' '.join(i) for i in cellparam_splits]
        cell_parameters = '\n'.join(join_line)
            
        #scf_inputs["cellparam_coordinate_type"] = cellparam_coordinate_type
        #scf_inputs["cell_parameters"] = cell_parameters
        
        scf_inputs = {**scf_inputs, 'cellparam_coordinate_type': cellparam_coordinate_type, 'cell_parameters': cell_parameters}
    
    
    create_pw_file_for_scf(scf_inputs, material_prefix.value)  
    pw_scf_simulation(self)
    create_pw_file_for_nscf(scf_inputs, material_prefix.value)
    pw_nscf_simulation(self)

# %%


# %%


# %%


# %%

def create_pw_file_for_scf(scf_inputs, material_prefix):
    #print('CREATING PW SCF INPUT FILE')
    #print('now changing progress number')
    #call function to update progress bar
    global_update_1()
    status_1()
    
    #print(f'asdlfkjhasd;lfhas;dlfkjhas')
    
    #Build pw input file
    scf_name = 'scf-%s.in' % material_prefix #assigns pw input file to variable pw_name 

    if ibrav.value == 0:
        input_file = '''
        &control                                                                                                                                                                                     
        calculation = 'scf'                                                                                                                                                                         
        prefix = '{material_prefix}'                                                                                                                                                                                

        restart_mode = '{restart_mode}'                                                                                                                                                               
        wf_collect = .true.                                                                                                                                                                          
        verbosity = 'high'                                                                                                                                                                           
        outdir = './'
        max_seconds = {max_seconds}
        etot_conv_thr = {etot_conv_thr}
        pseudo_dir = './'
        /                                                                                                                                                                                            
        &system                                                                                                                                                                                      
        ibrav = {ibrav}                                                                                                                                                                                    
        celldm(1) = {celldm_1}                                                                                                                                                                           
        nat = {nat}                                                                                                                                                                                      
        ntyp = {ntyp}                                                                                                                                                                                                                                                                                                                                                                          
        ecutwfc = {ecutwfc} 
        nosym = {nosym}
        /                                                                                                                                                                                            
        &electrons                                                                                                                                                                                   
        diagonalization = '{diagonalization}'                                                                                                                                                                  
        mixing_beta = {mixing_beta}                                                                                                                                                                            
        conv_thr = {conv_thr}                                                                                                                                                                          
        /                                                                                                                                                                                            
        ATOMIC_SPECIES                                                                                                                                                                               
        {atom} {atomic_mass} {psName}                                                                                                                                                            
        ATOMIC_POSITIONS {coordinate_type}                                                                                                                                                                        
        {atomic_positions}
        K_POINTS automatic                                                                                                                                                                           
        {kptx} {kpty} {kptz} {offsetx} {offsety} {offsetz}
        CELL_PARAMETERS {cellparam_coordinate_type}
        {cell_parameters}
        '''.format(**scf_inputs) #assigns information in ''' ''' to variable inputfile

    else:
        input_file = '''
        &control                                                                                                                                                                                     
        calculation = 'scf'                                                                                                                                                                         
        prefix = '{material_prefix}'                                                                                                                                                                                

        restart_mode = '{restart_mode}'                                                                                                                                                               
        wf_collect = .true.                                                                                                                                                                          
        verbosity = 'high'                                                                                                                                                                           
        outdir = './'
        max_seconds = {max_seconds}
        etot_conv_thr = {etot_conv_thr}
        pseudo_dir = './'
        /                                                                                                                                                                                            
        &system                                                                                                                                                                                      
        ibrav = {ibrav}                                                                                                                                                                                    
        celldm(1) = {celldm_1}                                                                                                                                                                           
        nat = {nat}                                                                                                                                                                                      
        ntyp = {ntyp}                                                                                                                                                                                                                                                                                                                                                                          
        ecutwfc = {ecutwfc} 
        nosym = {nosym}
        /                                                                                                                                                                                            
        &electrons                                                                                                                                                                                   
        diagonalization = '{diagonalization}'                                                                                                                                                                  
        mixing_beta = {mixing_beta}                                                                                                                                                                            
        conv_thr = {conv_thr}                                                                                                                                                                          
        /                                                                                                                                                                                            
        ATOMIC_SPECIES                                                                                                                                                                               
        {atom} {atomic_mass} {psName}                                                                                                                                                              
        ATOMIC_POSITIONS {coordinate_type}                                                                                                                                                                        
        {atomic_positions}
        K_POINTS automatic                                                                                                                                                                           
        {kptx} {kpty} {kptz} {offsetx} {offsety} {offsetz}
        '''.format(**scf_inputs) #assigns information in ''' ''' to variable inputfile

    with open(scf_name, "w") as f: #opens file pw_name
        f.write(input_file) #writes inputfile to file pw_name
        

# %%
%use espresso-6.8
    
def pw_scf_simulation(self):
    #print('STARTED PW SCF SIMULATION')
    #print('now changing progress number again')
    #call function to update progress bar
    import time
    time.sleep(2)
    global_update_2()
    status_2()
    
    
    mat_prefix = material_prefix.value

    save_dir = f'{mat_prefix}.save'
    xml_file = f'{mat_prefix}.xml'
    outscf = f'{mat_prefix}scf'   
    psFile = pseudo_filename.value.lstrip('file://')
    
    scf_name = 'scf-%s.in' % mat_prefix #assigns pw input file to variable pw_name 
    
    
    
    
    
    #print(f'STARTED PW SIMULATION WITH {scf_name}, wall time {pw_walltime}, output file asdfsadf and cores per node {cores_per_node}')
    
    
    
    
    
    #Run pw.x simulation
    
    !submit  -w $pw_walltime -n $pw_cores -N $cores_per_node --runName=$outscf -i $psFile espresso-6.6_pw < $scf_name 
    

# %% [markdown]
# ## DFT NSCF Step (PW.X)

# %%
#nscf pw.x phase with kmesh.pl
#Build pw input file

def create_pw_file_for_nscf(scf_inputs, material_prefix):
    #print('CREATING PW NSCF INPUT FILE')
    time.sleep(2)
    global_update_3()
    status_3()
    
    
    
    nscf_name = 'nscf-%s.in' % material_prefix #assigns pw input file to variable pw_name 

    if ibrav == 0:

        input_file = '''
        &control                                                                                                                                                                                     
        calculation = 'nscf'                                                                                                                                                                         
        prefix = '{material_prefix}'                                                                                                                                                                                

        restart_mode = '{restart_mode}'                                                                                                                                                               
        wf_collect = .true.                                                                                                                                                                          
        verbosity = 'high'                                                                                                                                                                           
        outdir = './'
        pseudo_dir = '{pseudo_dir}'
        /                                                                                                                                                                                            
        &system                                                                                                                                                                                      
        ibrav = {ibrav}                                                                                                                                                                                    
        celldm(1) = {celldm_1}                                                                                                                                                                           
        nat = {nat}                                                                                                                                                                                      
        ntyp = {ntyp}                                                                                                                                                                                                                                                                                                                                                                          
        ecutwfc = {ecutwfc}
        /                                                                                                                                                                                            
        &electrons                                                                                                                                                                                   
        diagonalization = '{diagonalization}'                                                                                                                                                                  
        mixing_beta = {mixing_beta}                                                                                                                                                                            
        conv_thr = {conv_thr}                                                                                                                                                                          
        /                                                                                                                                                                                            
        ATOMIC_SPECIES                                                                                                                                                                               
        {atom} {atomic_mass} {psName}                                                                                                                                                              
        ATOMIC_POSITIONS {coordinate_type}                                                                                                                                                                        
        {atomic_positions}
        CELL_PARAMETERS {cellparam_coordinate_type}
        {cell_parameters}
        '''.format(**scf_inputs) #assigns information in ''' ''' to variable inputfile

    else:
        input_file = '''
        &control                                                                                                                                                                                     
        calculation = 'nscf'                                                                                                                                                                         
        prefix = '{material_prefix}'                                                                                                                                                                                

        restart_mode = '{restart_mode}'                                                                                                                                                               
        wf_collect = .true.                                                                                                                                                                          
        verbosity = 'high'                                                                                                                                                                           
        outdir = './'
        pseudo_dir = '{pseudo_dir}'
        /                                                                                                                                                                                            
        &system                                                                                                                                                                                      
        ibrav = {ibrav}                                                                                                                                                                                    
        celldm(1) = {celldm_1}                                                                                                                                                                           
        nat = {nat}                                                                                                                                                                                      
        ntyp = {ntyp}                                                                                                                                                                                                                                                                                                                                                                          
        ecutwfc = {ecutwfc}
        /                                                                                                                                                                                            
        &electrons                                                                                                                                                                                   
        diagonalization = '{diagonalization}'                                                                                                                                                                  
        mixing_beta = {mixing_beta}                                                                                                                                                                            
        conv_thr = {conv_thr}                                                                                                                                                                          
        /                                                                                                                                                                                            
        ATOMIC_SPECIES                                                                                                                                                                               
        {atom} {atomic_mass} {psName}                                                                                                                                                              
        ATOMIC_POSITIONS {coordinate_type}                                                                                                                                                                        
        {atomic_positions}
        '''.format(**scf_inputs) #assigns information in ''' ''' to variable inputfile

    with open(nscf_name, "w") as f: #opens file pw_name
        f.write(input_file) #writes inputfile to file pw_na

# %%
EXTRA_FILES = ["kmesh.pl"]

# %%
#Run pw.x simulation 
%use espresso-6.8
    
def pw_nscf_simulation(self):
    #print('STARTED PW SCF SIMULATION')
    time.sleep(2)
    global_update_4()
    status_4()
    
    mat_prefix = material_prefix.value

    save_dir = f'{mat_prefix}.save'
    xml_file = f'{mat_prefix}.xml'
    
    outnscf = f'{mat_prefix}nscf'  
    psFile = pseudo_filename.value.lstrip('file://')
    
    nscf_name = 'nscf-%s.in' % mat_prefix #assigns pw input file to variable pw_name 
    
    #Run pw.x simulation  
    kptx_var = kptx.value
    kpty_var = kpty.value
    kptz_var = kptz.value
    
    !perl kmesh.pl $kptx_var $kpty_var $kptz_var >> $nscf_name
    
    if nbnd_flag == 0:
        
        !submit -w $pw_walltime -n $pw_cores -N $cores_per_node --runName=$outnscf -i $psFile -i $save_dir -i $xml_file espresso-6.6_pw -npool $pw_cores < $nscf_name
        
    if nbnd_flag == 1:
        scf_inputs['nbnd'] = nbnd

    else:
        #find appropriate number of bands and eventual Wannier functions for optical calculation, rule of thumb being 2*# of valence bands
        with open(f'{pwd}/{outnscf}.stdout') as file:
            for line in file:
                j = line.split()
                for i in j:
                    if i == 'Kohn-Sham':
                        KSstates = float(j[len(j)-1])
                        print(KSstates)
                        nbnd = KSstates*2
        scf_inputs['nbnd'] = int(nbnd)      

# %%


# %% [markdown]
# ### Functions to Update Progress Bar

# %%
##('CREATING PW SCF INPUT FILE')
def global_update_1():
    global progressNum
    progressNum=1

    

# %%
##('STARTED PW SCF SIMULATION')
def global_update_2():
    global progressNum
    progressNum=2


# %%
##('CREATING PW NSCF INPUT FILE')
def global_update_3():
    global progressNum
    progressNum=3


# %%
##('STARTED PW SCF SIMULATION')
def global_update_4():
    global progressNum
    progressNum=4


# %% [markdown]
# ### Functions to Display and Update Status Messages

# %%
##('CREATING PW SCF INPUT FILE')
def status_1():
    global x
    x=1
    global flag
    flag=1
    

# %%
##('STARTED PW SCF SIMULATION')
def status_2():
    global x
    x=2
    global flag
    flag=1
    

# %%
##('CREATING PW NSCF INPUT FILE')
def status_3():
    global x
    x=3
    global flag
    flag=1
    

# %%
##STARTED PW SCF SIMULATION
def status_4():
    global x
    x=4
    global flag
    flag=1
    


