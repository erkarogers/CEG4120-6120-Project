# %% [markdown]
# # DFPT Phonon Step (PH.X)

##
# @file ph_x.py
# 
# @brief PH_X step input widget creation
# 
# @page PH_X
# @section description_ph Description
# @ref ph_x.py is is a jupyter notebook that contains the input widgets for the PH_X step
# 
# @section todo_ph TODO
# 
# @section notes_ph Notes
# @ref Quantum Espresso Input_PH documentation: https://www.quantum-espresso.org/Doc/INPUT_PH.html#idm218
#
# @section libraries_ph Libraries/Modules
# - ipywidgets  : https://ipywidgets.readthedocs.io/en/latest/
#  - Style      : https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20Styling.html#Styling-Widgets
#  - Box        : https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20List.html#Box
#  - Layout     : https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20List.html#Layout
#  - Label      : https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20List.html#Label
#  - IntText    : https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20List.html#IntText
#  - FloatText  : https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20List.html#FloatText
#  - Dropdown   : https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20List.html#Dropdown
# - os          : https://docs.python.org/3/library/os.html
# - nanohub     : https://nanohub.org/tools/quantumespresso/
# - styles      : styles.py
#  - Runs the styles.ipynb notebook 
# 


# %%
### Import libraries
from ipywidgets import Style,HBox, Box, Layout, HTML, Button, Label, IntText, FloatText, Dropdown
import hublib.use
import time
%run styles.ipynb

# %% [markdown]
# ## Documentation Link

# %% [markdown]
# Check out the documentation on inputs variables needed for this step with this link:
# <a href="https://www.quantum-espresso.org/Doc/INPUT_PH.html#idm218" target="_blank">INPUT_PW.html</a>

# %%
documentation_link = HTML(value='<a href="https://www.quantum-espresso.org/Doc/INPUT_PH.html#idm218" target="_blank">INPUT_PH.html</a>')

# %% [markdown]
# ## PH.X Inputs

# %%
# ph.x inputs

#qpoints in 100 crystal direction, for sampling periodic cell in reciprocal space in the phonon calculation
nq1 = IntText(name="nq1", value= 3)

#qpoints in 010 crystal direction, for sampling periodic cell in reciprocal space in the phonon calculation
nq2 = IntText(name="nq2", value= 3)
# qpoints in 001 crystal direction, for sampling periodic cell in reciprocal space in the phonon calculation
nq3 = IntText(name="nq3", value= 3)

# convergence threshold for phonon calculation
self_consistency_threshold = FloatText(name="Self-Consistency Threshold", value=1.0e-12)

# option for recovering from previous failed/timed out calculation
recover = Dropdown(name="recover", value="false", options=["true","false"])


####
# Creates the info buttons and adds in the description when you hover your cursor over the button

icon_nq1 = Button( icon='fa-info-circle',
                                tooltip='Parameters of the Monkhorst-Pack grid (no offset) used' '\n'
                                        'when ldisp=.true. Same meaning as for nk1, nk2, nk3 in the input of pw.x.' '\n'
                                         '\n''DEFAULT: 0,0,0' , 
                                layout = Layout(width='40px'))

icon_nq2 = Button( icon='fa-info-circle',
                                tooltip='Parameters of the Monkhorst-Pack grid (no offset) used' '\n'
                                        'when ldisp=.true. Same meaning as for nk1, nk2, nk3 in the input of pw.x.''\n'
                                          '\n''DEFAULT: 0,0,0', 
                                layout = Layout(width='40px'))

icon_nq3 = Button( icon='fa-info-circle',
                                tooltip='Parameters of the Monkhorst-Pack grid (no offset) used' '\n'
                                        'when ldisp=.true. Same meaning as for nk1, nk2, nk3 in the input of pw.x.''\n'
                                         '\n' 'DEFAULT: 0,0,0', 
                                layout = Layout(width='40px'))

icon_self_consistency_threshold = Button( icon='fa-info-circle',
                                tooltip='Threshold for self-consistency.' '\n'
                                         '\n''DEFAULT:1e-12', 
                                layout = Layout(width='40px'))

icon_recover = Button( icon='fa-info-circle',
                                tooltip=' If .true. restart from an interrupted run.' '\n'
                                          '\n''DEFAULT: false', 
                                layout = Layout(width='40px'))





form_items = [
    
    HBox([Label(value='nq1'),
          Box([nq1, icon_nq1])
         ], layout = form_item_layout()),
        
    HBox([Label(value='nq2'), 
          Box([nq2, icon_nq2])
         ],layout=form_item_layout()),
    
     HBox([Label(value='nq3'), 
          Box([nq3, icon_nq3])
         ],layout=form_item_layout()),
    
     HBox([Label(value='Self-Consistency Threshold'), 
          Box([self_consistency_threshold, icon_self_consistency_threshold])
         ],layout=form_item_layout()),
    
      HBox([Label(value='Recover'), 
          Box([recover, icon_recover])
         ],layout=form_item_layout()),   
    
    HBox([Label(value='Documentation URL:'), documentation_link], layout=Layout())
    
    #Box([Label(value='nq1'), nq1], layout=form_item_layout()),
    #Box([Label(value='nq2'), nq2], layout=form_item_layout()),
    #Box([Label(value='nq3'), nq3], layout=form_item_layout()),
    #Box([Label(value='Self-Consistency Threshold'), self_consistency_threshold], layout=form_item_layout()),
    #Box([Label(value='Recover'), recover], layout=form_item_layout())
]

ph_x_box = Box(form_items, layout=box_layout(45))

#ph_x_box


# %% [markdown]
# ## NON USER INPUTS

# %%
#input parameters
#&inputph                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
search_sym = '.false.'   
ph_cores = 20
ph_walltime = '00:05:00'

# %% [markdown]
# ## Bind inputs to outputs

# %%
def bind_PH_X_inputs(self):
    
    ph_inputs = {
        'material_prefix': material_prefix.value,
        'nq1': nq1.value,
        'nq2': nq2.value,
        'nq3': nq3.value,
        'tr2_ph': self_consistency_threshold.value,
        'recover': recover.value,
        'search_sym': search_sym
    }

    
    create_ph_file(ph_inputs, material_prefix.value)
    ph_simulation(self)
    get_qpoints(self)
    runPP(self)

# %%
def create_ph_file(ph_inputs, material_prefix):
    #print('CREATING PH INPUT FILE')
    
    time.sleep(2)
    global_update_5()
    status_5()
    
    
    
    
    
    ph_name = 'ph-%s.in' % material_prefix #assigns ph input file to variable ph_name 

    input_file = '''
    --
    &inputph
    prefix = '{material_prefix}'
    fildvscf = 'dvscf'
    ldisp = .true.
    fildyn = '{material_prefix}.dyn'
    nq1 = {nq1}
    nq2 = {nq2}
    nq3 = {nq3}
    tr2_ph = {tr2_ph}
    recover = {recover}
    search_sym = {search_sym}
    /
    '''.format(**ph_inputs) #assigns information in ''' ''' to variable inputfile

    with open(ph_name, "w") as f: #opens file pw_name
        f.write(input_file) #writes inputfile to file ph

# %%
%use espresso-6.8
    
def ph_simulation(self):
    #print('STARTED PH SIMULATION')
    
    time.sleep(2)
    global_update_6()
    status_6()
    
    
    
    mat_prefix = material_prefix.value

    save_dir = f'{mat_prefix}.save'
    xml_file = f'{mat_prefix}.xml'

    outph = f'{mat_prefix}ph'    
    psFile = pseudo_filename.value.lstrip('file://')
    
    ph_name = 'ph-%s.in' % mat_prefix #assigns ph input file to variable ph_name 
    
    #print(f'STARTED PH SIMULATION WITH {ph_name}')
    #Run ph.x simulation
   
    !submit -w $ph_walltime -n $ph_cores -N $cores_per_node --runName=$outph -i $psFile -i $save_dir -i $xml_file espresso-6.6_ph < $ph_name
    

# %% [markdown]
# ## PH Workflow

# %% [markdown]
# ### List of q-points

# %%
qpoint_list_append = ""
def get_qpoints(self):
    global qpoint_list_append
    #obtain list of q-points from {material_prefix}.dyn0, to be appended to the end of the eventual epw.in file
    with open(f'{pwd}/{material_prefix.value}.dyn0') as file:
        qpoint_list = file.readlines()
        del qpoint_list[0]
        num_irr_q = str(qpoint_list[0]).strip()
        qpoint_list[0]= num_irr_q+' cartesian'
        #print(qpoint_list)

    i = 0 
    while i < len(qpoint_list):
        qpoint_list_append = qpoint_list_append+qpoint_list[i]+'\n'
        i += 1

# %% [markdown]
# ### Run pp.y

# %%
#run pp.py
# Post-processing script from of PH data in format used by EPW
# 14/07/2015 - Creation of the script - Samuel Ponce
# 14/03/2018 - Automatically reads the number of q-points - Michael Waters
# 14/03/2018 - Detect if SOC is included in the calculation - Samuel Ponce
# 05/06/2019 - Removed SOC for xml detection instead - Felix Goudreault
#




    
    
    

from __future__ import print_function
try:
    from builtins import input
except ImportError:
    print('Install future. e.g. "pip install --user future"')
# import numpy as np

import os
from xml.dom import minidom


# Return the number of q-points in the IBZ
def get_nqpt(prefix):
    
    #STARTED POST PROCESSING PH SIMULATION RESULTS
    time.sleep(2)
    global_update_7()
    status_7()
    
    
    
    
    
    
    fname = '_ph0/' + prefix + '.phsave/control_ph.xml'

    fid = open(fname, 'r')
    lines = fid.readlines()
    # these files are relatively small so reading the whole thing shouldn't
    # be an issue
    fid.close()

    line_number_of_nqpt = 0
    while 'NUMBER_OF_Q_POINTS' not in lines[line_number_of_nqpt]:
        # increment to line of interest
        line_number_of_nqpt += 1
    line_number_of_nqpt += 1  # its on the next line after that text

    nqpt = int(lines[line_number_of_nqpt])

    return nqpt


# Check if the calculation include SOC
def hasSOC(prefix):
    fname = prefix+'.save/data-file-schema.xml'

    xmldoc = minidom.parse(fname)
    item = xmldoc.getElementsByTagName('spinorbit')[0]
    lSOC = item.childNodes[0].data

    return lSOC


# Check if the calculation includes PAW
def hasPAW(prefix):
    fname = prefix+'.save/data-file-schema.xml'

    xmldoc = minidom.parse(fname)
    item = xmldoc.getElementsByTagName('paw')[0]
    lPAW = (item.childNodes[0].data == 'true')

    return lPAW


# Check if the calculation used .fc or .fc.xml files
def hasfc(prefix):
    fname = str(prefix)+'.fc.xml'
    if (os.path.isfile(fname)):
        lfc = True
    else:
        fname_no_xml = fname.strip(".xml")
        if (os.path.isfile(fname_no_xml)):
            lfc = True
        else:
            lfc = False

    return lfc


# check if calculation used xml files (irrelevant of presence of SOC)
def hasXML(prefix):
    # check for a file named prefix.dyn1.xml
    # if it exists => return True else return False
    fname = os.path.join(prefix + ".dyn1.xml")
    if os.path.isfile(fname):
        return True
    # check if the other without .xml extension exists
    # if not raise an error
    fname_no_xml = fname.strip(".xml")

    class FileNotFoundError(Exception):
        pass
    if not os.path.isfile(fname_no_xml):
        raise FileNotFoundError(
                "No dyn0 file found cannot tell if xml format was used.")
    return False


# Check if the calculation was done in sequential
def isSEQ(prefix):
    fname = '_ph0/'+str(prefix)+'.dvscf'
    if (os.path.isfile(fname)):
        lseq = True
    else:
        lseq = False

    return lseq


# Enter the number of irr. q-points
#user_input = input(
#        'Enter the prefix used for PH calculations (e.g. diam)\n')


def runPP(self):
    print('RUNNING pp.py')
    prefix = material_prefix.value

    # # Test if SOC
    # SOC = hasSOC(prefix)
    # Test if '.xml' files are used
    XML = hasXML(prefix)

    # Test if PAW
    PAW = hasPAW(prefix)

    # Test if fc
    fc = hasfc(prefix)

    # Test if seq. or parallel run
    SEQ = isSEQ(prefix)

    if True:  # this gets the nqpt from the outputfiles
        nqpt = get_nqpt(prefix)

    else:
        # Enter the number of irr. q-points
        user_input = input(
                'Enter the number of irreducible q-points\n')
        nqpt = user_input
        try:
            nqpt = int(user_input)
        except ValueError:
            raise Exception('The value you enter is not an integer!')

    os.system('mkdir save 2>/dev/null')

    for iqpt in range(1, nqpt+1):
        label = str(iqpt)

        # Case calculation in seq.
        if SEQ:
            # Case with XML files
            if XML:
                os.system('cp '+prefix+'.dyn0 '+prefix+'.dyn0.xml')
                os.system('cp '+prefix+'.dyn'+str(iqpt)+'.xml save/'+prefix
                          + '.dyn_q'+label+'.xml')
                if (iqpt == 1):
                    os.system('cp _ph0/'+prefix+'.dvscf* save/'+prefix+'.dvscf_q'
                              + label)
                    os.system('cp -r _ph0/'+prefix+'.phsave save/')
                    if fc:
                        os.system('cp '+prefix+'.fc.xml save/ifc.q2r.xml')
                    if PAW:
                        os.system('cp _ph0/'+prefix+'.dvscf_paw* save/'+prefix +
                                  '.dvscf_paw_q'+label)
                else:
                    os.system('cp _ph0/'+prefix+'.q_'+str(iqpt)+'/'+prefix +
                              '.dvscf* save/'+prefix+'.dvscf_q'+label)
                    os.system('rm _ph0/'+prefix+'.q_'+str(iqpt)+'/*wfc*')
                    if PAW:
                        os.system('cp _ph0/'+prefix+'.q_'+str(iqpt)+'/'+prefix +
                                  '.dvscf_paw* save/'+prefix+'.dvscf_paw_q'+label)
            # Case without XML files
            else:
                os.system('cp '+prefix+'.dyn'+str(iqpt)+' save/'+prefix+'.dyn_q' +
                          label)
                if (iqpt == 1):
                    os.system('cp _ph0/'+prefix+'.dvscf save/'+prefix+'.dvscf_q' +
                              label)
                    os.system('cp -r _ph0/'+prefix+'.phsave save/')
                    if fc:
                        os.system('cp '+prefix+'.fc save/ifc.q2r')
                    if PAW:
                        os.system('cp _ph0/'+prefix+'.dvscf_paw save/'+prefix +
                                  '.dvscf_paw_q'+label)
                else:
                    os.system('cp _ph0/'+prefix+'.q_'+str(iqpt)+'/'+prefix +
                              '.dvscf save/'+prefix+'.dvscf_q'+label)
                    os.system('rm _ph0/'+prefix+'.q_'+str(iqpt)+'/*wfc*')
                    if PAW:
                        os.system('cp _ph0/'+prefix+'.q_'+str(iqpt)+'/'+prefix +
                                  '.dvscf_paw save/'+prefix+'.dvscf_paw_q'+label)
        else:
            # Case with XML format
            if XML:
                os.system('cp '+prefix+'.dyn0 '+prefix+'.dyn0.xml')
                os.system('cp '+prefix+'.dyn'+str(iqpt)+'.xml save/'+prefix +
                          '.dyn_q'+label+'.xml')
                if (iqpt == 1):
                    os.system('cp _ph0/'+prefix+'.dvscf1 save/'+prefix+'.dvscf_q' +
                              label)
                    os.system('cp -r _ph0/'+prefix+'.phsave save/')
                    if fc:
                        os.system('cp '+prefix+'.fc.xml save/ifc.q2r.xml')
                    if PAW:
                        os.system('cp _ph0/'+prefix+'.dvscf_paw1 save/'+prefix +
                                  '.dvscf_paw_q'+label)
                else:
                    os.system('cp _ph0/'+prefix+'.q_'+str(iqpt)+'/'+prefix +
                              '.dvscf1 save/'+prefix+'.dvscf_q'+label)
                    os.system('rm _ph0/'+prefix+'.q_'+str(iqpt)+'/*wfc*')
                    if PAW:
                        os.system('cp _ph0/'+prefix+'.q_'+str(iqpt)+'/'+prefix +
                                  '.dvscf_paw1 save/'+prefix+'.dvscf_paw_q'+label)
            # Case without XML format
            else:
                os.system('cp '+prefix+'.dyn'+str(iqpt)+' save/'+prefix+'.dyn_q' +
                          label)
                if (iqpt == 1):
                    os.system('cp _ph0/'+prefix+'.dvscf1 save/'+prefix+'.dvscf_q' +
                              label)
                    os.system('cp -r _ph0/'+prefix+'.phsave save/')
                    if fc:
                        os.system('cp '+prefix+'.fc save/ifc.q2r')
                    if PAW:
                        os.system('cp _ph0/'+prefix+'.dvscf_paw1 save/'+prefix +
                                  '.dvscf_paw_q'+label)
                else:
                    os.system('cp _ph0/'+prefix+'.q_'+str(iqpt)+'/'+prefix +
                              '.dvscf1 save/'+prefix+'.dvscf_q'+label)
                    os.system('rm _ph0/'+prefix+'.q_'+str(iqpt)+'/*wfc*')
                    if PAW:
                        os.system('cp _ph0/'+prefix+'.q_'+str(iqpt)+'/'+prefix +
                                  '.dvscf_paw1 save/'+prefix+'.dvscf_paw_q'+label)

# %%


# %% [markdown]
# ### Functions to Update Progress Bar

# %%
##('CREATING PH INPUT FILE')
def global_update_5():
    global progressNum
    progressNum=5

    

# %%
##'STARTED PH SIMULATION'
def global_update_6():
    global progressNum
    progressNum=6

    

# %%
##STARTED POST PROCESSING OF PH SIMULATION RESULTS
def global_update_7():
    global progressNum
    progressNum=7

    

# %% [markdown]
# ### Functions to Display and Update Status Messages

# %%
##('CREATING PH INPUT FILE')
def status_5():
    global x
    x=5
    global flag
    flag=1
    

# %%
##('STARTED PH SIMULATION')
def status_6():
    global x
    x=6
    global flag
    flag=1
    

# %%
##('STARTED POST PROCESSING OF PH SIMULATION RESULTS')
def status_7():
    global x
    x=7
    global flag
    flag=1
    


