# %% [markdown]
# # DFPT Phonon Step (PH.X)

# %%
### Import libraries
from ipywidgets import Style, Box, Layout, Label, IntText, FloatText, Dropdown
%run styles.ipynb

# %% [markdown]
# ## Documentation Link

# %%
%%html
<a href="https://www.quantum-espresso.org/Doc/INPUT_PH.html#idm218" target="_blank">INPUT_PW.html</a>

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

ph_x_box = Box(form_items, layout=box_layout(40))

ph_x_box


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

    build_PH_Input(ph_inputs, material_prefix.value)

# %%
def build_PH_Input(ph_inputs, material_prefix):
    print('BUILDING PH INPUT FILE')
    ph_name = 'ph-%s.in' % material_prefix #assigns ph input file to variable ph_name 

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


