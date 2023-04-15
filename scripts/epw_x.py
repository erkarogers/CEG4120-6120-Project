# %% [markdown]
# # EPW Step (EPW.x)

##
# @file epw_x.py
# 
# @brief EPW_X step input widget creation
# 
# @page EPW_X
# @section description_epw Description
# @ref EPW_X.py is a Jupyter Notebook that creates the input widgets for the EPW step of the simulation tool.  The widgets are created using the ipywidgets library.  The widgets are then displayed in a tabbed format using the AppLayout widget.  The widgets are then bound to the output widget.  The output widget is then displayed in the AppLayout widget.  The AppLayout widget is then displayed in the Jupyter Notebook.
# 
# @section todo_epw TODO
# 
# @section notes_epw Notes
# @ref Quantum Espresso EPW documentation: https://docs.epw-code.org/doc/Inputs.html#filkf
# 
# @section libraries_epw Libraries/Modules
# - ipywidgets  : https://ipywidgets.readthedocs.io/en/latest/
#  - Tab         : https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20List.html#Tab
#  - Box         : https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20List.html#Box
#  - VBox        : https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20List.html#VBox
#  - GridBox     : https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20List.html#GridBox
#  - Layout      : https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20List.html#Layout
#  - Label       : https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20List.html#Label
#  - IntText     : https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20List.html#IntText
#  - FloatText   : https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20List.html#FloatText
#  - Dropdown    : https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20List.html#Dropdown
#  - Text        : https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20List.html#Text
#  - Textarea    : https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20List.html#Textarea
# - os          : https://docs.python.org/3/library/os.html
# - styles      : styles.ipynb
# - Runs the styles.ipynb notebook
# 


# %%
### Import libraries
from ipywidgets import Tab, Box, VBox, GridBox, Layout
from ipywidgets import Label, IntText, FloatText, Dropdown, Text, Textarea
import os
%run styles.ipynb

# %% [markdown]
# ## Documentation Link

# %%
%%html
<a href="https://docs.epw-code.org/doc/Inputs.html#filkf" target="_blank">INPUT_EPW.html</a>

# %%
documentation_link = HTML(value='<a href="https://docs.epw-code.org/doc/Inputs.html#filkf" target="_blank">INPUT_EPW.html</a>')

# %% [markdown]
# ## Setup

# %%
# Setup

## 
# @brief Option for turning on/off reading in kmaps from file
# @remark Creates the kmaps input dropdown for the EPW.x step
kmaps = Dropdown(name="kmaps", value=".false.", options={"false": ".false.", "true": ".true."})

## 
# @brief Option for turning on/off writing coarse bloch space electron phonon matrix elements (.epb files) to file
# @remark Creates the epbwrite input dropdown for the EPW.x step
# 
epbwrite = Dropdown(name="epbwrite", value=".true.", options={"false": ".false.", "true": ".true."})
    
## 
# @brief Option for turning on/off reading coarse bloch space electron phonon matrix elements (.epb files) from file
# @remark Creates the epbread input dropdown for the EPW.x step
epbread = Dropdown(name="epbread", value=".false.", options={"false": ".false.", "true": ".true."})

## 
# @brief Option for turning on/off writing coarse wannier space electron phonon matrix elements (.epw files) to file
# @remark Creates the epwwrite input dropdown for the EPW.x step
epwwrite = Dropdown(name="epwwrite", value=".true.", options={"false": ".false.", "true": ".true."})

## 
# @brief epwread option for turning on/off reading coarse wannier space electron phonon matrix elements (.epw files) from file
# @remark Creates the epwread input dropdown for the EPW.x step
epwread = Dropdown(name="epwread", value=".false.", options={"false": ".false.", "true": ".true."})

## 
# @brief Option for polar material correction to Wannier interpolation
# @remark Creates the lpolar input dropdown for the EPW.x step
lpolar = Dropdown(name="Correct for polar materials",value="false",options=["false", "true"])

## 
# @brief Create Wannier conversion 
# @remark Creates the wannierize input dropdown for the EPW.x step
wannierize = Dropdown(name="Center Wannier Functions",value="true",options=["false", "true"])


## 
# 
# Creates the info buttons and adds in the description when you hover your cursor over the button

icon_kmaps = Button( icon='fa-info-circle',
                                tooltip='Generate the map k+q –> k for folding the rotation matrix U(k+q).' '\n''\n'
                                        'If .true., the program reads ‘prefix.kmap’ and ‘prefix.kgmap’ from file. If .false., they are calculated.' '\n''\n'
                                        'Note that for a restart with epwread =.true., kmaps also needs to be set to true (since the information to potentially calculate kgmaps is not generated in a restart run).''\n'
                                        'However, the files “prefix.kmap” and “prefix.kgmap” themselves are actually not used if epwread=.true. and hence need not actually be there.', 
                                layout = Layout(width='40px'))

icon_epbwrite = Button( icon='fa-info-circle',
                                tooltip='If epbwrite = .true., the electron-phonon matrix elements in the coarse Bloch representation and relevant data (dyn matrices) are written to disk.''\n''\n' 
                                           'If epbread = .true. the above quantities are read from the ‘prefix.epb’ files. Pool dependent files.', 
                                layout = Layout(width='40px'))

icon_epbread = Button( icon='fa-info-circle',
                                tooltip='If epbwrite = .true., the electron-phonon matrix elements in the coarse Bloch representation and relevant data (dyn matrices) are written to disk.' '\n''\n' 
                                           'If epbread = .true. the above quantities are read from the ‘prefix.epb’ files. Pool dependent files.', 
                                layout = Layout(width='40px'))

icon_epwwrite = Button( icon='fa-info-circle',
                                tooltip='If epwwrite = .true., the electron-phonon matrix elements in the coarse Wannier representation and relevant data (dyn matrices) are written to disk.''\n'
                                           'Each pool reads the same file.', 
                                layout = Layout(width='40px'))

icon_epwread = Button( icon='fa-info-circle',
                                tooltip='If epwread = .true., the electron-phonon matrix elements in the coarse Wannier representation are read from the ‘epwdata.fmt’ and ‘XX.epmatwpX’ files.''\n'
                                          'Each pool reads the same file. It is used for a restart calculation and requires kmaps = .true. A prior calculation with epwwrite = .true is also required.', 
                                layout = Layout(width='40px'))


icon_lpolar = Button( icon='fa-info-circle',
                                tooltip='If .true. enable the correct Wannier interpolation in the case of polar material.''\n'
                                         'DEFAULT: false', 
                                layout = Layout(width='40px'))


icon_wannierize = Button( icon='fa-info-circle',
                                tooltip=' Calculate the Wannier functions using W90 library calls and write rotation matrix to file ‘filukk’. If .false., filukk is read from disk.' '\n'
                                         'DEFAULT: false', 
                                layout = Layout(width='40px'))



form_items = [
    
    
     HBox([Label(value='kmaps'),
          Box([kmaps, icon_kmaps])
         ], layout = form_item_layout()),
        
    HBox([Label(value='epbwrite'), 
          Box([epbwrite, icon_epbwrite])
         ],layout=form_item_layout()),
    
    HBox([Label(value='epbread'), 
          Box([epbread, icon_epbread])
         ],layout=form_item_layout()),
    
     HBox([Label(value='epwwrite'), 
          Box([epwwrite, icon_epwwrite])
         ],layout=form_item_layout()),
    
     HBox([Label(value='epwread'), 
          Box([epwread, icon_epwread])
         ],layout=form_item_layout()),
    
     HBox([Label(value='lpolar'), 
          Box([lpolar, icon_lpolar])
         ],layout=form_item_layout()),
    
     HBox([Label(value='wannierize'), 
          Box([wannierize, icon_wannierize])
         ],layout=form_item_layout()),
        
    
    HBox([Label(value='Documentation URL:'), documentation_link], layout=Layout())
    
    
    
    #Box([Label(value='kmaps'), kmaps], layout=form_item_layout()),
    #Box([Label(value='epbwrite'), epbwrite], layout=form_item_layout()),
    #Box([Label(value='epbread'), epbread], layout=form_item_layout()),
    #Box([Label(value='epwwrite'), epwwrite], layout=form_item_layout()),
    #Box([Label(value='Correct for polar materials'), lpolar], layout=form_item_layout()),
    #Box([Label(value='Center Wannier Functions'), wannierize], layout=form_item_layout())
]

setup_box = Box(form_items, layout=box_layout(40))
setup_box

# %% [markdown]
# ## Wannier

# %%
# Wannier

## number of iterations for creating wannier function representations
num_iter = IntText(name="Number of Iterations", value=1500)

## verbosity of EPW output file
iprint = IntText(name="Verbosity Level", value=2)

## Maximum value of the disentanglement window. See wannier90 documentation.
dis_win_max = FloatText(name="Disentaglement Window Max", value=18.0)

## Window which includes frozen states for Wannier90. See wannier90 documentation.
dis_froz_max = FloatText(name="Window Max for Frozen States",value=8.5)

## Initial wannier projections, to be passed to Wannier90. These must agree with nbnd if specified other than the default 'random' setting (number of desired wannier projections = number of computed bands)
projections = Textarea(name="Projections?", value="proj(1) = 'random'")



####
## Creates the info buttons and adds in the description when you hover your cursor over the button

icon_num_iter = Button( icon='fa-info-circle',
                                tooltip='Number of iterations to produce maximally localized wannier functions', 
                                layout = Layout(width='40px'))

icon_iprint = Button( icon='fa-info-circle',
                                tooltip='This indicates the level of verbosity of the output from 0 (“low”), the bare minimum, to 3 (“high”), which corresponds to full debugging output.' '\n' '\n'
                                        'The default value is 1.', 
                                layout = Layout(width='40px'))

icon_dis_win_max = Button( icon='fa-info-circle',
                                tooltip='The upper bound of the outer energy window for the disentanglement procedure. Units are eV.''\n' '\n'
                                        'The default is the highest eigenvalue in the given states (i.e., all states are included in the disentanglement procedure).', 
                                layout = Layout(width='40px'))

icon_dis_froz_max = Button( icon='fa-info-circle',
                                tooltip='The upper bound of the inner (frozen) energy window for the disentanglement procedure.''\n'
                                       'If dis_froz_max is not specified, then there are no frozen states. Units are eV.' '\n' '\n'
                                        'No default.', 
                                layout = Layout(width='40px'))


icon_projections = Button( icon='fa-info-circle',
                                tooltip='', 
                                layout = Layout(width='40px'))

form_items = [
    
    
    HBox([Label(value='num_iter'),
          Box([num_iter, icon_num_iter])
         ], layout = form_item_layout()),
    
     HBox([Label(value='iprint'),
          Box([iprint, icon_iprint])
         ], layout = form_item_layout()),   
    
      HBox([Label(value='dis_win_max'),
          Box([dis_win_max, icon_dis_win_max])
         ], layout = form_item_layout()),    
    
        HBox([Label(value='dis_froz_max'),
          Box([dis_froz_max, icon_dis_froz_max])
         ], layout = form_item_layout()),    
      
         HBox([Label(value='projections'),
          Box([projections, icon_projections])
         ], layout = form_item_layout()),       
    
    HBox([Label(value='Documentation URL:'), documentation_link], layout=Layout())
    
    
    #Box([Label(value='Number of Iterations'), num_iter], layout=form_item_layout()),
    #Box([Label(value='Verbosity Level'), iprint], layout=form_item_layout()),
    #Box([Label(value='Disentaglement Window Max'), dis_win_max], layout=form_item_layout()),
    #Box([Label(value='Window Max for Frozen States'), dis_froz_max], layout=form_item_layout()),
    #Box([Label(value='Projections?'), projections], layout=form_item_layout())
]

wannier_box = Box(form_items, layout=box_layout(40))

wannier_box


# %% [markdown]
# ## Misc

# %%
# Misc

## option to print all electron-phonon coupling elements to output file
## will create a masssive file that will take a long time to write, use caution when setting to .true.
prgtkk = Dropdown(name="Print Electron-Photon Vertexs", value=".false.", 
                  options={"false": ".false.", "true": ".true."})

## option to calculate optical absorption function
lindabs = Dropdown(name="Calculate Optical Parameters", value=".true.", 
                   options={"false": ".false.", "true": ".true."})
                                                                                                                                                                           
## scissor shift to correct for DFT Bandgap narrowing
scissor = FloatText(name="Scissor Shift", value=0)

## lowest optical frequency of interest, in eV
omegamin = FloatText(name="Min Photon Energy", value=0.05 )

## highest optical frequency of interest, in eV
omegamax = FloatText( name="Max Photon Energy", value=0.05)

## increment for sweeping optical frequencies, in eV
omegastep = FloatText(name="Steps in Photon Energy", value=0.00)
                                                                                                                                                                            
## material's index of refraction
n_r = FloatText(name="Refractive Index", value=3.4)

## Width of the Fermi surface window to take into account states in the self-energy delta functions in eV. Narrowing this value reduces the FloatText of bands included in the selfenergy calculations.
fsthick = FloatText(name="Band Gap Width", value=4.0)
 
## system temperature in Kelvin
temps = FloatText(name="Temperature (K)", value=300)

## Smearing in the energy-conserving delta functions in eV
degaussw = FloatText(name="Step Function Broadening Parameter", value=0.005)

degaussq = FloatText(name="User Specific Fermi Energy", value=0.05)

## option to specify fermi energy, such as from prior nscf step
efermi_read = Dropdown(name="Fermi Energy", value=".true.",
                       options={"false": ".false.", "true": ".true."})



####
## Creates the info buttons and adds in the description when you hover your cursor over the button

icon_prgtkk = Button( icon='fa-info-circle',
                                tooltip='Allows to print the electron-phonon vertex |g| (in meV) for each q-point, k-point, i-band, j-band and modes.''\n'
                                        'Note: Average over degenerate i-band, j-band and modes is performed but not on degenerate k or q-points.''\n'
                                        'Warning: this produces huge text data in the main output file and considerably slows down the calculation.''\n'
                                        'Suggestion: Use only 1 k-point (like Gamma).' '\n' '\n'
                                         'DEFAULT: false', 
                                layout = Layout(width='40px'))

icon_lindabs = Button( icon='fa-info-circle',
                                tooltip='If .true. computes indirect phonon absorption. See the input variables omegamax, omegamin, omegastep and n_r.''\n''\n'
                                          'DEFUALT: false', 
                                layout = Layout(width='40px'))

icon_scissor = Button( icon='fa-info-circle',
                                tooltip='Gives the value of the scissor shift of the gap (in eV).''\n''\n'
                                          'DEFAULT: 0.0', 
                                layout = Layout(width='40px'))

icon_omegamin = Button( icon='fa-info-circle',
                                tooltip='Photon energy minimum (in eV) when lindabs = .true.''\n''\n'
                                           'DEFAULT: 0', 
                                layout = Layout(width='40px'))

icon_omegamax = Button( icon='fa-info-circle',
                                tooltip='Photon energy maximum (in eV) when lindabs = .true.''\n''\n'
                                       'DEFAULT: 10', 
                                layout = Layout(width='40px'))

icon_omegastep = Button( icon='fa-info-circle',
                                tooltip='Steps in photon energy (in eV) when lindabs = .true.''\n''\n'
                                        'DEFAULT: 10', 
                                layout = Layout(width='40px'))

icon_n_r = Button( icon='fa-info-circle',
                                tooltip='index of refraction needed to convert the imaginary dielectric function (directly computed by EPW) to an absorption coefficient.''\n'
                                          'Look up from experimental measurement if unknown', 
                                layout = Layout(width='40px'))

icon_fsthick = Button( icon='fa-info-circle',
                                tooltip='Width of the Fermi surface window to take into account states in the self-energy delta functions in [eV].''\n'
                                          'Narrowing this value reduces the number of bands included in the selfenergy calculations.' '\n''\n'
                                          'DEFAULT: 1e10', 
                                layout = Layout(width='40px'))

icon_temps = Button( icon='fa-info-circle',
                                tooltip='Temperature values used in superconductivitiy, transport, indabs, etc. in kelvin unit. If no temps are provided, temps=300 and nstemp =1.''\n'
                                        'If two temps are provided, with temps(1)<temps(2) and nstemp >2, then temps is transformed into an evenly spaced grid with nstemp points, including temps(1) and temps(2) as the minimum and maximum values, respectively [Ex) nstemp      =   5 temps       = 300 500].''\n'
                                        'In this case, points are spaced according to (temps(2) - temps(1)) / (nstemp-1). Otherwise, temps is treated as a list, with the given temperatures used directly [Ex) temps    = 17 20 30].''\n'
                                        'No more than 50 temperatures can be supplied in this way.'
                                        'DEFAULT: 300', 
                                layout = Layout(width='40px'))

icon_degaussw = Button( icon='fa-info-circle',
                                tooltip='Smearing in the energy-conserving delta functions in [eV]''\n''\n'
                                           'DEFAULT: 0.025', 
                                layout = Layout(width='40px'))

icon_degaussq = Button( icon='fa-info-circle',
                                tooltip='Smearing for sum over q in the e-ph coupling in [meV]''\n''\n'
                                        'DEFAULT: 0.05', 
                                layout = Layout(width='40px'))

icon_efermi_read = Button( icon='fa-info-circle',
                                tooltip='If .true. the Fermi energy is read from the input file.''\n''\n'
                                          'DEFAULT: false', 
                                layout = Layout(width='40px'))

form_items = [
    
    
     HBox([Label(value='prgtkk'),
          Box([prgtkk, icon_prgtkk])
         ], layout = form_item_layout()),   
    
      HBox([Label(value='lindabs'),
          Box([lindabs, icon_lindabs])
         ], layout = form_item_layout()),
    
      HBox([Label(value='scissor'),
          Box([scissor, icon_scissor])
         ], layout = form_item_layout()),    
    
      HBox([Label(value='omegamin'),
          Box([omegamin, icon_omegamin])
         ], layout = form_item_layout()),      
    
       HBox([Label(value='omegamax'),
          Box([omegamax, icon_omegamax])
         ], layout = form_item_layout()),     
    
       HBox([Label(value='omegastep'),
          Box([omegastep, icon_omegastep])
         ], layout = form_item_layout()),      
    
          HBox([Label(value='n_r'),
          Box([n_r, icon_n_r])
         ], layout = form_item_layout()),   
    
       HBox([Label(value='fsthick'),
          Box([fsthick, icon_fsthick])
         ], layout = form_item_layout()),      

       HBox([Label(value='temps'),
          Box([temps, icon_temps])
         ], layout = form_item_layout()),   
    
       HBox([Label(value='degaussw'),
          Box([degaussw, icon_degaussw])
         ], layout = form_item_layout()),      
    
        HBox([Label(value='degaussq'),
          Box([degaussq, icon_degaussq])
         ], layout = form_item_layout()),    
    
     HBox([Label(value='efermi_read'),
          Box([efermi_read, icon_efermi_read])
         ], layout = form_item_layout()), 
    
    
        HBox([Label(value='Documentation URL:'), documentation_link], layout=Layout())
    
    #Box([Label(value='Print Electron-Photon Vertexs'), prgtkk], layout=form_item_layout()),
    #Box([Label(value='Calculate Optical Parameters'), lindabs], layout=form_item_layout()),
    #Box([Label(value='Scissor Shift'), scissor], layout=form_item_layout()),
    #Box([Label(value='Min Photon Energy'), omegamin], layout=form_item_layout()),
    #Box([Label(value='Max Photon Energy'), omegamax], layout=form_item_layout()),
    #Box([Label(value='Steps Photon Energy'), omegastep], layout=form_item_layout()),
    #Box([Label(value='Refractive Index'), n_r], layout=form_item_layout()),
    #Box([Label(value='Band Gap Width'), fsthick], layout=form_item_layout()),
    #Box([Label(value='Temperature (K)'), temps], layout=form_item_layout()),
    #Box([Label(value='Step Function Broadening Parameter'), degaussw], layout=form_item_layout()),
    #Box([Label(value='User Specific Fermi Energy'), degaussq], layout=form_item_layout()),
    #Box([Label(value='Fermi Energy'), efermi_read], layout=form_item_layout())
]

misc_box = Box(form_items, layout=box_layout(35))
misc_box


# %% [markdown]
# ## Mesh Sampling

# %%
## 

## 
# @var nkf1
# kpoints in 100 crystal direction, for sampling periodic cell in reciprocal space in the electronic structure part of the final wannier interpolation
nkf1 = IntText(name="nkf1",cvalue=20, layout=input_layout(30))       

## 
# @var nkf2
# kpoints in 010 crystal direction, for sampling periodic cell in reciprocal space in the electronic structure part of the final wannier interpolation
nkf2 = IntText(name="nkf2", value=20, layout=input_layout(30))  

## 
# @var nkf3
# kpoints in 001 crystal direction, for sampling periodic cell in reciprocal space in the electronic structure part of the final wannier interpolation
nkf3 = IntText(name="nkf3", value=20, layout=input_layout(30)) 

## 
# @var nqf1
# fine mesh kpoints in 100 crystal direction, for sampling periodic cell in reciprocal space in the electronic structure part of the final wannier interpolation
nqf1 = IntText(name="nqf1", value=20, layout=input_layout(30)) 

## 
# @var nqf2
# fine kpoints in 010 crystal direction, for sampling periodic cell in reciprocal space in the electronic structure part of the final wannier interpolation
nqf2 = IntText(name="nqf2", value=20, layout=input_layout(30)) 

## @var nqf3
# fine kpoints in 001 crystal direction, for sampling periodic cell in reciprocal space in the electronic structure part of the final wannier interpolation
nqf3 = IntText(name="nqf3", value=20, layout=input_layout(30)) 


## 
# @section tool_tip_button_creation
# Creates the info (hint) buttons and adds in the description when you hover your cursor over the button

icon_nkf1 = Button( icon='fa-info-circle',
                                tooltip='Dimensions of the fine electron grid''\n''\n'
                                           'DEFAULT: 0', 
                                layout = Layout(width='40px'))

icon_nkf2 = Button( icon='fa-info-circle',
                                tooltip='Dimensions of the fine electron grid''\n''\n'
                                           'DEFAULT: 0', 
                                layout = Layout(width='40px'))

icon_nkf3 = Button( icon='fa-info-circle',
                                tooltip='Dimensions of the fine electron grid''\n''\n'
                                           'DEFAULT: 0', 
                                layout = Layout(width='40px'))

icon_nqf1 = Button( icon='fa-info-circle',
                                tooltip='Dimensions of the fine phonon grid''\n''\n'
                                       'DEFAULT: 0', 
                                layout = Layout(width='40px'))

icon_nqf2 = Button( icon='fa-info-circle',
                                tooltip='Dimensions of the fine phonon grid''\n''\n'
                                       'DEFAULT: 0', 
                                layout = Layout(width='40px'))

icon_nqf3 = Button( icon='fa-info-circle',
                                tooltip='Dimensions of the fine phonon grid''\n''\n'
                                       'DEFAULT: 0', 
                                layout = Layout(width='40px'))

form_items = [
    
      HBox([Label(value='nkf1'),
          Box([nkf1, icon_nkf1])
         ], layout = form_item_layout()),   
    
       HBox([Label(value='nkf2'),
          Box([nkf2, icon_nkf2])
         ], layout = form_item_layout()),     
    
       HBox([Label(value='nkf3'),
          Box([nkf3, icon_nkf3])
         ], layout = form_item_layout()),  
    
        HBox([Label(value='nqf1'),
          Box([nqf1, icon_nqf1])
         ], layout = form_item_layout()), 
        
        HBox([Label(value='nqf2'),
          Box([nqf2, icon_nqf2])
         ], layout = form_item_layout()), 
    
        HBox([Label(value='nqf3'),
          Box([nqf3, icon_nqf3])
         ], layout = form_item_layout()), 
    
    #Box([Label(value='nkf1'), nkf1], layout=form_item_layout()),
    #Box([Label(value='nkf2'), nkf2], layout=form_item_layout()),
    #Box([Label(value='nkf3'), nkf3], layout=form_item_layout()),
    #Box([Label(value='nqf1'), nqf3], layout=form_item_layout()),
    #Box([Label(value='nqf2'), nqf2], layout=form_item_layout()),
    #Box([Label(value='nqf3'), nqf3], layout=form_item_layout()),
]

mesh_box = Box(form_items, layout=box_layout(20))

mesh_box


# %% [markdown]
# ## Combine UI Forms for this section

# %%
## combine forms into tabs
# @section combine_forms
# Combines the forms into tabs
# @var epw_x_tabs
# Tabs for the epw.x input file
# @var tab_contents
# Contents of the tabs
# @var children
# Children of the tabs

epw_x_tabs = Tab()
tab_contents = [setup_box, wannier_box, misc_box, mesh_box]
children = [content for content in tab_contents]
epw_x_tabs.children = children

epw_x_tabs.set_title(0, "Setup")
epw_x_tabs.set_title(1, "Wannier")
epw_x_tabs.set_title(2, "Misc")
epw_x_tabs.set_title(3, "Mesh Sampling")

epw_x_tabs

# %%
## NON USER INPUTS

# %%
## option for centering Wannier function at center of Wigner-Seitz cell
use_ws = '.false.'

w90_data = '''
wdata(1) = 'bands_plot = .true.'
wdata(2) = 'begin kpoint_path'
wdata(3) = 'L 0.50 0.00 0.00 G 0.00 0.00 0.00'
wdata(4) = 'G 0.00 0.00 0.00 X 0.50 0.50 0.00'
wdata(5) = 'end kpoint_path'
wdata(6) = 'bands_plot_format = gnuplot'
wdata(7) = 'guiding_centres = .true.'
wdata(8) = 'dis_num_iter = 500'
wdata(9) = 'num_print_cycles = 10'
wdata(10) = 'dis_mix_ratio = 1.0'
wdata(11) = 'use_ws_distance = T'
'''

epw_cores = 40
epw_walltime = '00:30:00'

# %% [markdown]
# ## Bind Inputs to Outputs

# %%
## binds EPW inputs to outputs
# @section bind_epw_inputs Inputs for EPW
def bind_EPW_X_inputs(self):
        
    epw_inputs = {
    'material_prefix': material_prefix.value,
    'atomic_mass': atomic_mass.value,
    'outdir': outdir.value,
    'kmaps': kmaps.value,
    'epbwrite': epbwrite.value,
    'epbread': epbread.value,
    'epwwrite': epwwrite.value,
    'epwread': epwread.value,
    'lpolar': lpolar.value,
    'use_ws': use_ws,
    'nbndsub': nbnd, 
    'wannierize': wannierize.value,
    'num_iter': num_iter.value,
    'iprint': iprint.value,
    'dis_win_max': dis_win_max.value,
    'dis_froz_max': dis_froz_max.value,
    'proj': projections.value,
    'wdata': w90_data,    
    'prtgkk': prgtkk.value,
    'lindabs': lindabs.value,
    'scissor': scissor.value,
    'omegamin': omegamin.value,                                                                                                                                                                         
    'omegamax': omegamax.value,
    'omegastep': omegastep.value,
    'n_r': n_r.value,  
    'fsthick': fsthick.value,
    'temps': temps.value,
    'degaussw': degaussw.value,
    'degaussq': degaussq.value,
    'efermi_read': efermi_read.value,
    'fermi_energy': "", #fermi_energy
    'nk1': kptx.value,
    'nk2': kpty.value,
    'nk3': kptz.value,
    'nq1': nq1.value,
    'nq2': nq2.value,
    'nq3': nq3.value,
    'nkf1': nkf1.value,
    'nkf2': nkf2.value,
    'nkf3': nkf3.value,
    'nqf1': nqf1.value,
    'nqf2': nqf2.value,
    'nqf3': nqf2.value,
    'qpoint_list': "" #qpoint_list_append
    }
    
    build_EPW_Input(epw_inputs, material_prefix.value)
        

# %%
## builds EPW input file
# @section build_epw_input Build EPW input
# @param epw_inputs, material_prefix
# 
def build_EPW_Input(epw_inputs, material_prefix):
    print('BUILDING EPW INPUT FILE')
    #Build epw input file
    
    ## @var epw_name
    epw_name = 'epw-%s.in' % material_prefix #assigns epw input file to variable epw_name 
    ## @var input_file
    input_file = '''
    --
    &inputepw
    prefix = '{material_prefix}'
    amass(1) = {atomic_mass}
    outdir = '{outdir}'
    iverbosity = 1
    dvscf_dir = './save'

    elph = .true.
    kmaps = {kmaps}
    epbwrite = {epbwrite}
    epbread = {epbread}
    epwwrite = {epwwrite}
    epwread = {epwread}
    etf_mem = 1

    lpolar = {lpolar}
    use_ws = {use_ws}

    nbndsub = {nbndsub}

    wannierize = {wannierize}
    num_iter = {num_iter}
    iprint = {iprint}
    dis_win_max = {dis_win_max}
    dis_froz_max= {dis_froz_max}
    {proj}
    {wdata}

    elecselfen = .false.
    phonselfen = .false.
    a2f = .false.
    prtgkk = {prtgkk}

    efermi_read = {efermi_read}
    fermi_energy = {fermi_energy}

    lindabs = {lindabs}
    scissor = {scissor}                                                                                                                                                                           
    omegamin = {omegamin}                                                                                                                                                                            
    omegamax = {omegamax}                                                                                                                                                                         
    omegastep= {omegastep}                                                                                                                                                                         
    n_r = {n_r}


    fsthick = {fsthick}
    temps = {temps}
    degaussw = {degaussw}
    degaussq = {degaussq}

    nk1 = {nk1}
    nk2 = {nk2}
    nk3 = {nk3}
    nq1 = {nq1}
    nq2 = {nq2}
    nq3 = {nq3}

    nkf1 = {nkf1}
    nkf2 = {nkf2}
    nkf3 = {nkf3}
    nqf1 = {nqf1}
    nqf2 = {nqf2}
    nqf3 = {nqf3}
    /
    {qpoint_list}
    '''.format(**epw_inputs) ## assigns information in ''' ''' to variable inputfile

    with open(epw_name, "w") as f: ## opens file epw_name
        f.write(input_file) ## writes inputfile to file epw_na


